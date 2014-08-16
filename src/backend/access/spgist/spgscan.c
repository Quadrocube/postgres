/*-------------------------------------------------------------------------
 *
 * spgscan.c
 *	  routines for scanning SP-GiST indexes
 *
 *
 * Portions Copyright (c) 1996-2014, PostgreSQL Global Development Group
 * Portions Copyright (c) 1994, Regents of the University of California
 *
 * IDENTIFICATION
 *			src/backend/access/spgist/spgscan.c
 *
 *-------------------------------------------------------------------------
 */

#include "postgres.h"

#include "access/spgist_private.h"
#include "access/relscan.h"
#include "miscadmin.h"
#include "storage/bufmgr.h"
#include "utils/datum.h"
#include "utils/memutils.h"
#include "utils/rel.h"
#include "access/spgist_proc.h"

extern double get_float8_infinity();

typedef void (*storeRes_func) (SpGistScanOpaque so, ItemPointer heapPtr,
								 Datum leafValue, bool isnull, bool recheck);

static void
freeSearchTreeItem(SpGistScanOpaque so, SpGistSearchItem *item)
{
    elog(WARNING, "suppLen == %d", so->state.config.suppLen);
	if (so->state.config.suppLen > 0
            && DatumGetPointer(item->suppValue) != NULL
            && item->itemState == INNER) {
        elog(WARNING, "pfree == %d", item->suppValue);
		pfree(DatumGetPointer(item->suppValue));
    }
    if (!so->state.attType.attbyval &&
            DatumGetPointer(item->value) != NULL) {
        elog(WARNING, "pfree == %d", item->value);
        pfree(DatumGetPointer(item->value));
    }

    elog(WARNING, "pfree == %d", item);
	pfree(item);
}

/*
 * Initialize queue to search the root page, resetting
 * any previously active scan
 */
static void
resetSpGistScanOpaque(IndexScanDesc scan)
{
	SpGistScanOpaque so = (SpGistScanOpaque) scan->opaque;
	SpGistSearchItem *startEntry;

	MemoryContext oldCtx;
	oldCtx = MemoryContextSwitchTo(so->queueCxt);
	memset(so->distances, 0, sizeof(double) * scan->numberOfOrderBys);
	

	if (so->searchNulls)
	{
		/* Add a work item to scan the null index entries */
		startEntry = (SpGistSearchItem *) palloc0(sizeof(SpGistSearchItem));
		ItemPointerSet(&startEntry->heap, SPGIST_NULL_BLKNO, FirstOffsetNumber);
		startEntry->itemState = INNER;
		startEntry->level = 0;
		addSearchItemToQueue(scan, startEntry, so->distances);
	}

	if (so->searchNonNulls)
	{
		/* Add a work item to scan the non-null index entries */
		startEntry = (SpGistSearchItem *) palloc0(sizeof(SpGistSearchItem));
		ItemPointerSet(&startEntry->heap, SPGIST_ROOT_BLKNO, FirstOffsetNumber);
		startEntry->itemState = INNER;
		startEntry->level = 0;
		addSearchItemToQueue(scan, startEntry, so->distances);
	}
	
	MemoryContextSwitchTo(oldCtx);

	if (so->want_itup)
	{
		/* Must pfree IndexTuples to avoid memory leak */
		int			i;

		for (i = 0; i < so->nPtrs; i++)
			pfree(so->indexTups[i]);
	}
	so->iPtr = so->nPtrs = 0;
}

/*
 * Prepare scan keys in SpGistScanOpaque from caller-given scan keys
 *
 * Sets searchNulls, searchNonNulls, numberOfKeys, keyData fields of *so.
 *
 * The point here is to eliminate null-related considerations from what the
 * opclass consistent functions need to deal with.  We assume all SPGiST-
 * indexable operators are strict, so any null RHS value makes the scan
 * condition unsatisfiable.  We also pull out any IS NULL/IS NOT NULL
 * conditions; their effect is reflected into searchNulls/searchNonNulls.
 */
static void
spgPrepareScanKeys(IndexScanDesc scan)
{
	SpGistScanOpaque so = (SpGistScanOpaque) scan->opaque;
	bool		qual_ok;
	bool		haveIsNull;
	bool		haveNotNull;
	int			nkeys;
	int			i;

	if (scan->numberOfKeys <= 0)
	{
		/* If no quals, whole-index scan is required */
		so->searchNulls = true;
		so->searchNonNulls = true;
		so->numberOfKeys = 0;
		return;
	}

	/* Examine the given quals */
	qual_ok = true;
	haveIsNull = haveNotNull = false;
	nkeys = 0;
	for (i = 0; i < scan->numberOfKeys; i++)
	{
		ScanKey		skey = &scan->keyData[i];

		if (skey->sk_flags & SK_SEARCHNULL)
			haveIsNull = true;
		else if (skey->sk_flags & SK_SEARCHNOTNULL)
			haveNotNull = true;
		else if (skey->sk_flags & SK_ISNULL)
		{
			/* ordinary qual with null argument - unsatisfiable */
			qual_ok = false;
			break;
		}
		else
		{
			/* ordinary qual, propagate into so->keyData */
			so->keyData[nkeys++] = *skey;
			/* this effectively creates a not-null requirement */
			haveNotNull = true;
		}
	}

	/* IS NULL in combination with something else is unsatisfiable */
	if (haveIsNull && haveNotNull)
		qual_ok = false;

	/* Emit results */
	if (qual_ok)
	{
		so->searchNulls = haveIsNull;
		so->searchNonNulls = haveNotNull;
		so->numberOfKeys = nkeys;
	}
	else
	{
		so->searchNulls = false;
		so->searchNonNulls = false;
		so->numberOfKeys = 0;
	}
}

Datum
spgbeginscan(PG_FUNCTION_ARGS)
{
    elog(WARNING, "beginscan fired");
	Relation	rel = (Relation) PG_GETARG_POINTER(0);
	int			keysz = PG_GETARG_INT32(1);
	int			orderbys = PG_GETARG_INT32(2);

	IndexScanDesc scan;
	SpGistScanOpaque so;

	scan = RelationGetIndexScan(rel, keysz, orderbys);

	so = (SpGistScanOpaque) palloc0(sizeof(SpGistScanOpaqueData));
	if (keysz > 0)
		so->keyData = (ScanKey) palloc(sizeof(ScanKeyData) * keysz);
	else
		so->keyData = NULL;
	initSpGistState(&so->state, scan->indexRelation);
	
	so->tempCxt = AllocSetContextCreate(CurrentMemoryContext,
										"SP-GiST search temporary context",
										ALLOCSET_DEFAULT_MINSIZE,
										ALLOCSET_DEFAULT_INITSIZE,
										ALLOCSET_DEFAULT_MAXSIZE);

	/* Set up indexTupDesc and xs_itupdesc in case it's an index-only scan */
	so->indexTupDesc = scan->xs_itupdesc = RelationGetDescr(rel);
	so->tmpTreeItem = palloc(SPGISTHDRSZ + sizeof(double) * scan->numberOfOrderBys);
	so->distances = palloc(sizeof(double) * scan->numberOfOrderBys);

	scan->opaque = so;

	so->queueCxt = AllocSetContextCreate(CurrentMemoryContext,
			"SP-GiST queue context",
			ALLOCSET_DEFAULT_MINSIZE,
			ALLOCSET_DEFAULT_INITSIZE,
			ALLOCSET_DEFAULT_MAXSIZE);

	PG_RETURN_POINTER(scan);
}

Datum
spgrescan(PG_FUNCTION_ARGS)
{
    elog(WARNING, "rescan fired");
	IndexScanDesc scan = (IndexScanDesc) PG_GETARG_POINTER(0);
	SpGistScanOpaque so = (SpGistScanOpaque) scan->opaque;
	ScanKey		scankeys = (ScanKey) PG_GETARG_POINTER(1);
	ScanKey		orderbys = (ScanKey) PG_GETARG_POINTER(3);
	
	MemoryContext oldCxt;

	/* copy scankeys into local storage */
	if (scankeys && scan->numberOfKeys > 0)
	{
		memmove(scan->keyData, scankeys,
				scan->numberOfKeys * sizeof(ScanKeyData));
	}

	/* preprocess scankeys, set up the representation in *so */
	spgPrepareScanKeys(scan);

	MemoryContextReset(so->queueCxt);
	oldCxt = MemoryContextSwitchTo(so->queueCxt);
	so->queue = rb_create(SPGISTHDRSZ + sizeof (double) * scan->numberOfOrderBys,
			SpGistSearchTreeItemComparator,
			SpGistSearchTreeItemCombiner,
			SpGistSearchTreeItemAllocator,
			SpGistSearchTreeItemDeleter,
			scan);
	MemoryContextSwitchTo(oldCxt);
	
	/* set up starting queue entries */
	resetSpGistScanOpaque(scan);
    so->curTreeItem = NULL;
	
	if (orderbys && scan->numberOfOrderBys > 0)
	{
		memmove(scan->orderByData, orderbys,
				scan->numberOfOrderBys * sizeof(ScanKeyData));
	}

	PG_RETURN_VOID();
}

Datum
spgendscan(PG_FUNCTION_ARGS)
{
	IndexScanDesc scan = (IndexScanDesc) PG_GETARG_POINTER(0);
	SpGistScanOpaque so = (SpGistScanOpaque) scan->opaque;

	MemoryContextDelete(so->tempCxt);

	PG_RETURN_VOID();
}

Datum
spgmarkpos(PG_FUNCTION_ARGS)
{
	elog(ERROR, "SPGiST does not support mark/restore");
	PG_RETURN_VOID();
}

Datum
spgrestrpos(PG_FUNCTION_ARGS)
{
	elog(ERROR, "SPGiST does not support mark/restore");
	PG_RETURN_VOID();
}

/*
 * Test whether a leaf tuple satisfies all the scan keys
 *
 * *leafValue is set to the reconstructed datum, if provided
 * *recheck is set true if any of the operators are lossy
 *
 * If scan keys are satisfied, fill up so->distances in accordance to so->
 */
static bool
spgLeafTest(Relation index, IndexScanDesc scan,
			SpGistLeafTuple leafTuple, bool isnull,
			int level, Datum reconstructedValue,
			bool *reportedSome, storeRes_func storeRes)
{
    elog(WARNING, "Leaf test: entering");
	bool		result;
	Datum		leafDatum;
	spgLeafConsistentIn in;
	spgLeafConsistentOut out;
	FmgrInfo   *procinfo;
	MemoryContext oldCtx;
	SpGistScanOpaque so = scan->opaque;
	Datum leafValue;
	bool recheck;

	if (isnull)
	{
		/* Should not have arrived on a nulls page unless nulls are wanted */
		Assert(so->searchNulls);
		leafValue = (Datum) 0;
		recheck = false;
        result = true;
	} else {
        leafDatum = SGLTDATUM(leafTuple, &so->state);
    }

	/* use temp context for calling leaf_consistent */
	oldCtx = MemoryContextSwitchTo(so->tempCxt);
    if (isnull) goto report;

	in.scankeys = so->keyData;
	in.nkeys = so->numberOfKeys;
	in.reconstructedValue = reconstructedValue;
	in.level = level;
	in.returnData = so->want_itup;
	in.leafDatum = leafDatum;
	in.orderbykeys = scan->orderByData;
	in.norderbys = scan->numberOfOrderBys;

	out.leafValue = (Datum) 0;
	out.recheck = false;

	procinfo = index_getprocinfo(index, 1, SPGIST_LEAF_CONSISTENT_PROC);
	result = DatumGetBool(FunctionCall2Coll(procinfo,
											index->rd_indcollation[0],
											PointerGetDatum(&in),
											PointerGetDatum(&out)));
	recheck = out.recheck;
	leafValue = out.leafValue;
	
report:
	if (result) {
        elog(WARNING, "Leaf test: success, storing");
		if (scan->numberOfOrderBys > 0) {
			MemoryContextSwitchTo(so->queueCxt);
			addSearchItemToQueue(scan,
				newHeapItem(so, level, leafTuple->heapPtr, leafValue, recheck), 
				out.distances);
			MemoryContextSwitchTo(oldCtx);
		} else {
			MemoryContextSwitchTo(oldCtx);
			storeRes(so, &leafTuple->heapPtr, leafValue, isnull, recheck);
			*reportedSome = true;
		}
        elog(WARNING, "Leaf test: success, stored");
	} else {
        elog(WARNING, "check failed -> no store =(");
    }

	return result;
}

void inner_consistent_input_init(spgInnerConsistentIn *in, IndexScanDesc scan, 
			SpGistSearchItem *item, SpGistInnerTuple innerTuple) {
	SpGistScanOpaque so = scan->opaque;
	in->scankeys = so->keyData;
	in->nkeys = so->numberOfKeys;
	in->reconstructedValue = item->value;
	in->level = item->level;
	in->returnData = so->want_itup;
	in->allTheSame = innerTuple->allTheSame;
	in->hasPrefix = (innerTuple->prefixSize > 0);
	in->prefixDatum = SGITDATUM(innerTuple, &so->state);
	in->nNodes = innerTuple->nNodes;
	in->nodeLabels = spgExtractNodeLabels(&so->state, innerTuple);
	in->norderbys = scan->numberOfOrderBys;
    in->orderbyKeys = scan->orderByData;
	in->suppValue = item->suppValue;
}

/*
 * Walk the tree and report all tuples passing the scan quals to the storeRes
 * subroutine.
 *
 * If scanWholeIndex is true, we'll do just that.  If not, we'll stop at the
 * next page boundary once we have reported at least one tuple.
 */
static void
spgWalk(Relation index, IndexScanDesc scan, bool scanWholeIndex,
		storeRes_func storeRes)
{
	Buffer		buffer = InvalidBuffer;
	bool		reportedSome = false;
	MemoryContext oldCxt;
	SpGistScanOpaque so = (SpGistScanOpaque) scan->opaque;
	
	while (scanWholeIndex || !reportedSome) {
        elog(WARNING, "SpgWalk: newturn");
		BlockNumber blkno;
		OffsetNumber offset;
		Page		page;
		bool		isnull;

		SpGistSearchItem *item;

		oldCxt = MemoryContextSwitchTo(so->queueCxt);
init:
		/* Update curTreeItem if we don't have one */
		if (so->curTreeItem == NULL)
		{
			so->curTreeItem = (SpGistSearchTreeItem *) rb_leftmost(so->queue);
			/* Done when tree is empty */
			if (so->curTreeItem == NULL)
				return;
		}

		item = so->curTreeItem->head;
		if (item != NULL)
		{
			/* Delink item from chain */
			so->curTreeItem->head = item->next;
			if (item == so->curTreeItem->lastHeap)
				so->curTreeItem->lastHeap = NULL;
		} else {
		    /* curTreeItem is exhausted, so remove it from rbtree */
			rb_delete(so->queue, (RBNode *) so->curTreeItem);
			so->curTreeItem = NULL;
			goto init;
		}
		

		MemoryContextSwitchTo(oldCxt);
redirect:
		/* Check for interrupts, just in case of infinite loop */
		CHECK_FOR_INTERRUPTS();

		blkno = ItemPointerGetBlockNumber(&item->heap);
		offset = ItemPointerGetOffsetNumber(&item->heap);

		if (buffer == InvalidBuffer)
		{
			buffer = ReadBuffer(index, blkno);
			LockBuffer(buffer, BUFFER_LOCK_SHARE);
		}
		else if (blkno != BufferGetBlockNumber(buffer))
		{
			UnlockReleaseBuffer(buffer);
			buffer = ReadBuffer(index, blkno);
			LockBuffer(buffer, BUFFER_LOCK_SHARE);
		}
		
		/* else new pointer points to the same page, no work needed */

		page = BufferGetPage(buffer);

		isnull = SpGistPageStoresNulls(page) ? true : false;

		if (SPGISTSearchItemIsHeap(*item)) {
            elog(WARNING, "SpgWalk: Got heap item, reporting");
			/* We store heap items in the queue only in case of ordered search */
			Assert(scan->numberOfOrderBys > 0);
			/* Heap items can only be stored on leaf pages */
			Assert(SpGistPageIsLeaf(page));
			storeRes(so, &item->heap, item->value, isnull, 
				item->itemState == HEAP_RECHECK);
			reportedSome = true;
		}
		else if (SpGistPageIsLeaf(page))
		{
			/* Page is a leaf - that is, all it's tuples are heap items */
			SpGistLeafTuple leafTuple;
			OffsetNumber max = PageGetMaxOffsetNumber(page);

			if (SpGistBlockIsRoot(blkno))
			{
                elog(WARNING, "SpgWalk: Got leaf root");
				/* When root is a leaf, examine all its tuples */
				for (offset = FirstOffsetNumber; offset <= max; offset++)
				{
                    elog(WARNING, "SpgWalk: leafroot: new item");
					leafTuple = (SpGistLeafTuple)
						PageGetItem(page, PageGetItemId(page, offset));
					if (leafTuple->tupstate != SPGIST_LIVE)
					{
						/* all tuples on root should be live */
						elog(ERROR, "unexpected SPGiST tuple state: %d",
							 leafTuple->tupstate);
					}

					Assert(ItemPointerIsValid(&leafTuple->heapPtr));
					spgLeafTest(index, scan, leafTuple, isnull, item->level,
									item->value, &reportedSome, storeRes);
				}
			}
			else
			{
                elog(WARNING, "SpgWalk: Got leaf page");
				/* Normal case: just examine the chain we arrived at */
				while (offset != InvalidOffsetNumber)
				{
                    elog(WARNING, "SpgWalk: leafpage: new item");
					Assert(offset >= FirstOffsetNumber && offset <= max);
					leafTuple = (SpGistLeafTuple)
						PageGetItem(page, PageGetItemId(page, offset));
					if (leafTuple->tupstate != SPGIST_LIVE)
					{
						if (leafTuple->tupstate == SPGIST_REDIRECT)
						{
                            elog(WARNING, "SpgWalk: leafpage: redirect");
							/* redirection tuple should be first in chain */
							Assert(offset == ItemPointerGetOffsetNumber(&item->heap));
							/* transfer attention to redirect point */
							item->heap = ((SpGistDeadTuple) leafTuple)->pointer;
							Assert(ItemPointerGetBlockNumber(&item->heap) != SPGIST_METAPAGE_BLKNO);
							goto redirect;
						}
						if (leafTuple->tupstate == SPGIST_DEAD)
						{
                            elog(WARNING, "SpgWalk: leafpage: dead");
							/* dead tuple should be first in chain */
							Assert(offset == ItemPointerGetOffsetNumber(&item->heap));
							/* No live entries on this page */
							Assert(leafTuple->nextOffset == InvalidOffsetNumber);
							break;
						}
						/* We should not arrive at a placeholder */
						elog(ERROR, "unexpected SPGiST tuple state: %d",
							 leafTuple->tupstate);
					}

					Assert(ItemPointerIsValid(&leafTuple->heapPtr));
					spgLeafTest(index, scan, leafTuple, isnull, item->level,
									item->value, &reportedSome, storeRes);
                    elog(WARNING, "SpgWalk: leafpage: leaftuple == %d", leafTuple);
                    elog(WARNING, "SpgWalk: leafpage: nextoffset == %d", leafTuple->nextOffset);
					offset = leafTuple->nextOffset;
				}
			}
		}
		else	/* page is inner */
		{
            elog(WARNING, "SpgWalk: Got inner page");
			SpGistInnerTuple innerTuple;
			spgInnerConsistentIn in;
			spgInnerConsistentOut out;
			FmgrInfo   *consistent_procinfo;
			SpGistNodeTuple *nodes;
			SpGistNodeTuple node;
			int	i;

			innerTuple = (SpGistInnerTuple) PageGetItem(page,
												PageGetItemId(page, offset));

			if (innerTuple->tupstate != SPGIST_LIVE)
			{
				if (innerTuple->tupstate == SPGIST_REDIRECT)
				{
					/* transfer attention to redirect point */
					item->heap = ((SpGistDeadTuple) innerTuple)->pointer;
					Assert(ItemPointerGetBlockNumber(&item->heap) != SPGIST_METAPAGE_BLKNO);
					goto redirect;
				}
				elog(ERROR, "unexpected SPGiST tuple state: %d",
					 innerTuple->tupstate);
			}

			/* use temp context for calling inner_consistent */
			oldCxt = MemoryContextSwitchTo(so->tempCxt);
			
            elog(WARNING, "SpgWalk: inner_consistent_input_init");
			inner_consistent_input_init(&in, scan, item, innerTuple);
			double *inf_distances = palloc(scan->numberOfOrderBys * sizeof (double));
			for (i = 0; i < scan->numberOfOrderBys; ++i)
				inf_distances[i] = get_float8_infinity();

			/* collect node pointers */
			nodes = (SpGistNodeTuple *) palloc(sizeof(SpGistNodeTuple) * in.nNodes);
			SGITITERATE(innerTuple, i, node)
			{
				nodes[i] = node;
			}

			memset(&out, 0, sizeof(out));

            elog(WARNING, "SpgWalk: call consistent method");
			if (!isnull)
			{
				/* use user-defined inner consistent method */
				consistent_procinfo = index_getprocinfo(index, 1, SPGIST_INNER_CONSISTENT_PROC);
				FunctionCall2Coll(consistent_procinfo,
								  index->rd_indcollation[0],
								  PointerGetDatum(&in),
								  PointerGetDatum(&out));
			}
			else
			{
				/* force all children to be visited */
				out.nNodes = in.nNodes;
				out.nodeNumbers = (int *) palloc(sizeof(int) * in.nNodes);
				for (i = 0; i < in.nNodes; i++) {
					out.nodeNumbers[i] = i;
				}
			}

			MemoryContextSwitchTo(so->queueCxt);

			/* If allTheSame, they should all or none of 'em match */
			if (innerTuple->allTheSame)
				if (out.nNodes != 0 && out.nNodes != in.nNodes)
					elog(ERROR, "inconsistent inner_consistent results for allTheSame inner tuple");

            elog(WARNING, "SpgWalk: get result from consistent");
			for (i = 0; i < out.nNodes; i++)
			{
                elog(WARNING, "SpgWalk: get result from consistent for i == %d", i);
				int			nodeN = out.nodeNumbers[i];

				Assert(nodeN >= 0 && nodeN < in.nNodes);
				if (ItemPointerIsValid(&nodes[nodeN]->t_tid))
				{
					double *distances;
					SpGistSearchItem *newItem;

					/* Create new work item for this node */
					newItem = palloc(sizeof(SpGistSearchItem));
					newItem->heap = nodes[nodeN]->t_tid;
					if (out.levelAdds)
						newItem->level = item->level + out.levelAdds[i];
					else
						newItem->level = item->level;
					/* Must copy value out of temp context */
                    elog(WARNING, "SpgWalk: copy reconvalues, i == %d", i);
					if (out.reconstructedValues) {
						newItem->value =
							datumCopy(out.reconstructedValues[i],
									  so->state.attType.attbyval,
									  so->state.attType.attlen);
					} else {
						newItem->value = (Datum) 0;
					}
					
                    elog(WARNING, "SpgWalk: copy suppvalues, i == %d", i);
					if (out.suppValues) {
						newItem->suppValue = 
							datumCopy(out.suppValues[nodeN], false,
									  so->state.config.suppLen);
					} else {
						newItem->suppValue = (Datum) 0;
					}
					
                    elog(WARNING, "SpgWalk: copy distances, i == %d", i);
					if (out.distances != NULL) {
						distances = out.distances[nodeN];
					} else {
						distances = inf_distances;
					}
					newItem->itemState = INNER;
                    elog(WARNING, "SpgWalk: addSearchItem to queue, i == %d", i);
					addSearchItemToQueue(scan, newItem, distances);
				}
			}
			MemoryContextSwitchTo(oldCxt);
		}

        elog(WARNING, "SpgWalk: freeing queue entry");
		/* done with this scan entry */
		oldCxt = MemoryContextSwitchTo(so->queueCxt);
		freeSearchTreeItem(so, item);
		MemoryContextSwitchTo(oldCxt);
		/* clear temp context before proceeding to the next one */
		MemoryContextReset(so->tempCxt);
        elog(WARNING, "SpgWalk: endrun");
	}

	if (buffer != InvalidBuffer)
		UnlockReleaseBuffer(buffer);
    elog(WARNING, "SpgWalk: exit");
}

/* storeRes subroutine for getbitmap case */
static void
storeBitmap(SpGistScanOpaque so, ItemPointer heapPtr,
			Datum leafValue, bool isnull, bool recheck)
{
    elog(WARNING, "storeBitmap: entring");
	tbm_add_tuples(so->tbm, heapPtr, 1, recheck);
	so->ntids++;
}

Datum
spggetbitmap(PG_FUNCTION_ARGS)
{
	IndexScanDesc scan = (IndexScanDesc) PG_GETARG_POINTER(0);
	TIDBitmap  *tbm = (TIDBitmap *) PG_GETARG_POINTER(1);
	SpGistScanOpaque so = (SpGistScanOpaque) scan->opaque;

	/* Copy want_itup to *so so we don't need to pass it around separately */
	so->want_itup = false;

	so->tbm = tbm;
	so->ntids = 0;

	spgWalk(scan->indexRelation, scan, true, storeBitmap);

	PG_RETURN_INT64(so->ntids);
}

/* storeRes subroutine for gettuple case */
static void
storeGettuple(SpGistScanOpaque so, ItemPointer heapPtr,
			  Datum leafValue, bool isnull, bool recheck)
{
    elog(WARNING, "storeGettuple: entering");
	Assert(so->nPtrs < MaxIndexTuplesPerPage);
	so->heapPtrs[so->nPtrs] = *heapPtr;
	so->recheck[so->nPtrs] = recheck;
    elog(WARNING, "storeGettuple: recheck == %d", recheck);
	if (so->want_itup)
	{
		/*
		 * Reconstruct desired IndexTuple.  We have to copy the datum out of
		 * the temp context anyway, so we may as well create the tuple here.
		 */
        elog(WARNING, "storeGettuple: reconstructing itup, mc == %d, isnull == %d", CurrentMemoryContext, isnull);
		so->indexTups[so->nPtrs] = index_form_tuple(so->indexTupDesc,
													&leafValue,
													&isnull);
	}
	so->nPtrs++;
    elog(WARNING, "storeGettuple: leaving");
}

Datum
spggettuple(PG_FUNCTION_ARGS)
{
	IndexScanDesc scan = (IndexScanDesc) PG_GETARG_POINTER(0);
	ScanDirection dir = (ScanDirection) PG_GETARG_INT32(1);
	SpGistScanOpaque so = (SpGistScanOpaque) scan->opaque;

	if (dir != ForwardScanDirection)
		elog(ERROR, "SP-GiST only supports forward scan direction");

	/* Copy want_itup to *so so we don't need to pass it around separately */
	so->want_itup = scan->xs_want_itup;

	for (;;)
	{
        elog(WARNING, "spggettuple: newturn");
		if (so->iPtr < so->nPtrs)
		{
			/* continuing to return reported tuples */
			scan->xs_ctup.t_self = so->heapPtrs[so->iPtr];
			scan->xs_recheck = so->recheck[so->iPtr];
			scan->xs_itup = so->indexTups[so->iPtr];
			so->iPtr++;
            elog(WARNING, "spggettuple: returning true");
			PG_RETURN_BOOL(true);
		}

		if (so->want_itup)
		{
			/* Must pfree IndexTuples to avoid memory leak */
			int			i;

            elog(WARNING, "spggettuple: freeing index tuples, cmc == %d", CurrentMemoryContext);
			for (i = 0; i < so->nPtrs; i++) {
                elog(WARNING, "pfree ptr[i] == %d", so->indexTups[i]);
				pfree(so->indexTups[i]);
            }
            elog(WARNING, "spggettuple: completed freeing index tuples");
		}
		so->iPtr = so->nPtrs = 0;

		spgWalk(scan->indexRelation, scan, false, storeGettuple);

        elog(WARNING, "spggettuple: endturn");
		if (so->nPtrs == 0)
			break;				/* must have completed scan */
	}

    elog(WARNING, "spggettuple: returning false");
	PG_RETURN_BOOL(false);
}

Datum
spgcanreturn(PG_FUNCTION_ARGS)
{
	Relation	index = (Relation) PG_GETARG_POINTER(0);
	SpGistCache *cache;

	/* We can do it if the opclass config function says so */
	cache = spgGetCache(index);

	PG_RETURN_BOOL(cache->config.canReturnData);
}
