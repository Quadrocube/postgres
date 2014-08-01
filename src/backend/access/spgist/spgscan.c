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
#include "utils/debug_spg.h"

extern double get_float8_infinity();

typedef void (*storeRes_func) (SpGistScanOpaque so, ItemPointerData heapPtr,
								 Datum leafValue, bool isnull, bool recheck);

static void
freeSearchTreeItem(SpGistScanOpaque so, SpGistSearchItem *item)
{
	if (so->state.config.suppLen > 0
            && DatumGetPointer(item->value) != NULL
            && item->itemState == INNER) {
		pfree(DatumGetPointer(item->value));
    }
	pfree(item);
}

/*
 * Initialize scanStack to search the root page, resetting
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
		/* Stack a work item to scan the null index entries */
		startEntry = (SpGistSearchItem *) palloc0(sizeof(SpGistSearchItem));
		ItemPointerSet(&startEntry->heap, SPGIST_NULL_BLKNO, FirstOffsetNumber);
		startEntry->itemState = INNER;
		startEntry->level = 0;
		addSearchItemToQueue(scan, startEntry, so->distances);
	}

	if (so->searchNonNulls)
	{
		/* Stack a work item to scan the non-null index entries */
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
	
	/* set up starting stack entries */
	resetSpGistScanOpaque(scan);
	
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
	bool		result;
	Datum		leafDatum;
	spgLeafConsistentIn in;
	spgLeafConsistentOut out;
	FmgrInfo   *procinfo;
	MemoryContext oldCtx = CurrentMemoryContext;
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
		goto report;
	}

	leafDatum = SGLTDATUM(leafTuple, &so->state);

	/* use temp context for calling leaf_consistent */
	oldCtx = MemoryContextSwitchTo(so->tempCxt);

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
		if (scan->numberOfOrderBys > 0) {
			MemoryContextSwitchTo(so->queueCxt);
			addSearchItemToQueue(scan,
				newHeapItem(so, level, leafTuple->heapPtr, leafValue, recheck), 
				out.distances);
			MemoryContextSwitchTo(oldCtx);
		} else {
			MemoryContextSwitchTo(oldCtx);
			storeRes(so, leafTuple->heapPtr, leafValue, isnull, recheck);
			*reportedSome = true;
		}
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
 *
 * -- searchForKNearest - set in case of performing KNN-Search, denotes the parameter K
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
			/* We store heap items in the queue only in case of ordered search */
			Assert(scan->numberOfOrderBys > 0);
			/* Heap items can only be stored on leaf pages */
			Assert(SpGistPageIsLeaf(page));
			storeRes(so, item->heap, item->value, isnull, 
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
				/* When root is a leaf, examine all its tuples */
				for (offset = FirstOffsetNumber; offset <= max; offset++)
				{
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
				/* Normal case: just examine the chain we arrived at */
				while (offset != InvalidOffsetNumber)
				{
					Assert(offset >= FirstOffsetNumber && offset <= max);
					leafTuple = (SpGistLeafTuple)
						PageGetItem(page, PageGetItemId(page, offset));
					if (leafTuple->tupstate != SPGIST_LIVE)
					{
						if (leafTuple->tupstate == SPGIST_REDIRECT)
						{
							/* redirection tuple should be first in chain */
							Assert(offset == ItemPointerGetOffsetNumber(&stackEntry->ptr));
							/* transfer attention to redirect point */
							item->heap = ((SpGistDeadTuple) leafTuple)->pointer;
							Assert(ItemPointerGetBlockNumber(&stackEntry->ptr) != SPGIST_METAPAGE_BLKNO);
							goto redirect;
						}
						if (leafTuple->tupstate == SPGIST_DEAD)
						{
							/* dead tuple should be first in chain */
							Assert(offset == ItemPointerGetOffsetNumber(&stackEntry->ptr));
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
					offset = leafTuple->nextOffset;
				}
			}
		}
		else	/* page is inner */
		{
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
					Assert(ItemPointerGetBlockNumber(&stackEntry->ptr) != SPGIST_METAPAGE_BLKNO);
					goto redirect;
				}
				elog(ERROR, "unexpected SPGiST tuple state: %d",
					 innerTuple->tupstate);
			}

			/* use temp context for calling inner_consistent */
			oldCxt = MemoryContextSwitchTo(so->tempCxt);
			
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

			for (i = 0; i < out.nNodes; i++)
			{
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
						newItem->level = item->level + out.levelAdds[nodeN];
					else
						newItem->level = item->level;
					/* Must copy value out of temp context */
					if (out.reconstructedValues) {
						newItem->value =
							datumCopy(out.reconstructedValues[nodeN],
									  so->state.attType.attbyval,
									  so->state.attType.attlen);
					} else {
						newItem->value = (Datum) 0;
					}
					
					if (out.suppValues) {
						newItem->suppValue = 
							datumCopy(out.suppValues[nodeN], false,
									  so->state.config.suppLen);
					} else {
						newItem->suppValue = (Datum) 0;
					}
					
					if (out.distances != NULL) {
						distances = out.distances[nodeN];
					} else {
						distances = inf_distances;
					}
					newItem->itemState = INNER;
					addSearchItemToQueue(scan, newItem, distances);
				}
			}
			MemoryContextSwitchTo(oldCxt);
		}

		/* done with this scan entry */
		oldCxt = MemoryContextSwitchTo(so->queueCxt);
		freeSearchTreeItem(so, item);
		MemoryContextSwitchTo(oldCxt);
		/* clear temp context before proceeding to the next one */
		MemoryContextReset(so->tempCxt);
	}

	if (buffer != InvalidBuffer)
		UnlockReleaseBuffer(buffer);
}

/* storeRes subroutine for getbitmap case */
static void
storeBitmap(SpGistScanOpaque so, ItemPointerData heapPtr,
			Datum leafValue, bool isnull, bool recheck)
{
	tbm_add_tuples(so->tbm, &heapPtr, 1, recheck);
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
storeGettuple(SpGistScanOpaque so, ItemPointerData heapPtr,
			  Datum leafValue, bool isnull, bool recheck)
{
	Assert(so->nPtrs < MaxIndexTuplesPerPage);
	so->heapPtrs[so->nPtrs] = heapPtr;
	so->recheck[so->nPtrs] = recheck;
	if (so->want_itup)
	{
		/*
		 * Reconstruct desired IndexTuple.  We have to copy the datum out of
		 * the temp context anyway, so we may as well create the tuple here.
		 */
		so->indexTups[so->nPtrs] = index_form_tuple(so->indexTupDesc,
													&leafValue,
													&isnull);
	}
	so->nPtrs++;
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
		if (so->iPtr < so->nPtrs)
		{
			/* continuing to return reported tuples */
			scan->xs_ctup.t_self = so->heapPtrs[so->iPtr];
			scan->xs_recheck = so->recheck[so->iPtr];
			scan->xs_itup = so->indexTups[so->iPtr];
			so->iPtr++;
			PG_RETURN_BOOL(true);
		}

		if (so->want_itup)
		{
			/* Must pfree IndexTuples to avoid memory leak */
			int			i;

			for (i = 0; i < so->nPtrs; i++)
				pfree(so->indexTups[i]);
		}
		so->iPtr = so->nPtrs = 0;

		spgWalk(scan->indexRelation, scan, false, storeGettuple);

		if (so->nPtrs == 0)
			break;				/* must have completed scan */
	}

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
