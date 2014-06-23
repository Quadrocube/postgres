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

#include "access/relscan.h"
#include "access/spgist_private.h"
#include "miscadmin.h"
#include "storage/bufmgr.h"
#include "utils/datum.h"
#include "utils/memutils.h"
#include "utils/rel.h"


typedef void (*storeRes_func) (SpGistScanOpaque so, ItemPointer heapPtr,
								 Datum leafValue, bool isnull, bool recheck);

typedef struct ScanStackEntry
{
	Datum		reconstructedValue;		/* value reconstructed from parent */
	int			level;			/* level of items on this page */
	ItemPointerData ptr;		/* block and offset to scan from */
} ScanStackEntry;


/* Free a ScanStackEntry */
static void
freeScanStackEntry(SpGistScanOpaque so, ScanStackEntry *stackEntry)
{
	if (!so->state.attType.attbyval &&
		DatumGetPointer(stackEntry->reconstructedValue) != NULL)
		pfree(DatumGetPointer(stackEntry->reconstructedValue));
	pfree(stackEntry);
}

/* Free the entire stack */
static void
freeScanStack(SpGistScanOpaque so)
{
	ListCell   *lc;

	foreach(lc, so->scanStack)
	{
		freeScanStackEntry(so, (ScanStackEntry *) lfirst(lc));
	}
	list_free(so->scanStack);
	so->scanStack = NIL;
}

/*
 * Initialize scanStack to search the root page, resetting
 * any previously active scan
 */
static void
resetSpGistScanOpaque(SpGistScanOpaque so)
{
	ScanStackEntry *startEntry;

	freeScanStack(so);

	if (so->searchNulls)
	{
		/* Stack a work item to scan the null index entries */
		startEntry = (ScanStackEntry *) palloc0(sizeof(ScanStackEntry));
		ItemPointerSet(&startEntry->ptr, SPGIST_NULL_BLKNO, FirstOffsetNumber);
		so->scanStack = lappend(so->scanStack, startEntry);
	}

	if (so->searchNonNulls)
	{
		/* Stack a work item to scan the non-null index entries */
		startEntry = (ScanStackEntry *) palloc0(sizeof(ScanStackEntry));
		ItemPointerSet(&startEntry->ptr, SPGIST_ROOT_BLKNO, FirstOffsetNumber);
		so->scanStack = lappend(so->scanStack, startEntry);
	}

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
	int norderbys = PG_GETARG_INT32(2);
	int knearest = PG_GETARG_INT32(3);

	IndexScanDesc scan;
	SpGistScanOpaque so;

	scan = RelationGetIndexScan(rel, keysz, 0);

	so = (SpGistScanOpaque) palloc0(sizeof(SpGistScanOpaqueData));
	if (keysz > 0)
		so->keyData = (ScanKey) palloc(sizeof(ScanKeyData) * keysz);
	else
		so->keyData = NULL;
	initSpGistState(&so->state, scan->indexRelation);
	scan->numberOfOrderBys = norderbys;
	so->nnCount = knearest;
	so->tempCxt = AllocSetContextCreate(CurrentMemoryContext,
										"SP-GiST search temporary context",
										ALLOCSET_DEFAULT_MINSIZE,
										ALLOCSET_DEFAULT_INITSIZE,
										ALLOCSET_DEFAULT_MAXSIZE);

	/* Set up indexTupDesc and xs_itupdesc in case it's an index-only scan */
	so->indexTupDesc = scan->xs_itupdesc = RelationGetDescr(rel);

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
	ScanKey		scankey = (ScanKey) PG_GETARG_POINTER(1);
	int norderbys = PG_GETARG_INT32(2);
	int knearest = PG_GETARG_INT32(3);

	/* copy scankeys into local storage */
	if (scankey && scan->numberOfKeys > 0)
	{
		memmove(scan->keyData, scankey,
				scan->numberOfKeys * sizeof(ScanKeyData));
	}

	/* preprocess scankeys, set up the representation in *so */
	spgPrepareScanKeys(scan);

	/* set up starting stack entries */
	resetSpGistScanOpaque(so);

	/* TODO: Complete initialization of order distances */
	scan->numberOfOrderBys = norderbys;
	so->nnCount = knearest;

	MemoryContext oldCxt;
	oldCxt = MemoryContextSwitchTo(so->queueCxt);
	so->queue = rb_create(SPGISTHDRSZ + sizeof (double) * scan->numberOfOrderBys,
			SpGistSearchTreeItemComparator,
			SpGistSearchTreeItemCombiner,
			SpGistSearchTreeItemAllocator,
			SpGistSearchTreeItemDeleter,
			scan);
	MemoryContextSwitchTo(oldCxt);

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
			Datum *leafValue, bool *recheck)
{
	bool		result;
	Datum		leafDatum;
	spgLeafConsistentIn in;
	spgLeafConsistentOut out;
	FmgrInfo   *procinfo;
	MemoryContext oldCtx;
	SpGistScanOpaque so = scan->opaque;

	if (isnull)
	{
		/* Should not have arrived on a nulls page unless nulls are wanted */
		Assert(so->searchNulls);
		*leafValue = (Datum) 0;
		*recheck = false;
		return true;
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

	out.leafValue = (Datum) 0;
	out.recheck = false;

	procinfo = index_getprocinfo(index, 1, SPGIST_LEAF_CONSISTENT_PROC);
	result = DatumGetBool(FunctionCall2Coll(procinfo,
											index->rd_indcollation[0],
											PointerGetDatum(&in),
											PointerGetDatum(&out)));

	*leafValue = out.leafValue;
	*recheck = out.recheck;

	if (!result) return false;

	/* TODO: Fix distance computation */
	ScanKey key = scan->orderByData;
	int keySize = scan->numberOfOrderBys;
	double *distance_p = so->distances;
	while (keySize > 0) {
		Datum datum;
		bool isNull;

		datum = index_getattr(tuple,
				key->sk_attno,
				giststate->tupdesc,
				&isNull);

		if ((key->sk_flags & SK_ISNULL) || isNull) {
			/* Assume distance computes as null and sorts to the end */
			*distance_p = get_float8_infinity();
		} else {
			Datum dist;
			GISTENTRY de;

			gistdentryinit(giststate, key->sk_attno - 1, &de,
					datum, r, page, offset,
					FALSE, isNull);

			/*
			 * Call the Distance function to evaluate the distance.  The
			 * arguments are the index datum (as a GISTENTRY*), the comparison
			 * datum, and the ordering operator's strategy number and subtype
			 * from pg_amop.
			 */
			dist = FunctionCall4Coll(&key->sk_func,
					key->sk_collation,
					PointerGetDatum(&de),
					key->sk_argument,
					Int32GetDatum(key->sk_strategy),
					ObjectIdGetDatum(key->sk_subtype));

			*distance_p = DatumGetFloat8(dist);
		}

		key++;
		distance_p++;
		keySize--;
	}

	MemoryContextSwitchTo(oldCtx);

	return result;
}

typedef struct SpGistSearchItem {
	struct SpGistSearchItem *next; /* list link */
	ItemPointerData heapPtr;
	Datum leafValue;

} SpGistSearchItem;

typedef struct SpGistSearchTreeItem {
	SpGistSearchItem *head; /* first chain member */
	double distances[1];
} SpGistSearchTreeItem;

#define SPGISTHDRSZ offsetof(GISTSearchTreeItem, distances)

static int
SpGistSearchTreeItemComparator(const RBNode *a, const RBNode *b, void *arg) {
	const SpGistSearchTreeItem *sa = (const SpGistSearchTreeItem *) a;
	const SpGistSearchTreeItem *sb = (const SpGistSearchTreeItem *) b;
	IndexScanDesc scan = (IndexScanDesc) arg;
	int i;

	/* Order according to distance comparison so that the farthest element is the leftmost */
	for (i = 0; i < scan->numberOfOrderBys; i++) {
		if (sa->distances[i] != sb->distances[i])
			return (sa->distances[i] > sb->distances[i]) ? -1 : 1;
	}

	return 0;
}

static void
SpGistSearchTreeItemCombiner(RBNode *existing, const RBNode *newrb, void *arg) {
	SpGistSearchTreeItem *scurrent = (SpGistSearchTreeItem *) existing;
	const SpGistSearchTreeItem *snew = (const SpGistSearchTreeItem *) newrb;
	SpGistSearchItem *newitem = snew->head;

	/* snew should have just one item in its chain */
	Assert(newitem && newitem->next == NULL);

	newitem->next = scurrent->head;
	scurrent->head = newitem;
}

static RBNode *
SpGistSearchTreeItemAllocator(void *arg) {
	IndexScanDesc scan = (IndexScanDesc) arg;

	return palloc(SPGISTHDRSZ + sizeof (double) * scan->numberOfOrderBys);
}

static void
SpGistSearchTreeItemDeleter(RBNode *rb, void *arg) {
	pfree(rb);
}

/*
 * Create new item to insert in rb-tree
 * Invoked in heap context
 */
static void makeTreeItem(ItemPointerData heapPtr, Datum leafValue,
		IndexScanDesc scan, SpGistSearchTreeItem *tmpItem) {
	SpGistSearchItem *item;
	SpGistScanOpaque so = (SpGistScanOpaque) scan->opaque;
	item = palloc(sizeof (SpGistSearchItem));
	item->next = NULL;
	item->heapPtr = heapPtr;
	item->leafValue = leafValue;

	tmpItem->head = item;
	memcpy(tmpItem->distances, so->distances,
			sizeof (double) * scan->numberOfOrderBys);
}

static void updateQueue(IndexScanDesc scan, ItemPointer heapPtr,
		int *foundNNearest, Datum leafValue, bool isnull) {
	if (isnull) return;
	bool isNew;
	SpGistScanOpaque so = (SpGistScanOpaque) scan->opaque;
	if (*foundNNearest < so->nnCount) {
		/* We search for K nearest and still have not filled the K-bucket up
			 => Add the item to the nearestQueue */
		makeTreeItem(heapPtr, leafValue, scan, so->tmpTreeItem);
		(void) rb_insert(so->queue, (RBNode *) so->tmpTreeItem, &isNew);
		++(*foundNNearest);
	}
	else {
		/* We already have seen at least K items
		   If current one is better -> remove the worst and add current to nearestQueue */
		SpGistSearchTreeItem currentWorst = (SpGistSearchTreeItem) rb_leftmost(so->queue);
		for (int orderby = 0; orderby < scan->numberOfOrderBys; ++orderby) {
			if (so->distances[orderby] < currentWorst.distances[orderby]) {
				SpGistSearchItem listHead = currentWorst->head;
				if (listHead->next == NULL) {
					(void) rb_delete(so->queue, (RBNode) currentWorst);
				} else {
					currentWorst->head = listHead->next;
				}
				makeTreeItem(heapPtr, leafValue, scan, so->tmpTreeItem);
				(void) rb_insert(so->queue, (RBNode *) so->tmpTreeItem, &isNew);
				break;
			}
		}
	}
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
	MemoryContext oldcxt;
	SpGistScanOpaque so = (SpGistScanOpaque) scan->opaque;
	RBTree *nearestQueue = so->queue;
	const int searchForKNearest = (const int) so->nnCount;
	int foundNNearest = 0;

	while (scanWholeIndex || !reportedSome || searchForKNearest) {
		ScanStackEntry *stackEntry;
		BlockNumber blkno;
		OffsetNumber offset;
		Page		page;
		bool		isnull;

		/* Pull next to-do item from the list */
		if (so->scanStack == NIL)
			break;				/* there are no more pages to scan */

		stackEntry = (ScanStackEntry *) linitial(so->scanStack);
		so->scanStack = list_delete_first(so->scanStack);

redirect:
		/* Check for interrupts, just in case of infinite loop */
		CHECK_FOR_INTERRUPTS();

		blkno = ItemPointerGetBlockNumber(&stackEntry->ptr);
		offset = ItemPointerGetOffsetNumber(&stackEntry->ptr);

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

		if (SpGistPageIsLeaf(page))
		{
			SpGistLeafTuple leafTuple;
			OffsetNumber max = PageGetMaxOffsetNumber(page);
			Datum		leafValue = (Datum) 0;
			bool		recheck = false;

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
					if (spgLeafTest(index, so,
									leafTuple, isnull,
									stackEntry->level,
									stackEntry->reconstructedValue,
									&leafValue,
									&recheck))
					{
						if (searchForKNearest) {
							oldcxt = MemoryContextSwitchTo(so->queueCxt);
							updateQueue(so, &leafTuple->heapPtr,
									leafValue, &foundNNearest, isnull);
							MemoryContextSwitchTo(oldcxt);
						} else {
						storeRes(so, &leafTuple->heapPtr,
								 leafValue, isnull, recheck);
						}
						reportedSome = true;
					}
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
							stackEntry->ptr = ((SpGistDeadTuple) leafTuple)->pointer;
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
					if (spgLeafTest(index, so,
									leafTuple, isnull,
									stackEntry->level,
									stackEntry->reconstructedValue,
									&leafValue,
									&recheck))
					{
						if (searchForKNearest) {
							oldcxt = MemoryContextSwitchTo(so->queueCxt);
							updateQueue(so, &leafTuple->heapPtr,
									leafValue, &foundNNearest, isnull);
							MemoryContextSwitchTo(oldcxt);
						} else {
						storeRes(so, &leafTuple->heapPtr,
								 leafValue, isnull, recheck);
						}
						reportedSome = true;
					}

					offset = leafTuple->nextOffset;
				}
			}
		}
		else	/* page is inner */
		{
			SpGistInnerTuple innerTuple;
			spgInnerConsistentIn in;
			spgInnerConsistentOut out;
			FmgrInfo   *procinfo;
			SpGistNodeTuple *nodes;
			SpGistNodeTuple node;
			int			i;
			MemoryContext oldCtx;

			innerTuple = (SpGistInnerTuple) PageGetItem(page,
												PageGetItemId(page, offset));

			if (innerTuple->tupstate != SPGIST_LIVE)
			{
				if (innerTuple->tupstate == SPGIST_REDIRECT)
				{
					/* transfer attention to redirect point */
					stackEntry->ptr = ((SpGistDeadTuple) innerTuple)->pointer;
					Assert(ItemPointerGetBlockNumber(&stackEntry->ptr) != SPGIST_METAPAGE_BLKNO);
					goto redirect;
				}
				elog(ERROR, "unexpected SPGiST tuple state: %d",
					 innerTuple->tupstate);
			}

			if (searchForKNearest && foundNNearest == searchForKNearest) {
				oldcxt = MemoryContextSwitchTo(so->queueCxt);
				SpGistSearchTreeItem currentWorst = (SpGistSearchTreeItem) rb_leftmost(nearestQueue);
				if (so->areaFilter(currentWorst->distances)) {
					freeScanStackEntry(so, stackEntry);
					continue;
				}
				MemoryContextSwitchTo(oldcxt);
			}

			/* use temp context for calling inner_consistent */
			oldCtx = MemoryContextSwitchTo(so->tempCxt);

			in.scankeys = so->keyData;
			in.nkeys = so->numberOfKeys;
			in.reconstructedValue = stackEntry->reconstructedValue;
			in.level = stackEntry->level;
			in.returnData = so->want_itup;
			in.allTheSame = innerTuple->allTheSame;
			in.hasPrefix = (innerTuple->prefixSize > 0);
			in.prefixDatum = SGITDATUM(innerTuple, &so->state);
			in.nNodes = innerTuple->nNodes;
			in.nodeLabels = spgExtractNodeLabels(&so->state, innerTuple);

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
				procinfo = index_getprocinfo(index, 1, SPGIST_INNER_CONSISTENT_PROC);
				FunctionCall2Coll(procinfo,
								  index->rd_indcollation[0],
								  PointerGetDatum(&in),
								  PointerGetDatum(&out));
			}
			else
			{
				/* force all children to be visited */
				out.nNodes = in.nNodes;
				out.nodeNumbers = (int *) palloc(sizeof(int) * in.nNodes);
				for (i = 0; i < in.nNodes; i++)
					out.nodeNumbers[i] = i;
			}

			MemoryContextSwitchTo(oldCtx);

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
					ScanStackEntry *newEntry;

					/* Create new work item for this node */
					newEntry = palloc(sizeof(ScanStackEntry));
					newEntry->ptr = nodes[nodeN]->t_tid;
					if (out.levelAdds)
						newEntry->level = stackEntry->level + out.levelAdds[i];
					else
						newEntry->level = stackEntry->level;
					/* Must copy value out of temp context */
					if (out.reconstructedValues)
						newEntry->reconstructedValue =
							datumCopy(out.reconstructedValues[i],
									  so->state.attType.attbyval,
									  so->state.attType.attlen);
					else
						newEntry->reconstructedValue = (Datum) 0;

					so->scanStack = lcons(newEntry, so->scanStack);
				}
			}
		}

		/* done with this scan stack entry */
		freeScanStackEntry(so, stackEntry);
		/* clear temp context before proceeding to the next one */
		MemoryContextReset(so->tempCxt);
	}

	if (searchForKNearest) {
		oldcxt = MemoryContextSwitchTo(so->queueCxt);
		rb_begin_iterate(nearestQueue, DirectWalk);
		while (foundNNearest > 0) {
			SpGistSearchTreeItem titem = (SpGistSearchTreeItem) rb_iterate(nearestQueue);
			SpGistSearchItem item = titem.head;
			while (item.next != NULL) {
				storeRes(so, &item.heapPtr, item.leafValue, false, false);
				foundNNearest--;
				item = item.next;
			}
			storeRes(so, &item.heapPtr, item.leafValue, false, false);
			foundNNearest--;
		}
		MemoryContextSwitchTo(oldcxt);
	}

	if (buffer != InvalidBuffer)
		UnlockReleaseBuffer(buffer);
}

/* storeRes subroutine for getbitmap case */
static void
storeBitmap(SpGistScanOpaque so, ItemPointer heapPtr,
			Datum leafValue, bool isnull, bool recheck)
{
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
	Assert(so->nPtrs < MaxIndexTuplesPerPage);
	so->heapPtrs[so->nPtrs] = *heapPtr;
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
			/* continuing to return tuples from a leaf page */
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
