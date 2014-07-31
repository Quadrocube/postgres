#ifndef SPGIST_SEARCH_H
#define	SPGIST_SEARCH_H
#include "postgres.h"

#include "access/relscan.h"
#include "access/spgist_private.h"
#include "miscadmin.h"
#include "storage/bufmgr.h"
#include "utils/datum.h"
#include "utils/memutils.h"
#include "utils/rel.h"

#include "utils/geo_decls.h"

#define SPGISTHDRSZ offsetof(SpGistSearchTreeItem, distances)
#define SPGISTSearchItemIsHeap(item)	((item).itemState == HEAP_RECHECK \
                                      || (item).itemState == HEAP_NORECHECK)

static int
SpGistSearchTreeItemComparator(const RBNode *a, const RBNode *b, void *arg)
{
	const SpGistSearchTreeItem *sa = (const SpGistSearchTreeItem *) a;
	const SpGistSearchTreeItem *sb = (const SpGistSearchTreeItem *) b;
	IndexScanDesc scan = (IndexScanDesc) arg;
	int			i;

	/* Order according to distance comparison */
	for (i = 0; i < scan->numberOfOrderBys; i++)
	{
		if (sa->distances[i] != sb->distances[i])
			return (sa->distances[i] > sb->distances[i]) ? 1 : -1;
	}

	return 0;
}

static void
SpGistSearchTreeItemCombiner(RBNode *existing, const RBNode *newrb, void *arg)
{
	SpGistSearchTreeItem *scurrent = (SpGistSearchTreeItem *) existing;
	const SpGistSearchTreeItem *snew = (const SpGistSearchTreeItem *) newrb;
	SpGistSearchItem *newitem = snew->head;

	/* snew should have just one item in its chain */
	Assert(newitem && newitem->next == NULL);

	/*
	 * If new item is heap tuple, it goes to front of chain; otherwise insert
	 * it before the first index-page item, so that index pages are visited in
	 * LIFO order, ensuring depth-first search of index pages.  See comments
	 * in gist_private.h.
	 */
	if (SPGISTSearchItemIsHeap(*newitem))
	{
		newitem->next = scurrent->head;
		scurrent->head = newitem;
		if (scurrent->lastHeap == NULL)
			scurrent->lastHeap = newitem;
	}
	else if (scurrent->lastHeap == NULL)
	{
		newitem->next = scurrent->head;
		scurrent->head = newitem;
	}
	else
	{
		newitem->next = scurrent->lastHeap->next;
		scurrent->lastHeap->next = newitem;
	}
}

#define GSTIHDRSZ offsetof(SpGistSearchTreeItem, distances)

static RBNode *
SpGistSearchTreeItemAllocator(void *arg)
{
	IndexScanDesc scan = (IndexScanDesc) arg;

	return (RBNode *) palloc(GSTIHDRSZ + sizeof(double) * scan->numberOfOrderBys);
}

static void
SpGistSearchTreeItemDeleter(RBNode *rb, void *arg)
{
	pfree(rb);
}

static void 
addSearchItemToQueue(IndexScanDesc scan, SpGistSearchItem *item, double *distances) {
	bool isNew;
	SpGistScanOpaque so = (SpGistScanOpaque) scan->opaque;
	SpGistSearchTreeItem *newItem = so->tmpTreeItem;
	memcpy(newItem->distances, distances, scan->numberOfOrderBys * sizeof (double));
	newItem->head = item;
	item->next = NULL;
	newItem->lastHeap = SPGISTSearchItemIsHeap(*item) ? item : NULL;
	rb_insert(so->queue, (RBNode *) newItem, &isNew);
}

SpGistSearchItem *newHeapItem(SpGistScanOpaque so, int level, 
        ItemPointerData heapPtr, Datum leafValue, bool recheck) {
	SpGistSearchItem *newItem = (SpGistSearchItem *) palloc(sizeof(SpGistSearchItem));
	newItem->next = NULL;
	newItem->level = level;
	newItem->heap = heapPtr;
        /* copy value to queue cxt out of tmp cxt */
        newItem->value = datumCopy(leafValue, so->state.attType.attbyval, 
                so->state.attType.attlen);
	newItem->itemState = recheck ? HEAP_RECHECK : HEAP_NORECHECK;
	return newItem;
}

void
spg_point_distance(Datum to, int norderbys, 
        ScanKey orderbyKeys, double **distances, bool isLeaf) 
{
	int sk_num;
        *distances = malloc(norderbys * sizeof (double *));
        double *distance = *distances;
        for (sk_num = 0; sk_num < norderbys; ++sk_num) {
            Datum from_point = orderbyKeys[sk_num].sk_argument;
            if (isLeaf) {
                *distance = DatumGetFloat8 (
                    DirectFunctionCall2(point_distance, from_point, to) );
            } else {
                *distance = DatumGetFloat8 (
                    DirectFunctionCall2(dist_pb, from_point, to) );
            }
            distance++;
        }
}


#ifdef	__cplusplus
}
#endif

#endif	/* SPGIST_SEARCH_H */

