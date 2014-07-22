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

enum SPGistSEARCHITEMSTATE {
    HEAP_RECHECK,        /* SearchItem is heap item and rechek is needed before reporting */
    HEAP_NORECHECK,      /* SearchItem is heap item and no rechek is needed */
    INNER               /* SearchItem is inner tree item - rechek irrelevant */
};

typedef struct SpGistSearchItem
{
        SPGistSEARCHITEMSTATE itemState;        /* see above */
    	Datum   value;                          /* value reconstructed from parent or leafValue if heaptuple */
	int	level;                          /* level of items on this page */
	struct SpGistSearchItem *next;          /* list link */
        ItemPointerData heap;                   /* heap info, if heap tuple */
} GISTSearchItem;

typedef struct SpGistSearchTreeItem
{
	RBNode		rbnode;			/* this is an RBTree item */
	SpGistSearchItem *head;		/* first chain member */
	SpGistSearchItem *lastHeap;	/* last heap-tuple member, if any */
	double		distances[1];	/* array with numberOfOrderBys entries */
} GISTSearchTreeItem;

#define SPGISTHDRSZ offsetof(GISTSearchTreeItem, distances)
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
	GISTSearchTreeItem *scurrent = (GISTSearchTreeItem *) existing;
	const GISTSearchTreeItem *snew = (const GISTSearchTreeItem *) newrb;
	GISTSearchItem *newitem = snew->head;

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

#define GSTIHDRSZ offsetof(SpGistearchTreeItem, distances)

static RBNode *
SpGistSearchTreeItemAllocator(void *arg)
{
	IndexScanDesc scan = (IndexScanDesc) arg;

	return palloc(GSTIHDRSZ + sizeof(double) * scan->numberOfOrderBys);
}

static void
SpGistSearchTreeItemDeleter(RBNode *rb, void *arg)
{
	pfree(rb);
}

/*
 * Create new item to insert in rb-tree
 * Invoked in heap context
 */
static void initTreeItem(const SpGistSearchItem *source, const double *distances,
        SpGistSearchTreeItem *dest) {
	source->next = NULL;
	dest->head = source;
	dest->distances = distances == NULL ? dest->distances : distances;
}

static void 
addSearchItemToQueue(IndexScanDesc *scan, SpGistSearchItem *item, double *distances) {
	bool isNew;
	SpGistScanOpaque so = (SpGistScanOpaque) scan->opaque;
	SpGistSearchTreeItem *newItem = so->tmpTreeItem;
	memcpy(newItem->distances, distances, scan->numberOfOrderBys * sizeof (double));
	newItem->head = item;
	item->next = NULL;
	newItem->lastHeap = SPGISTSearchItemIsHeap(*item) ? item : NULL;
	rb_insert(so->queue, (RBNode *) newItem, &isNew);
}

GISTSearchItem *newHeapItem(int level, ItemPointerData heapPtr, 
		Datum leafValue, bool recheck) {
	GISTSearchItem *newItem = (GISTSearchItem *) palloc(sizeof(GISTSearchItem));
	newItem->next = NULL;
	newItem->level = level;
	newItem->heap = heapPtr;
	newItem->value = leafValue;
	newItem->itemState = recheck ? HEAP_RECHECK : HEAP_NORECHECK;
	
}

#ifdef	__cplusplus
}
#endif

#endif	/* SPGIST_SEARCH_H */

