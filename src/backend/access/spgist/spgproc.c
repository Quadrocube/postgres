#include "access/spgist_proc.h"

int
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

void
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

RBNode *
SpGistSearchTreeItemAllocator(void *arg)
{
	IndexScanDesc scan = (IndexScanDesc) arg;

	return (RBNode *) palloc(SPGISTHDRSZ + sizeof(double) * scan->numberOfOrderBys);
}

void
SpGistSearchTreeItemDeleter(RBNode *rb, void *arg)
{
	pfree(rb);
}

/* 
 * Construct SpGistSearchTreeItem storing item and add it to queue 
 * 
 * Called in queue context 
 */
void 
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

/*
 * Leaf SpGistSearchItem constructor, called in queue context
 */
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
	newItem->suppValue = (Datum) 0;
	return newItem;
}

void
freeSearchTreeItem(SpGistScanOpaque so, SpGistSearchItem *item)
{
	if (so->state.config.suppLen > 0
            && DatumGetPointer(item->suppValue) != NULL
            && item->itemState == INNER) {
		pfree(DatumGetPointer(item->suppValue));
    }
    if (!so->state.attType.attbyval &&
            DatumGetPointer(item->value) != NULL) {
        pfree(DatumGetPointer(item->value));
    }

	pfree(item);
}

/* Point-box distance in the assumption that box is aligned by axis */
double dist_pb_simplified(Datum p, Datum b) {
	Point *point = DatumGetPointP(p);
	BOX *box = DatumGetBoxP(b);
	double dx = 0.0, dy = 0.0;

	if (point->x < box->low.x)
		dx = box->low.x - point->x;
	if (point->x > box->high.x)
		dx = point->x - box->high.x;
	if (point->y < box->low.y)
		dy = box->low.y - point->y;
	if (point->y > box->high.y)
		dy = point->y - box->high.y;
	return HYPOT(dx,dy);
}

void
spg_point_distance(Datum to, int norderbys, 
        ScanKey orderbyKeys, double **distances, bool isLeaf) 
{
    int sk_num;
    double *distance;
    *distances = malloc(norderbys * sizeof (double *));
    *distance = *distances;
    for (sk_num = 0; sk_num < norderbys; ++sk_num) {
        Datum from_point = orderbyKeys[sk_num].sk_argument;
        if (isLeaf) {
            *distance = DatumGetFloat8 (
                    DirectFunctionCall2(point_distance, from_point, to) );
        } else {
            *distance = dist_pb_simplified(from_point, to);
        }
        distance++;
    }
}
