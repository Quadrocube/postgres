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

extern double get_float8_infinity();

extern int
SpGistSearchTreeItemComparator(const RBNode *a, const RBNode *b, void *arg);

extern void
SpGistSearchTreeItemCombiner(RBNode *existing, const RBNode *newrb, void *arg);

#define GSTIHDRSZ offsetof(SpGistSearchTreeItem, distances)

extern RBNode *
SpGistSearchTreeItemAllocator(void *arg);

extern void
SpGistSearchTreeItemDeleter(RBNode *rb, void *arg);

extern void 
addSearchItemToQueue(IndexScanDesc scan, SpGistSearchItem *item, double *distances);

extern SpGistSearchItem *newHeapItem(SpGistScanOpaque so, int level, 
        ItemPointerData heapPtr, Datum leafValue, bool recheck);

extern void
spg_point_distance(Datum to, int norderbys, 
        ScanKey orderbyKeys, double **distances, bool isLeaf); 

extern void
freeSearchTreeItem(SpGistScanOpaque so, SpGistSearchItem *item);

extern double pb_dist_simplified(Datum p, Datum b);


#ifdef	__cplusplus
}
#endif

#endif	/* SPGIST_SEARCH_H */

