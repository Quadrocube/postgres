/*-------------------------------------------------------------------------
 *
 * spgkdtreeproc.c
 *	  implementation of k-d tree over points for SP-GiST
 *
 *
 * Portions Copyright (c) 1996-2014, PostgreSQL Global Development Group
 * Portions Copyright (c) 1994, Regents of the University of California
 *
 * IDENTIFICATION
 *			src/backend/access/spgist/spgkdtreeproc.c
 *
 *-------------------------------------------------------------------------
 */

#include "postgres.h"

#include "access/gist.h"		/* for RTree strategy numbers */
#include "access/spgist.h"
#include "catalog/pg_type.h"
#include "utils/builtins.h"
#include "utils/geo_decls.h"
#include "access/spgist_proc.h"


Datum
spg_kd_config(PG_FUNCTION_ARGS)
{
	/* spgConfigIn *cfgin = (spgConfigIn *) PG_GETARG_POINTER(0); */
	spgConfigOut *cfg = (spgConfigOut *) PG_GETARG_POINTER(1);

	cfg->prefixType = FLOAT8OID;
	cfg->labelType = VOIDOID;	/* we don't need node labels */
	cfg->canReturnData = true;
	cfg->longValuesOK = false;
	cfg->suppLen = sizeof (BOX);
	PG_RETURN_VOID();
}

static int
getSide(double coord, bool isX, Point *tst)
{
	double		tstcoord = (isX) ? tst->x : tst->y;

	if (coord == tstcoord)
		return 0;
	else if (coord > tstcoord)
		return 1;
	else
		return -1;
}

Datum
spg_kd_choose(PG_FUNCTION_ARGS)
{
	spgChooseIn *in = (spgChooseIn *) PG_GETARG_POINTER(0);
	spgChooseOut *out = (spgChooseOut *) PG_GETARG_POINTER(1);
	Point	   *inPoint = DatumGetPointP(in->datum);
	double		coord;

	if (in->allTheSame)
		elog(ERROR, "allTheSame should not occur for k-d trees");

	Assert(in->hasPrefix);
	coord = DatumGetFloat8(in->prefixDatum);

	Assert(in->nNodes == 2);

	out->resultType = spgMatchNode;
	out->result.matchNode.nodeN =
		(getSide(coord, in->level % 2, inPoint) > 0) ? 0 : 1;
	out->result.matchNode.levelAdd = 1;
	out->result.matchNode.restDatum = PointPGetDatum(inPoint);

	PG_RETURN_VOID();
}

typedef struct SortedPoint
{
	Point	   *p;
	int			i;
} SortedPoint;

static int
x_cmp(const void *a, const void *b)
{
	SortedPoint *pa = (SortedPoint *) a;
	SortedPoint *pb = (SortedPoint *) b;

	if (pa->p->x == pb->p->x)
		return 0;
	return (pa->p->x > pb->p->x) ? 1 : -1;
}

static int
y_cmp(const void *a, const void *b)
{
	SortedPoint *pa = (SortedPoint *) a;
	SortedPoint *pb = (SortedPoint *) b;

	if (pa->p->y == pb->p->y)
		return 0;
	return (pa->p->y > pb->p->y) ? 1 : -1;
}


Datum
spg_kd_picksplit(PG_FUNCTION_ARGS)
{
	spgPickSplitIn *in = (spgPickSplitIn *) PG_GETARG_POINTER(0);
	spgPickSplitOut *out = (spgPickSplitOut *) PG_GETARG_POINTER(1);
	int			i;
	int			middle;
	SortedPoint *sorted;
	double		coord;

	sorted = palloc(sizeof(*sorted) * in->nTuples);
	for (i = 0; i < in->nTuples; i++)
	{
		sorted[i].p = DatumGetPointP(in->datums[i]);
		sorted[i].i = i;
	}

	qsort(sorted, in->nTuples, sizeof(*sorted),
		  (in->level % 2) ? x_cmp : y_cmp);
	middle = in->nTuples >> 1;
	coord = (in->level % 2) ? sorted[middle].p->x : sorted[middle].p->y;

	out->hasPrefix = true;
	out->prefixDatum = Float8GetDatum(coord);

	out->nNodes = 2;
	out->nodeLabels = NULL;		/* we don't need node labels */

	out->mapTuplesToNodes = palloc(sizeof(int) * in->nTuples);
	out->leafTupleDatums = palloc(sizeof(Datum) * in->nTuples);

	/*
	 * Note: points that have coordinates exactly equal to coord may get
	 * classified into either node, depending on where they happen to fall in
	 * the sorted list.  This is okay as long as the inner_consistent function
	 * descends into both sides for such cases.  This is better than the
	 * alternative of trying to have an exact boundary, because it keeps the
	 * tree balanced even when we have many instances of the same point value.
	 * So we should never trigger the allTheSame logic.
	 */
	for (i = 0; i < in->nTuples; i++)
	{
		Point	   *p = sorted[i].p;
		int			n = sorted[i].i;

		out->mapTuplesToNodes[n] = (i < middle) ? 0 : 1;
		out->leafTupleDatums[n] = PointPGetDatum(p);
	}

	PG_RETURN_VOID();
}

Datum
spg_kd_inner_consistent(PG_FUNCTION_ARGS)
{
	spgInnerConsistentIn *in = (spgInnerConsistentIn *) PG_GETARG_POINTER(0);
	spgInnerConsistentOut *out = (spgInnerConsistentOut *) PG_GETARG_POINTER(1);
	double		coord;
	int			which;
	int			i;
    Datum *boxes;

	Assert(in->hasPrefix);
	coord = DatumGetFloat8(in->prefixDatum);

	if (in->allTheSame)
		elog(ERROR, "allTheSame should not occur for k-d trees");

	Assert(in->nNodes == 2);

	if (in->norderbys > 0) {
		out->distances = palloc(in->nNodes * sizeof (double *));
	}
	
	/* "which" is a bitmask of children that satisfy all constraints */
	which = (1 << 1) | (1 << 2);

	for (i = 0; i < in->nkeys; i++)
	{
		Point	   *query = DatumGetPointP(in->scankeys[i].sk_argument);
		BOX		   *boxQuery;

		switch (in->scankeys[i].sk_strategy)
		{
			case RTLeftStrategyNumber:
				if ((in->level % 2) != 0 && FPlt(query->x, coord))
					which &= (1 << 1);
				break;
			case RTRightStrategyNumber:
				if ((in->level % 2) != 0 && FPgt(query->x, coord))
					which &= (1 << 2);
				break;
			case RTSameStrategyNumber:
				if ((in->level % 2) != 0)
				{
					if (FPlt(query->x, coord))
						which &= (1 << 1);
					else if (FPgt(query->x, coord))
						which &= (1 << 2);
				}
				else
				{
					if (FPlt(query->y, coord))
						which &= (1 << 1);
					else if (FPgt(query->y, coord))
						which &= (1 << 2);
				}
				break;
			case RTBelowStrategyNumber:
				if ((in->level % 2) == 0 && FPlt(query->y, coord))
					which &= (1 << 1);
				break;
			case RTAboveStrategyNumber:
				if ((in->level % 2) == 0 && FPgt(query->y, coord))
					which &= (1 << 2);
				break;
			case RTContainedByStrategyNumber:

				/*
				 * For this operator, the query is a box not a point.  We
				 * cheat to the extent of assuming that DatumGetPointP won't
				 * do anything that would be bad for a pointer-to-box.
				 */
				boxQuery = DatumGetBoxP(in->scankeys[i].sk_argument);

				if ((in->level % 2) != 0)
				{
					if (FPlt(boxQuery->high.x, coord))
						which &= (1 << 1);
					else if (FPgt(boxQuery->low.x, coord))
						which &= (1 << 2);
				}
				else
				{
					if (FPlt(boxQuery->high.y, coord))
						which &= (1 << 1);
					else if (FPgt(boxQuery->low.y, coord))
						which &= (1 << 2);
				}
				break;
			default:
				elog(ERROR, "unrecognized strategy number: %d",
					 in->scankeys[i].sk_strategy);
				break;
		}

		if (which == 0)
			break;				/* no need to consider remaining conditions */
	}

	/* We must descend into the children identified by which */
	out->nodeNumbers = (int *) palloc(sizeof(int) * 2);
	out->nNodes = 0;
	
	if (in->level == 0 && in->norderbys > 0) {
        BOX *newbox = palloc(sizeof(BOX));
        Point newp;
        newp.x = newp.y = get_float8_infinity();
        newbox->high = newp;
        newp.x = newp.y = -get_float8_infinity();
        newbox->low = newp;
        in->suppValue = PointerGetDatum(newbox);
    }
    if (DatumGetBoxP(in->suppValue) != NULL) {
		boxes = (Datum *) palloc(sizeof(Datum) * 2);
		out->suppValues = (Datum *) palloc(sizeof(Datum) * 2);
		BOX *area = DatumGetBoxP(in->suppValue);
		BOX *newbox1 = (BOX *) palloc0(sizeof(BOX));
		BOX *newbox2 = (BOX *) palloc0(sizeof(BOX));
		Point p1, p2;
		switch (in->level % 2) {
            case 0:
				p1.x = p2.x = coord;
				p1.y = area->high.y;
				p2.y = area->low.y;
				break;
			case 1:
				p1.y = p2.y = coord;
				p1.x = area->high.x;
				p2.x = area->low.x;
				break;
		}
		newbox1->low = area->low;
		newbox1->high = p1;
		boxes[0] = BoxPGetDatum(newbox1);
		newbox2->low = p2;
		newbox2->high = area->high;
		boxes[1] = BoxPGetDatum(newbox2);
	}
	
	for (i = 1; i <= 2; i++)
	{
		if (which & (1 << i)) {
			out->nodeNumbers[out->nNodes++] = i - 1;
            if (in->suppValue != 0) {
                out->suppValues[out->nNodes-1] = boxes[i];
            }
			if (in->norderbys > 0) {
				spg_point_distance(out->suppValues[out->nNodes-1],
					in->norderbys, in->orderbyKeys, &out->distances[out->nNodes-1], false);
			}
		}
	}

	/* Set up level increments, too */
	out->levelAdds = (int *) palloc(sizeof(int) * 2);
	out->levelAdds[0] = 1;
	out->levelAdds[1] = 1;

	PG_RETURN_VOID();
}

/*
 * spg_kd_leaf_consistent() is the same as spg_quad_leaf_consistent(),
 * since we support the same operators and the same leaf data type.
 * So we just borrow that function.
 */
