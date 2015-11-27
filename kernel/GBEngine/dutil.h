/**
  UTILITIES FOR DYNAMIC ALGORITHM
*/

#ifndef __DUTIL_H
#define __DUTIL_H

#include <kernel/GBEngine/kutil.h>
#include <omalloc/omalloc.h>
#include <kernel/polys.h>
#include <kernel/ideals.h>
#include <kernel/GBEngine/kstd1.h>
#include <kernel/GBEngine/khstd.h>
#include <polys/kbuckets.h>
#include <polys/prCopy.h>
#include <polys/weight.h>
#include <misc/intvec.h>
#include "dynamic_engine.h"
#include "skeleton.h"

/**
  Creates a new ring from oldR, using new_weight to define the ordering.
  This is necessary for the dynamic algorithm
*/
ring new_dynamic_ring_from_old(ring oldR, ray *new_weight);

/**
  Converts the elements of an ideal from an old ring to a new ring.
  This is necessary for the dynamic algorithm.
*/
void convertIdeal(ideal & I, ring oldRing, ring newRing);

/**
  Converts oldP=strat->S[i] to the new ring. Since oldP may also be used by
  strat->T and strat->L, it likewise searches through these lists to check
  whether they need changing. In the case of strat->T, some polynomials may
  be missed, so we maintain a list (whichTs) to record which elements of T
  have been changed here.
*/
void convert_stratS(kStrategy strat, bool * & whichTs,
                    ring oldR, ring newR, ring old_tailRing, ring new_tailRing,
                    pShallowCopyDeleteProc p_shallow_copy_delete,
                    omBin new_tailBin);

/**
  Converts strat->P data to new ring.
  newT is the new tailRing.
*/
//void convert_stratP(kStrategy &strat, ring oldR, ring newR, ring newT);
void convert_stratP(kStrategy strat,
                    ring new_tailRing,
                    ring oldR, ring newR,
                    pShallowCopyDeleteProc p_shallow_copy_delete);

/**
  Converts strat[i]->L data to new ring. This includes L.p, but not L.p1 or L.p2,
  as those chang in convert_stratS and convert_stratT.
  However, we do use the value of L.p1 and L.p2 to determine how to change L.p.

  Also converts L.lcm, L.sev, L.FDeg, and L.tailRing.
*/
//void convert_stratL(kStrategy &strat, ring oldR, ring newR);
void convert_stratL(kStrategy strat, ring oldR, ring newR, ring new_tailRing,
                    pShallowCopyDeleteProc p_shallow_copy_delete, poly newTail);

void move_smallest_L_to_back(kStrategy strat);

/**
  Converts those elements of strat->T which are listed as false by whichTs.
  Those listed as true were converted by convert_stratS.
*/
//void convert_stratT(kStrategy & strat, bool * & whichTs, ring oldR, ring newR);
void convert_stratT(kStrategy strat, bool * & whichTs,
          ring oldR, ring newR, ring new_tailRing,
          pShallowCopyDeleteProc p_shallow_copy_delete,
          omBin new_tailBin);

void convert_currentLPPs(vector<poly> & CurrentLPPs, ring oldR, ring newR);

ring new_dynamic_tail_ring(kStrategy strat, ring new_ring);

void dynamic_change_tailRing(kStrategy strat, ring new_tailRing);

#endif // __DUTIL_H