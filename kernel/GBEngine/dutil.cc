/**
  UTILITIES FOR DYNAMIC ALGORITHM
*/

#ifndef __DUTIL_CC_
#define __DUTIL_CC_

#include "dutil.h"

/**
  Creates a new ring from oldR, using new_weight to define the ordering.
*/
ring new_dynamic_ring_from_old(ring oldR, ray *new_weight)
{
  // set up ordering
  int n = oldR->N;
  int ** wvhdl = (int **)omAlloc0(3*sizeof(int_ptr));
  wvhdl[0] = (int *)omAlloc(n*sizeof(int));
  wvhdl[1] = wvhdl[2] = NULL;
  for (int i = 0; i < n; ++i)
    wvhdl[0][i] = (new_weight == NULL) ? 1 : (*new_weight)[i];
  int * order = (int *)omAlloc(2*sizeof(int *));
  int *block0 = (int *)omAlloc0(2*sizeof(int *));
  int *block1 = (int *)omAlloc0(2*sizeof(int *));
  order[0] = ringorder_wp;
  order[1] = ringorder_C;
  block0[0] = block0[1] = 1;
  block1[0] = block1[1] = n;
  order[2] = 0;
  // create ring & return
  ring result = rDefault(nCopyCoeff(oldR->cf), oldR->N, 
                         ((const char **)oldR->names), 4, order,
                         block0, block1, wvhdl);
  rComplete(result, true);
  return result;
}

/**
  Converts the elements of an ideal from an old ring to a new ring.
*/
void convertIdeal(ideal & I, ring oldRing, ring newRing)
{
  if (I != NULL) I = idrMoveR(I, oldRing, newRing);
}

/**
  Resorts a polynomial whose head is in oldR and whose tail is in oldTR
  to a polynomial whose head is in newR and whose tail is in newTR.
  If keep_head is true, this returns a copy of the head in newTR;
  otherwise, it returns NULL.
*/
poly sort_split_poly(poly *pp, ring oldR, ring newR, ring oldTR, ring newTR,
  bool keep_head = false)
{
  poly t = *pp;     // old head
  poly u = t->next; // old tail
  // detach old head, copy into oldTR, to make sorting easier
  t->next = NULL;   // detach
  poly oldP = t;    // remember head in oldR for deletion
  poly oldT = u;    // remember tail in oldTR for deletion
  t = prShallowCopyR(t, oldR, oldTR);
  p_Setm(t, oldTR);
  t->next = u;      // reattach
  // sort into newTR
  t = prShallowCopyR(t, oldTR, newTR);
  p_Setm(t, newTR);
  p_Norm(t, newTR); // normalize
  // detach head and copy into newR
  u = t->next;
  t->next = NULL;
  poly t_in_TR = t; // remember head in newTR
  t = prShallowCopyR(t, newTR, newR);
  p_Setm(t, newR);
  // reattach head
  t->next = u;
  *pp = t;
  // cleanup
  p_ShallowDelete(&oldP, oldR);
  p_ShallowDelete(&oldT, oldTR);
  if (keep_head) {
    t_in_TR->next = u;
    return t_in_TR;
  }
  else {
    p_ShallowDelete(&t_in_TR, newTR);
    return NULL;
  }
}

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
                    omBin new_tailBin)
{
  if (strat->tl >= 0) whichTs = new bool[strat->tl + 1]();
  for (int i = 0; i <= strat->sl; ++i)
  {
    // the data in oldP will soon be corrupted,
    // but we need to remember its address for the sake of strat->L and strat->T
    poly oldP = strat->S[i];
    sort_split_poly(&(strat->S[i]), oldR, newR, strat->tailRing, new_tailRing);
    poly newP = strat->S[i];
    /*poly oldP = strat->S[i];
    poly oldT = pNext(oldP);
    oldP->next = NULL;
    //poly newT = p_shallow_copy_delete(oldT, old_tailRing, new_tailRing, new_tailBin);
    poly newT = prShallowCopyR(oldT, old_tailRing, new_tailRing);
    p_Setm(newT, new_tailRing);
    poly newP = prShallowCopyR(oldP, oldR, newR);
    p_Setm(newP, newR);
    pNext(newP) = newT;
    strat->S[i] = newP;*/
    strat->sevS[i] = pGetShortExpVector(strat->S[i]);
    for (int j = 0; j <= strat->tl; ++j)
      if (strat->T[j].p == oldP)
      {
        whichTs[j] = true;
        strat->T[j].p = newP;
        strat->sevT[j] = pGetShortExpVector(strat->T[j].p);
        strat->T[j].tailRing = new_tailRing;
        strat->T[j].t_p = k_LmInit_currRing_2_tailRing(strat->T[j].p, new_tailRing);
        pNext(strat->T[j].t_p) = pNext(strat->T[j].p);      }
    for (int j = 0; j <= strat->Ll; ++j)
    {
      if (strat->L[j].p1 == oldP)
        strat->L[j].p1 = newP;
      if (strat->L[j].p2 == oldP)
        strat->L[j].p2 = newP;
    }
  }
}

/**
  Converts strat->P data to new ring.
  newT is the new tailRing.
*/
//void convert_stratP(kStrategy &strat, ring oldR, ring newR, ring newT)
void convert_stratP(kStrategy strat,
                    ring new_tailRing,
                    ring oldR, ring newR,
                    pShallowCopyDeleteProc p_shallow_copy_delete)
{
  if (strat->P.t_p != NULL)
  {
    if (strat->P.p != NULL) {
      strat->P.t_p = sort_split_poly(&(strat->P.p),
          oldR, newR, strat->tailRing, new_tailRing, true);
    }
    else {
      poly oldP = strat->P.t_p;
      strat->P.t_p = prShallowCopyR(
          strat->P.t_p, strat->tailRing, new_tailRing);
      p_ShallowDelete(&oldP, strat->tailRing);
    }
    // working with tail ring, data in strat->P.t_p
    /*poly oldT = strat->P.t_p;
    strat->P.t_p = prShallowCopyR(oldT, strat->tailRing, new_tailRing);
    p_ShallowDelete(&oldT, strat->tailRing);
    p_Setm(strat->P.t_p, new_tailRing);
    p_Norm(strat->P.t_p, new_tailRing);
    strat->P.p->next = NULL;
    p_ShallowDelete(&(strat->P.p), oldR);
    poly next = strat->P.t_p->next;
    strat->P.t_p->next = NULL;
    strat->P.p = prShallowCopyR(strat->P.t_p, new_tailRing, newR);
    strat->P.p->next = next;
    strat->P.t_p->next = next;*/
  }
  else if ((strat->P.p != NULL) && pNext(strat->P.p) != strat->tail)
  {
    // strat->P.t_p == NULL: no tail ring?
    poly oldP = strat->P.p;
    strat->P.p = prShallowCopyR(oldP, oldR, newR);
    p_Setm(strat->P.p, newR);
    p_ShallowDelete(&oldP, oldR);
    strat->P.tailRing = newR; // should be newR
    p_Norm(strat->P.p, newR);
    //strat->P.ShallowCopyDelete(new_tailRing, p_shallow_copy_delete);
  }
  /*poly oldP = strat->P.p;
  strat->P.p = prShallowCopyR(oldP, oldR, newR);
  p_Setm(strat->P.p, newR);
  p_ShallowDelete(&oldP, oldR);*/
  strat->P.FDeg = strat->P.pFDeg();
  strat->P.sev = p_GetShortExpVector(strat->P.p, newR);
  strat->P.tailRing = new_tailRing;
  /*strat->P.tailRing = newR;
  p_Norm(strat->P.p, newR);*/
}

/**
  Converts strat[i]->L data to new ring. This includes L.p, but not L.p1 or L.p2,
  as those chang in convert_stratS and convert_stratT.
  However, we do use the value of L.p1 and L.p2 to determine how to change L.p.

  Also converts L.lcm, L.sev, L.FDeg, and L.tailRing.
*/
//void convert_stratL(kStrategy &strat, ring oldR, ring newR)
void convert_stratL(kStrategy strat, ring oldR, ring newR, ring new_tailRing,
                    pShallowCopyDeleteProc p_shallow_copy_delete, poly newTail)
{
  for (int i = 0; i <= strat->Ll; ++i)
  {
    if(strat->L[i].p != NULL)
    {
      // strat->L[i].p1 and strat->L[i].p2 are taken care of with strat->S
      // but if they are nonzero we need a short spoly
      if (strat->L[i].p1 != NULL)
      {
        p_LmDelete(&(strat->L[i].p), oldR);
        strat->L[i].p = ksCreateShortSpoly(
            strat->L[i].p1, strat->L[i].p2, new_tailRing);
        pNext(strat->L[i].p) = newTail;
      }
      else {
        if (strat->L[i].t_p != NULL)
          strat->L[i].t_p = sort_split_poly(&(strat->L[i].p),
              oldR, newR, strat->tailRing, new_tailRing, true);
        else {
          poly oldP = strat->L[i].p;
          strat->L[i].p = prShallowCopyR(oldP, oldR, newR);
          p_Setm(strat->L[i].p, newR);
          p_ShallowDelete(&oldP, oldR);
        }
      }
      /*poly oldP = strat->L[i].t_p;
      strat->L[i].t_p = prShallowCopyR(oldP, strat->tailRing, new_tailRing);
      p_Setm(strat->L[i].t_p, new_tailRing);
      p_ShallowDelete(&oldP, strat->tailRing);
      //resort_poly_in_new_ring(&(strat->L[i].t_p), new_tailRing);
      strat->L[i].p->next = NULL;
      p_ShallowDelete(&(strat->L[i].p), oldR);
      poly t = strat->L[i].t_p;
      poly u = strat->L[i].t_p->next;
      t->next = NULL;
      strat->L[i].p = prShallowCopyR(t, new_tailRing, newR);
      p_Setm(strat->L[i].p, newR);
      t->next = u;
      pNext(strat->L[i].p) = u;
      //p_SortMerge(strat->L[i].p, new_tailRing, true);*/
      (strat->L)[i].FDeg = (strat->L)[i].pFDeg();
      (strat->L)[i].sev = p_GetShortExpVector((strat->L)[i].p,newR);
      strat->L[i].tailRing = new_tailRing;
      /*if (strat->L[i].t_p != NULL)
      {
        poly old_tp = strat->L[i].t_p;
        strat->L[i].t_p = prShallowCopyR(old_tp, strat->tailRing, new_tailRing);
        p_ShallowDelete(&old_tp, strat->tailRing);
      }*/
    }
    /*if ((strat->L)[i].p != NULL)
    {
      poly oldP = (strat->L)[i].p;
      if ((strat->L)[i].p1 != NULL)
      {
        p_LmDelete(&oldP, oldR);
        strat->L[i].p = ksCreateShortSpoly((strat->L)[i].p1, (strat->L)[i].p2, newR);
        pNext(strat->L[i].p) = strat->tail;
      }
      else
      {
        (strat->L)[i].p = prShallowCopyR(oldP, oldR, newR);
        p_ShallowDelete(&oldP, oldR);
        p_Setm((strat->L)[i].p, newR);
      }
      (strat->L)[i].FDeg = (strat->L)[i].pFDeg();
      (strat->L)[i].sev = p_GetShortExpVector((strat->L)[i].p,newR);
      (strat->L)[i].tailRing = newR;
    }*/
    if ((strat->L)[i].lcm != NULL)
    {
      poly oldP = (strat->L)[i].lcm;
      if (!p_GetCoeff(oldP, oldR)) { p_SetCoeff(oldP, n_Init(1,newR), newR); }
      (strat->L)[i].lcm = prShallowCopyR(oldP, oldR, newR);
      p_ShallowDelete(&oldP, oldR);
      p_Setm((strat->L)[i].lcm, newR);
      /*(strat->L)[i].weighted_sugar = p_WDegree((strat->L)[i].lcm,newR);*/
    }
    else
    {
      //(strat->L)[i].weighted_sugar = p_WDegree((strat->L)[i].p,newR);
    }
  }
  strat->L->tailRing = new_tailRing;
}

/**
  Converts those elements of strat->T which are listed as false by whichTs.
  Those listed as true were converted by convert_stratS.
*/
//void convert_stratT(kStrategy & strat, bool * & whichTs, ring oldR, ring newR)
void convert_stratT(kStrategy strat, bool * & whichTs,
          ring oldR, ring newR, ring new_tailRing,
          pShallowCopyDeleteProc p_shallow_copy_delete,
          omBin new_tailBin)
{
  for (int i = 0; i <= strat->tl; ++i)
  {
    bool found = false;
    if ( whichTs != NULL && whichTs[i])
    {
      // this should be taken care of with the corresponding strat->S
      /*strat->T[i].tailRing = new_tailRing;
      strat->T[i].t_p = k_LmInit_currRing_2_tailRing(strat->T[i].p, new_tailRing);
      pNext(strat->T[i].t_p) = pNext(strat->T[i].p);
      */
    }
    else
    {
      poly oldP = (strat->T)[i].p;
      poly oldT = pNext(oldP);
      strat->T[i].t_p = sort_split_poly(&strat->T[i].p,
          oldR, newR, strat->tailRing, new_tailRing, true);
      /*oldP->next = NULL;
      strat->T[i].p = prShallowCopyR(oldP, oldR, newR);
      p_Setm(strat->T[i].p, newR);
      strat->T[i].p->next = prShallowCopyR(oldT, strat->tailRing, new_tailRing);
      p_Setm(strat->T[i].p->next, new_tailRing);
      resort_poly_in_new_ring(&(strat->T[i].p->next), new_tailRing);
      sort_split_poly(&(strat->T[i].p), newR, new_tailRing);*/
      for (int j = 0; j <= strat->Ll; ++j)
      {
        if (strat->L[j].p1 == oldP)
          strat->L[j].p1 = strat->T[i].p;
        if (strat->L[j].p2 == oldP)
          strat->L[j].p2  = strat->T[i].p;
      }
      //oldP->next = NULL;
      //p_ShallowDelete(&oldP, oldR);
      //p_ShallowDelete(&oldT, strat->tailRing);
      /*
      poly oldP = (strat->T)[i].p;
      (strat->T)[i].p = prShallowCopyR(oldP, oldR, newR);
      p_Setm((strat->T)[i].p, newR);
      strat->sevT[i] = pGetShortExpVector(strat->T[i].p);
      (strat->T)[i].t_p = NULL;
      for (int j = 0; j <= strat->Ll; ++j)
      {
        if (strat->L[j].p1 == oldP)
        {
          //cout << "changing pair poly 1  of " << j << ' '; pWrite(oldP);
          strat->L[j].p1 = strat->T[i].p;
        }
        if (strat->L[j].p2 == oldP)
        {
          //cout << "changing pair poly 2 of " << j << ' '; pWrite(oldP);
          strat->L[j].p2 = strat->T[i].p;
        }
      }
      p_ShallowDelete(&oldP, oldR);
      */
      //cout << "becomes "; pWrite((strat->T)[i].p);
    }
    //(strat->T)[i].tailRing = newR;
    (strat->T)[i].tailRing = new_tailRing;
    (strat->T)[i].FDeg = (strat->T)[i].pFDeg();
  }
  //if (strat->tl >= 0)
  //  omFree(whichTs);
  omMergeStickyBinIntoBin(strat->tailBin, strat->tailRing->PolyBin);
  //if (strat->tailRing != oldR) rKillModifiedRing(strat->tailRing);
  strat->tailRing = new_tailRing;
  strat->tailBin = new_tailBin;
  strat->p_shallow_copy_delete = p_shallow_copy_delete;
  if (whichTs != NULL) delete [] whichTs;
  kTest_TS(strat);
}

void convert_currentLPPs(vector<poly> & CurrentLPPs, ring oldR, ring newR)
{
  for (vector<poly>::iterator pp = CurrentLPPs.begin(); pp != CurrentLPPs.end(); ++pp)
  {
    poly oldP = *pp;
    *pp = prShallowCopyR(oldP, oldR, newR);
    p_Setm(*pp, newR);
    p_ShallowDelete(&oldP, oldR);
  }
}

ring new_dynamic_tail_ring(kStrategy strat, ring new_ring)
{
  /*ring tailRing = rCopy(strat->tailRing);
  int n = new_ring->N;
  for (int i = 0; i < n; ++i)
  {
    tailRing->wvhdl[0][i] = new_ring->wvhdl[0][i];
  }*/
  ring tailRing = rModifyRing_Wp(strat->tailRing, new_ring->wvhdl[0]);
  rComplete(tailRing, true);
  return tailRing;
}

ring old_dynamic_tail_ring(kStrategy strat, ring new_ring)
{
  // from kStratInitChangeTailRing(),
  unsigned long e = 0;
  int i;
  for (i = 0; i <= strat->Ll; ++i)
    e = p_GetMaxExpL(strat->L[i].p, new_ring, e);
  for (i = 0; i <= strat->tl; ++i)
    e = p_GetMaxExpL(strat->T[i].p, new_ring, e);
  e = p_GetMaxExp(e, new_ring);
  if (e <= 1) e = 2;

  // from rModifyRing()
  // as this is a special case, we set things up differently from rModifyRing()
  int n = new_ring->N;
  int **wvhdl = (int **)omAlloc0(3*sizeof(int_ptr));
  wvhdl[1] = wvhdl[2] = NULL;
  wvhdl[0] = (int *)omAlloc(n*sizeof(int));
  for (int i = 0; i < n; ++i)
    wvhdl[0][i] = (new_ring->wvhdl)[0][i];
  //wvhdl[0] = (new_ring->wvhdl)[0];
  int *order = (int *)omAlloc0(2*sizeof(int));
  int *block0 = (int *)omAlloc0(2*sizeof(int));
  int *block1 = (int *)omAlloc0(2*sizeof(int));
  order[0] = ringorder_wp;
  order[1] = ringorder_C;
  block0[0] = block0[1] = 1;
  block1[0] = block1[1] = n;
  order[2] = 0;
  ring result = (ring)omAlloc0Bin(sip_sring_bin);
  *result = *new_ring;
  result->wvhdl   = wvhdl;
  result->order   = order;
  result->block0  = block0;
  result->block1  = block1;
  result->bitmask = e;
  result->pFDegOrig = result->pFDeg = p_Deg;
  result->pLDegOrig = result->pLDeg = pLDegb;
  // back to rModifyRing()
  if (result->typ != NULL)
  {
    if (result->typ[0].ord_typ == ro_syz)
    {
      result->typ[0] = new_ring->typ[0];
      if (result->typ[0].data.syz.limit > 0)
      {
        result->typ[0].data.syz.syz_index
          = (int *)omAlloc((new_ring->typ[0].data.syz.limit + 1)*sizeof(int));
        memcpy(result->typ[0].data.syz.syz_index, new_ring->typ[0].data.syz.syz_index,
          (result->typ[0].data.syz.limit+1)*sizeof(int));
      }
    }
  }
  result->OrdSgn = new_ring->OrdSgn;
  rComplete(result, true); // force complete

  return result;
}

void dynamic_change_tailRing(kStrategy strat, ring new_tailRing)
{
  ring old_tailRing = (strat->tailRing == NULL) ? currRing : strat->tailRing;
  if (new_tailRing == old_tailRing) return;
  strat->pOrigFDeg_TailRing = new_tailRing->pFDeg;
  strat->pOrigLDeg_TailRing = new_tailRing->pLDeg;
  if (old_tailRing->pFDeg != old_tailRing->pFDegOrig)
  {
    new_tailRing->pFDeg = old_tailRing->pFDeg;
    new_tailRing->pLDeg = old_tailRing->pLDeg;
  }
  kTest_TS(strat);
  assume(new_tailRing != strat->tailRing);
  pShallowCopyDeleteProc
      p_shallow_copy_delete = pGetShallowCopyDeleteProc(strat->tailRing,
                                                        new_tailRing);
  omBin new_tailBin = omGetStickyBinOfBin(new_tailRing->PolyBin);
  for (int i = 0; i <= strat->tl; ++i)
    strat->T[i].ShallowCopyDelete(new_tailRing, new_tailBin,
                                  p_shallow_copy_delete);
  for (int i = 0; i <= strat->Ll; ++i)
  {
    assume(strat->L[i].p != NULL);
    if (pNext(strat->L[i].p) != strat->tail)
      strat->L[i].ShallowCopyDelete(new_tailRing, p_shallow_copy_delete);
  }
  if ((strat->P.t_p != NULL) ||
      ((strat->P.p != NULL) && pNext(strat->P.p) != strat->tail))
    strat->P.ShallowCopyDelete(new_tailRing, p_shallow_copy_delete);

  // no arguments L or T, so we skip corresponding part

  omMergeStickyBinIntoBin(strat->tailBin, strat->tailRing->PolyBin);
  if (strat->tailRing != old_tailRing) rKillModifiedRing(strat->tailRing);
  strat->tailRing = new_tailRing;
  strat->tailBin = new_tailBin;
  strat->p_shallow_copy_delete = pGetShallowCopyDeleteProc(old_tailRing,
                                                           new_tailRing);
  if (strat->kHEdge != NULL)
  {
    if (strat->t_kHEdge != NULL)
      p_LmFree(strat->t_kHEdge, strat->tailRing);
    strat->kHEdge = k_LmInit_currRing_2_tailRing(strat->kHEdge, new_tailRing);
  }
  if (strat->kNoether != NULL)
  {
    if (strat->t_kNoether != NULL)
      p_LmFree(strat->t_kNoether, strat->tailRing);
    strat->t_kNoether = k_LmInit_currRing_2_tailRing(strat->kNoether,
                                                     new_tailRing);
  }
  kTest_TS(strat);
}

#endif