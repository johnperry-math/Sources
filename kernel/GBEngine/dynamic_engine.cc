#ifndef __DYNAMIC_ENGINE_C__
#define __DYNAMIC_ENGINE_C__

#include "dynamic_engine.h"

#include <coeffs/coeffs.h>
#include <polys/coeffrings.h>
#include <polys/monomials/p_polys.h>
#include <polys/prCopy.h>

#include <iostream>
using std::cout; using std::endl;

#include <list>
using std::list;

ring QQt = NULL;

bool LessByHilbert (PPWithIdeal &a, PPWithIdeal &b)
{
  bool result;
  ring Rx = currRing;
  int n = Rx->N;
  // first check the coefficients of the Hilbert polynomial
  poly HPdiff = p_Minus_mm_Mult_qq(p_Copy(a.getHilbertPolynomial(), QQt), p_One(QQt), b.getHilbertPolynomial(), QQt);
  if (HPdiff != NULL)
    result = !(n_IsZero(p_GetCoeff(HPdiff, QQt), QQt) || n_GreaterZero(p_GetCoeff(HPdiff, QQt), QQt));
  else // use Hilbert series
  {
    intvec * h1 = a.getHilbertNumerator();
    intvec * h2 = b.getHilbertNumerator();
    int i = 1;
    for ( /* already initialized */ ;
         i < h1->length() && i < h2->length() && (*h1)[i] == (*h2)[i];
         i++)
    { /* taken care of in loop */ }
    if (i >= h1->length())
    {
      if (i >= h2->length())
      { // the numerators are equal; break tie via lex
        int i = 1;
        while (i <= n && p_GetExp(a.getPP(),i,Rx) == p_GetExp(b.getPP(),i,Rx))
          ++i;
        if (i == n)
          result = false;
        else
          result = p_GetExp(a.getPP(), i, Rx) > p_GetExp(b.getPP(), i, Rx);
      }
      else
        result = true;
    }
    else
    {
      if (i > h2->length())
        result = false;
      else
        result = (*h1)[i] < (*h2)[i];
    }
  }
  //cout << "\tfirst less than second? " << result << endl;
  return result;
};

bool LessBySmoothestDegrees (PPWithIdeal &a, PPWithIdeal &b)
{
  return a.getDifferenceInDegree() < b.getDifferenceInDegree();
};

bool LessByLargestMaxComponent (PPWithIdeal &a, PPWithIdeal &b)
{
  return a.getDifferenceInDegree() < b.getDifferenceInDegree();
}

void PPWithIdeal::computeNumberNewPairs()
{
  ring Rx = currRing;
  int n = Rx->N;
  num_new_pairs = min_deg = 0;
  bool * keepers = new bool [strat->sl+1];
  for (int i = 0; i < strat->sl + 1; ++i) keepers[i] = true;
  // first main loop: apply Buchberger's lcm criterion to new pairs
  for (int i = 0; i < strat->sl + 1; ++i)
  {
    // if gcd(t,Si) == 1 then skip i for the time being
    if (!pHasNotCF(t,strat->S[i]))
    {
      bool has_divisor = false;
      for (int j=0; (!has_divisor) && j < strat->sl + 1; ++j)
      {
        if (i != j && keepers[j])
        {
          // if some j satisfies lcm(t,Sj) | lcm(t,Si) then do not count i
          has_divisor = true;
          for (int k=1; has_divisor && k <= n; ++k)
            if // deg(lcm(t,lm(Si)))
                (((pGetExp(t,k) > pGetExp(strat->S[i],k)) ?
                  pGetExp(t,k) : pGetExp(strat->S[i],k))
               < // deg(lcm(t,lm(Sj)))
                ((pGetExp(t,k) > pGetExp(strat->S[j],k)) ?
                  pGetExp(t,k) : pGetExp(strat->S[j],k)))
              has_divisor = false;
        }
      }
    }
  }
  // second main loop: apply Buchberger's gcd criterion to new pairs, count survivors
  for (int i = 0; i < strat->sl + 1; ++i)
  {
    if (keepers[i] && !pHasNotCF(t,strat->S[i]))
    {
      int new_deg = 0;
      // determine deg(lcm(t,Si))
      for (int k=1; k <= n; ++k)
        new_deg += (pGetExp(t,k) > pGetExp(strat->S[i],k))
          ? pGetExp(t,k)
          : pGetExp(strat->S[i],k);
      /*for (int k=1; k <= n; ++k)
        new_deg += Rx->wvhdl[0][k] *
          ((pGetExp(t,k) > pGetExp(strat->S[i],k)) ?
              pGetExp(t,k) : pGetExp(strat->S[i],k));*/
      if (min_deg == 0 || min_deg > new_deg)
      {
        min_deg = new_deg; num_new_pairs = 1;
      }
      else if (min_deg == new_deg)
      {
        ++num_new_pairs;
        //cout << '\t'; pWrite(pHead(strat->S[i]));
      }
    }
  }
  delete [] keepers;
  // cout << "we get " << num_new_pairs << " from "; pWrite(t);
  // third main loop: apply Buchberger's lcm criterion to old pairs, UNLESS
  // all three lcm's are equal
  for (int i = 0; i < strat->Ll + 1; ++i)
  {
    if (strat->L[i].lcm != NULL)
    {
      poly u = strat->L[i].lcm;
      int new_deg = 0;
      for (int k=1; k <= n; ++k) new_deg += pGetExp(u,k);
      // for (int k=1; k <= n; ++k) new_deg += pGetExp(u,k) * Rx->wvhdl[0][k];
      // no point continuing if it wouldn't change min_deg
      if (min_deg == 0 || new_deg <= min_deg)
      {
        // see if new poly divides this one
        bool has_divisor = true;
        // see if Buchberger triple has same lcm
        bool all_equal = true;
        poly p1 = strat->L[i].p1;
        poly p2 = strat->L[i].p2;
        for (int k=1; has_divisor && all_equal && k <= n; ++k)
        {
          if (pGetExp(t,k) > pGetExp(u,k)) has_divisor = false;
          else
          {
            // check lcm(t,lm(p1)) == lcm(t,lm(p2)) == lcm(lm(p1),lm(p2)) in xk
            int a = (pGetExp(t,k) > pGetExp(p1,k)) ? pGetExp(t,k) : pGetExp(p1,k);
            int b = (pGetExp(t,k) > pGetExp(p2,k)) ? pGetExp(t,k) : pGetExp(p2,k);
            int c = (pGetExp(p1,k) > pGetExp(p2,k)) ? pGetExp(p1,k) : pGetExp(p2,k);
            all_equal = (a == c) && (b == c);
          }
        }
        if (!has_divisor || all_equal)
        {
          if (min_deg == 0 || min_deg > new_deg)
          {
            min_deg = new_deg; num_new_pairs = 1;
          }
          else // the only reason we'd be here is if min_deg == new_deg
          {
            ++num_new_pairs;
            //cout << '\t'; pWrite(u);
          }
        }
      }
    }
  }
  // cout << " which makes " << num_new_pairs << " pairs total at degree " << min_deg << ".\n";
}

bool LessByNumCritPairs (PPWithIdeal &a, PPWithIdeal &b)
{
  bool result;
  ring Rx = currRing;
  int n = Rx->N;
  // first check if the number of critical pairs has been computed
  if (a.howManyNewPairs() < 0) a.computeNumberNewPairs();
  if (b.howManyNewPairs() < 0) b.computeNumberNewPairs();
  if (a.degOfNewPairs() < b.degOfNewPairs())
    result = true;
  else if (a.degOfNewPairs() > b.degOfNewPairs())
    result = false;
  // at this point, the degrees of the new pairs will be equal
  else if (a.howManyNewPairs() > b.howManyNewPairs())
    result = true;
  else if (a.howManyNewPairs() < b.howManyNewPairs())
    result = false;
  else
  { // the numerators are equal; break tie via lex
    int i = 1;
    while (i <= n && p_GetExp(a.getPP(),i,Rx) == p_GetExp(b.getPP(),i,Rx))
      ++i;
    if (i == n)
      result = false;
    else
      result = p_GetExp(a.getPP(), i, Rx) > p_GetExp(b.getPP(), i, Rx);
  }
  //cout << "\tfirst less than second? " << result << endl;
  return result;
};

bool LessByHilbertThenDegree(PPWithIdeal &a, PPWithIdeal &b)
{
  bool result;
  ring Rx = currRing;
  int n = Rx->N;
  // first check the coefficients of the Hilbert polynomial
  poly HPdiff = p_Minus_mm_Mult_qq(p_Copy(a.getHilbertPolynomial(), QQt), p_One(QQt), b.getHilbertPolynomial(), QQt);
  if (HPdiff != NULL)
    result = !(n_IsZero(p_GetCoeff(HPdiff, QQt), QQt) || n_GreaterZero(p_GetCoeff(HPdiff, QQt), QQt));
  else // use Hilbert series
  {
    intvec * h1 = a.getHilbertNumerator();
    intvec * h2 = b.getHilbertNumerator();
    int i = 1;
    for ( /* already initialized */ ;
         i < h1->length() && i < h2->length() && (*h1)[i] == (*h2)[i];
         i++)
    { /* taken care of in loop */ }
    if (i >= h1->length())
    {
      if (i >= h2->length())
      {
        if (p_WDegree(a.getPP(), Rx) < p_WDegree(b.getPP(), Rx))
          result = true;
        else if (p_WDegree(a.getPP(), Rx) > p_WDegree(b.getPP(), Rx))
          result = false;
        else
        {
          int i = 1;
          while (i <= n && p_GetExp(a.getPP(),i,Rx) == p_GetExp(b.getPP(),i,Rx))
            ++i;
          if (i == n)
            result = false;
          else
            result = p_GetExp(a.getPP(), i, Rx) < p_GetExp(b.getPP(), i, Rx);
        }
      }
      else
        result = true;
    }
    else
    {
      if (i > h2->length())
        result = false;
      else
        result = (*h1)[i] < (*h2)[i];
    }
  }
  return result;
};

bool LessByGradHilbertThenDegree(PPWithIdeal &a, PPWithIdeal &b)
{
  bool result;
  ring Rx = currRing;
  int n = Rx->N;
  // first check the coefficients of the Hilbert polynomial
  poly HPdiff = p_Minus_mm_Mult_qq(p_Copy(a.getHilbertPolynomial(), QQt), p_One(QQt), b.getHilbertPolynomial(), QQt);
  if (HPdiff != NULL)
    result = !(n_IsZero(p_GetCoeff(HPdiff, QQt), QQt) || n_GreaterZero(p_GetCoeff(HPdiff, QQt), QQt));
  else // use Hilbert series
  {
    intvec * weights = new intvec(n+1);
    for (int i=0; i < n; ++i) (*weights)[i] = Rx->wvhdl[0][i]; 
    intvec * h1 = a.getHilbertNumerator(weights);
    intvec * h2 = b.getHilbertNumerator(weights);
    delete weights;
    int i = 1;
    for ( /* already initialized */ ;
         i < h1->length() && i < h2->length() && (*h1)[i] == (*h2)[i];
         i++)
    { /* taken care of in loop */ }
    if (i >= h1->length())
    {
      if (i >= h2->length())
      {
        if (p_WDegree(a.getPP(), Rx) < p_WDegree(b.getPP(), Rx))
          result = true;
        else if (p_WDegree(a.getPP(), Rx) > p_WDegree(b.getPP(), Rx))
          result = false;
        else
        {
          int i = 1;
          while (i <= n && p_GetExp(a.getPP(),i,Rx) == p_GetExp(b.getPP(),i,Rx))
            ++i;
          if (i == n)
            result = false;
          else
            result = p_GetExp(a.getPP(), i, Rx) < p_GetExp(b.getPP(), i, Rx);
        }
      }
      else
        result = true;
    }
    else
    {
      if (i > h2->length())
        result = false;
      else
        result = (*h1)[i] < (*h2)[i];
    }
  }
  return result;
};

bool LessByDegreeThenHilbert(PPWithIdeal &a, PPWithIdeal &b)
{
  bool result;
  ring Rx = currRing;
  int n = Rx->N;
  // first check the weighted degree
  if (p_WDegree(a.getPP(), Rx) < p_WDegree(b.getPP(), Rx))
    result = true;
  else if (p_WDegree(a.getPP(), Rx) > p_WDegree(b.getPP(), Rx))
    result = false;
  else {
    // now check the coefficients of the Hilbert polynomial
    poly HPdiff = p_Minus_mm_Mult_qq(p_Copy(a.getHilbertPolynomial(), QQt), p_One(QQt), b.getHilbertPolynomial(), QQt);
    if (HPdiff != NULL)
      result = !(n_IsZero(p_GetCoeff(HPdiff, QQt), QQt) || n_GreaterZero(p_GetCoeff(HPdiff, QQt), QQt));
    else // use Hilbert series
    {
      intvec * weights = new intvec(n+1);
      for (int i=0; i < n; ++i) (*weights)[i] = Rx->wvhdl[0][i];
      intvec * h1 = a.getHilbertNumerator(weights);
      //cout << "have h1 "; for (int i=0; i < h1->length(); ++i) cout << (*h1)[i] << ' '; cout << endl;
      intvec * h2 = b.getHilbertNumerator(weights);
      //cout << "have h2 "; for (int i=0; i < h2->length(); ++i) cout << (*h2)[i] << ' '; cout << endl;
      delete weights;
      //cout << "deleted weights\n";
      int i = 1;
      for ( /* already initialized */ ;
           i < h1->length() && i < h2->length() && (*h1)[i] == (*h2)[i];
           i++)
      { /* taken care of in loop */ }
      if (i >= h1->length())
      {
        if (i >= h2->length())
        {
          int i = 1;
          while (i <= n && p_GetExp(a.getPP(),i,Rx) == p_GetExp(b.getPP(),i,Rx))
            ++i;
          if (i == n)
            result = false;
          else
            result = p_GetExp(a.getPP(), i, Rx) < p_GetExp(b.getPP(), i, Rx);
        }
        else
          result = true;
      }
      else
      {
        if (i > h2->length())
          result = false;
        else
          result = (*h1)[i] < (*h2)[i];
      }
    }
  }
  return result;
};

bool LessByDegreeThenGradHilbert(PPWithIdeal &a, PPWithIdeal &b)
{
  bool result;
  ring Rx = currRing;
  int n = Rx->N;
  // first check the weighted degree
  if (p_WDegree(a.getPP(), Rx) < p_WDegree(b.getPP(), Rx))
    result = true;
  else if (p_WDegree(a.getPP(), Rx) > p_WDegree(b.getPP(), Rx))
    result = false;
  else {
    // now check the coefficients of the Hilbert polynomial
    poly HPdiff = p_Minus_mm_Mult_qq(p_Copy(a.getHilbertPolynomial(), QQt), p_One(QQt), b.getHilbertPolynomial(), QQt);
    if (HPdiff != NULL)
      result = !(n_IsZero(p_GetCoeff(HPdiff, QQt), QQt) || n_GreaterZero(p_GetCoeff(HPdiff, QQt), QQt));
    else // use Hilbert series
    {
      intvec * h1 = a.getHilbertNumerator();
      intvec * h2 = b.getHilbertNumerator();
      int i = 1;
      for ( /* already initialized */ ;
           i < h1->length() && i < h2->length() && (*h1)[i] == (*h2)[i];
           i++)
      { /* taken care of in loop */ }
      if (i >= h1->length())
      {
        if (i >= h2->length())
        {
          int i = 1;
          while (i <= n && p_GetExp(a.getPP(),i,Rx) == p_GetExp(b.getPP(),i,Rx))
            ++i;
          if (i == n)
            result = false;
          else
            result = p_GetExp(a.getPP(), i, Rx) < p_GetExp(b.getPP(), i, Rx);
        }
        else
          result = true;
      }
      else
      {
        if (i > h2->length())
          result = false;
        else
          result = (*h1)[i] < (*h2)[i];
      }
    }
  }
  return result;
};

/*poly monomial(int a, int *e, ring R)
{ poly p = p_One(R); p_SetExpV(p, e, R); p_SetCoeff(p, n_Init(a, R), R); p_Setm(p, R); return p; }

poly monomial_univariate(number a, ring R)
{ poly p = p_One(R); p_SetExp(p, 1, 0, R); p_SetCoeff(p, n_Init(a, R), R); p_Setm(p, R); return p;}*/

poly var_univariate(ring R)
{ poly p = p_One(R); p_SetExp(p, 1, 1, R); p_Setm(p, R); return p; }

poly constant(number a, ring R)
{ poly p = p_One(R); p_SetCoeff(p, a, R); p_Setm(p, R); return p; }

const char * tstring = "t";

/**
  Computes C(t+a,b).
*/
poly poly_binomial(int a, int b)
{
  poly p;
  if (QQt == NULL) { QQt = rDefault(0, 1, &tstring); }
  if (b == 0) { p = p_One(QQt); }
  else
  {
    // p = (t + a) / b
    poly t = var_univariate(QQt);
    p = p_Copy_noCheck(t, QQt);
    poly ca;
    if (a != 0) ca = constant(n_Init(a,QQt), QQt);
    else ca = NULL;
    p = p_Add_q(p, ca, QQt);
    p = p_Div_nn(p, n_Init(b,QQt), QQt);
    // p = p * (t + a - i) / i for i = 1, ..., b - 1
    poly q;
    for (int i = 1; i < b; ++i)
    {
      if (a != i) ca = constant(n_Init(a-i,QQt), QQt);
      else ca = NULL;
      q = p_Copy_noCheck(t, QQt); q = p_Add_q(q, ca, QQt);
      q = p_Div_nn(q, n_Init(i,QQt), QQt);
      p = p_Mult_q(p, q, QQt);
    }
  }
  return p;
}

poly HilbertPolynomial(intvec * hilbertNumerator, int dim)
{
  if (QQt == NULL)
  { QQt = rDefault(0, 1, &tstring); }
  //poly hp = constant(n_Init(0,QQt), QQt);
  poly hp = NULL;
  if (dim != 0)
  {
    int d = hilbertNumerator->length() - 2;
    poly q = poly_binomial(dim - 1 - d, dim - 1);
    poly r = (*hilbertNumerator)[d] == 0 ? NULL : constant(n_Init((*hilbertNumerator)[d], QQt), QQt);
    q = p_Mult_q(q, r, QQt);
    hp = p_Add_q(hp, q, QQt);
    for (int i = 1; i <= d; ++i)
    {
      q = poly_binomial(dim - 1 - d + i, dim - 1);
      r = (*hilbertNumerator)[d - i] == 0 ? NULL : constant(n_Init((*hilbertNumerator)[d - i], QQt), QQt);
      q = p_Mult_q(q, r, QQt);
      hp = p_Add_q(hp, q, QQt);
    }
    // compute factorial(dim) if you want integers, not otherwise
    /*int n = 1;
    for (int i = 2; i <= dim - 1; ++i) { n *= i; }
    q = constant(n_Init(n, QQt), QQt); 
    // multiply hp by this number
    hp = p_Mult_q(hp, q, QQt); */
  }
  return hp;
}

void compatiblePP(
  poly currentLPP,            // the current LPP
  const set<poly> &allPPs,    // the monomials to consider;
                                      // some will be removed
  const set<ray> &bndrys,     // known boundary vectors
  set<poly> &result,          // returned as PPs for Hilbert function
                                      // ("easy" (& efficient?) to extract exps
  set<poly> &boundary_mons,   // boundary monomials
  skeleton &skel              // used for alternate refinement
)
{
  // get the exponent vector of the current LPP, insert it
  ring Rx = currRing;
  unsigned long n = Rx->N;
  int * a = new int[n+1];
  unsigned long long * along = new unsigned long long [n];
  int * b = new int[n+1];
  p_GetExpV(currentLPP, a, Rx);
  for (unsigned long i = 1; i <= n; ++i) { along[i-1] = a[i]; }
  ray aray(n, along);
  set<poly> initial_candidates;
  initial_candidates.insert(currentLPP);
  // compare other monomials with LPP
  //cout << "initial test of monomials\n";
  /*for (
        set<poly>::iterator b_ptr = allPPs.begin();
        b_ptr != allPPs.end();
        ++b_ptr
      )
  {
    // take the dot product of each monomial's exponents with pp,
    // add the ones that pass muster to initial_candidates
    //cout << '\t'; p_Write(*b_ptr, Rx);
    p_GetExpV(*b_ptr, b, Rx);
    for (unsigned long i = 1; i <= n; ++i) { along[i-1] = b[i]; }
    ray bray(n, along);
    bool searching = true;
    for (
         set<ray>::iterator w_ptr = bndrys.begin();
         searching && w_ptr != bndrys.end();
         ++w_ptr
        )
    {
      // check b against a with all rays
      //cout << "checking " << bray << " against " << aray << " via " << *w_ptr << ": " << (*w_ptr) * bray << ' ' << (*w_ptr) * aray << endl;
      if ((*w_ptr) * bray > (*w_ptr) * aray)
      {
        // only need one success
        initial_candidates.insert(*b_ptr);
        cout << "succeeded with " << bray << endl;
        searching = false;
      }
    }
  }*/
  for (set<poly>::iterator b_ptr = allPPs.begin(); b_ptr != allPPs.end(); ++b_ptr)
  {
    if (*b_ptr != currentLPP && skel.makes_consistent_constraint(*b_ptr, currentLPP))
      initial_candidates.insert(*b_ptr);
  }
  result = initial_candidates;
  // (new) alternate refinement: compare remaining monomials against each other,
  // using skeleton to ensure consistency
  // this approach should be more efficient than the one below, yet equivalent to it
  for (set<poly>::iterator t = initial_candidates.begin();
       t != initial_candidates.end();
       ++t)
  {
    // cout << "testing for consistency: "; pWrite(*t);
    bool good_constraints = true;
    for (set<poly>::iterator u = initial_candidates.begin();
         good_constraints && u != initial_candidates.end();
         ++u)
      if (*t != *u)
        if (!skel.makes_consistent_constraint(*t,*u))
          good_constraints = false;
    if (good_constraints)
    {
      boundary_mons.insert(*t);
      // cout << "\tconsistent!\n";
    }
  }
  // second refinement: compare remaining monomials with each other
  /*for (set<poly>::iterator p_ptr = initial_candidates.begin(); p_ptr != initial_candidates.end(); ++p_ptr)
  {
    p_GetExpV(*p_ptr, a, Rx);
    for (unsigned long i = 1; i <= n; ++i) { along[i-1] = a[i]; }
    ray pray(n, along);
    bool cleared_the_hurdle = false;
    for (set<ray>::iterator w_ptr = bndrys.begin();
         !cleared_the_hurdle && w_ptr != bndrys.end(); ++w_ptr)
    {
      bool maybe_this_vector = true;
      ray w = *w_ptr;
      for (set<poly>::iterator q_ptr = initial_candidates.begin();
           maybe_this_vector && !cleared_the_hurdle && q_ptr != initial_candidates.end(); ++q_ptr)
        if (*p_ptr != *q_ptr)
        {
          p_GetExpV(*q_ptr, b, Rx);
          for (unsigned long i = 1; i <= n; ++i) { along[i-1] = b[i]; }
          ray qray(n, along);
          // guard against invalid exclusions
          unsigned long long wt = (*w_ptr) * qray;
          if ((((*w_ptr) * pray) <= wt) && wt != 0) { maybe_this_vector = false; }
        }
      cleared_the_hurdle = maybe_this_vector;
      if (!cleared_the_hurdle) boundary_mons.insert(*p_ptr);
    }
    if (cleared_the_hurdle) result.insert(*p_ptr);
  } */
  //cout << "boundary monomials:\n";
  //for (set<poly>::iterator piter = boundary_mons.begin(); piter != boundary_mons.end(); ++piter) pWrite(*piter);
  //cout << "compatible monomials:\n";
  //for (set<poly>::iterator piter = result.begin(); piter != result.end(); ++piter) pWrite(*piter);
  //result = initial_candidates;
  delete [] a; delete [] b; delete [] along;
}

bool verifyAndModifyIfNecessary(
  ring Tx,
  skeleton &skel,
  //const vector<poly> &currentPolys
  polyset currentPolys,
  int numPolys
)
{
  ring Rx = currRing;
  bool consistent = true; // innocent until proven guilty
  ray w = ray_sum(skel.get_rays()); // our tentative ordering
  unsigned long n = w.get_dimension();
  //ring Rx = currRing;
  int * pExp = new int[n + 1];
  int * qExp = new int[n + 1];
  unsigned long long *entries = new unsigned long long [n]; // used for entries for new rays
  long *coefficients = new long [n]; // used for coefficients for new constraints
  // loop through all polynomials; verify leading power product is unchanged
  for (
       //vector<const poly>::iterator piter = currentPolys.begin();
       //piter != currentPolys.end();
       //++piter
       int si = 0; si < numPolys; ++si
      )
  {
    // create a ray for the current LPP's exponents
    //poly pp = p_Head(*piter, Rx);
    poly pp = currentPolys[si];
    p_GetExpV(pp, pExp, Rx);
    for (unsigned long i = 0; i < n; ++i)
      entries[i] = pExp[i + 1];
    ray a(n, entries);
    // loop through the polynomial's remaining monomials
    //for (poly titer = (*piter)->next; consistent and titer != NULL; titer = titer->next)
    for (poly titer = pp->next; consistent and titer != NULL; titer = titer->next)
    {
      if (pp != titer) // don't compare with LPP; that would be Bad (TM)
      {
        // create a ray for the PP's exponents
        p_GetExpV(titer, qExp, Tx);
        for (unsigned long i = 0; i < n; ++i)
          entries[i] = qExp[i + 1];
        ray b(n, entries);
        // compare weights between a and b; if this fails,
        // recompute the skeleton with a new constraint
        if (a*w <= b*w)
        {
          if (coefficients == NULL) // ensure we have space
            coefficients = new long[n];
          for (unsigned long i = 0; i < n; ++i)
            coefficients[i] = a[i] - b[i];
          constraint new_constraint(n, coefficients);
          skeleton newskel(skel);
          consistent = newskel.ddm(new_constraint);
          w = ray_sum(skel.get_rays());
          // if we're consistent, we need to recompute the ordering
          if (consistent && a*w > b*w)
          {
            w = ray_sum(skel.get_rays());
            skel = newskel;
          } // if consistent
          else consistent = false;
        } // if LPP changed
      } // if PP != LPP
    } // loop through PPs
  } // loop through polys
  // cleanup
  delete [] entries;
  delete [] coefficients;
  delete [] pExp;
  delete [] qExp;
  // finally done
  return consistent;
}

void ConstraintsForNewPP(
  const PPWithIdeal &I,
  const set<poly> &monomialsForComparison,
  vector<constraint> &result
)
{
  // setup
  ring Rx = currRing;
  unsigned long n = Rx->N;
  int * a = new int[n+1];   // space for exponent vectors
  int * b = new int[n+1];
  long *c = new long[n+1];  // space for coefficients of constraint
  p_GetExpV(I.getPP(), a, Rx);  // exponent vector of candidate
  // loop through exponent vectors of other 
  for (
       set<poly>::iterator miter = monomialsForComparison.begin();
       miter != monomialsForComparison.end();
       ++miter
      )
  {
    // insert only different PPs (since I->t should also be in that set)
    //if ((*miter) != I.getPP())
    if (!p_EqualPolys(*miter, I.getPP(), Rx))
    {
      p_GetExpV(*miter, b, Rx);
      for (ulong i = 0; i < n; ++i)
        c[i] = a[i+1] - b[i+1];
      result.push_back(constraint(n,c));
    }
  }
  delete [] c;
  delete [] b;
  delete [] a;
}

void SelectMonomial(
    ring Tx,
    poly &r,                          // changes
    vector<poly> &CurrentLPPs,      // changes
    //const vector<poly> &CurrentPolys,
    polyset CurrentPolys,
    int numPolys,
    skeleton & currSkel,                // possibly changes
    bool &ordering_changed,
    skStrategy * strat,
    DynamicHeuristic method
)
{
  cout << "entering selmon\n";
  //cout << "current ring: "; rWrite(currRing); cout << endl;
  //cout << "new polynomial: "; pWrite(r); cout << endl;
  //cout << "skeleton before: " << currSkel << endl;
  ring Rx = currRing;
  ray w = ray_sum(currSkel.get_rays());
  vector<unsigned long long> ord;
  for (unsigned long i = 0; i < w.get_dimension(); ++i) { ord.push_back(w[i]); }
  poly currentLPP = p_Head(r, Rx);
  //cout << "comparing against: "; p_Write(currentLPP, Rx);
  set<poly> allPPs, boundaryPPs, compatiblePPs;
  allPPs.insert(currentLPP);
  // transform monomials into exponent vectors
  for (
       poly titer = pNext(r);
       titer != NULL;
       titer = titer->next
      )
  {
    //poly t = pHead(titer);
    poly t = prHeadR(titer, Tx, Rx);
    allPPs.insert(t);
    //allPPs.insert(p_Head(titer, Rx));
  }
  // loop through all exponent vectors
  cout << allPPs.size() << " possible monomials\n";
  compatiblePP(currentLPP, allPPs, currSkel.get_rays(), compatiblePPs, boundaryPPs, currSkel);
  cout << compatiblePPs.size() << " compatible monomials\n";
  //for (set<poly>::iterator piter = compatiblePPs.begin(); piter != compatiblePPs.end(); ++piter)
  //  p_Write(*piter, Rx);
  // list possible future ideals, sort by Hilbert Function
  // using a set sorts the ideals automagically by appropriate DynamicHeuristics
  list<PPWithIdeal> possibleIdealsBasic;
  for (
       set<poly>::iterator piter = compatiblePPs.begin();
       piter != compatiblePPs.end();
       ++piter
      )
  {
    PPWithIdeal newIdeal(*piter, CurrentLPPs, w, strat);
    possibleIdealsBasic.push_back(newIdeal);
  }
  switch(method)
  {
    case ORD_HILBERT_THEN_LEX:
      possibleIdealsBasic.sort(LessByHilbert);
      break;
    case ORD_HILBERT_THEN_DEG:
      possibleIdealsBasic.sort(LessByHilbertThenDegree);
      break;
    case DEG_THEN_ORD_HILBERT:
      possibleIdealsBasic.sort(LessByDegreeThenHilbert);
      break;
    case GRAD_HILB_THEN_DEG:
      possibleIdealsBasic.sort(LessByGradHilbertThenDegree);
      break;
    case DEG_THEN_GRAD_HILB:
      possibleIdealsBasic.sort(LessByDegreeThenGradHilbert);
      break;
    case SMOOTHEST_DEGREES:
      possibleIdealsBasic.sort(LessBySmoothestDegrees);
      break;
    case LARGEST_MAX_COMPONENT:
      possibleIdealsBasic.sort(LessByLargestMaxComponent);
      break;
    case MIN_CRIT_PAIRS:
      possibleIdealsBasic.sort(LessByNumCritPairs);
      break;
    default: possibleIdealsBasic.sort(LessByHilbert);
  }
  /*cout << "sorted list\n";
  for (list<PPWithIdeal>::iterator winner = possibleIdealsBasic.begin(); winner != possibleIdealsBasic.end(); ++winner)
  { cout << '\t'; p_Write(winner->getPP(), Rx); }*/
  list<PPWithIdeal>::iterator winner = possibleIdealsBasic.begin();
  bool searching = true;
  if (possibleIdealsBasic.size() != 1)
  {
    // test each combination of LPPs for consistency
    // one of them must work (current LPP, if nothing else -- see previous case) 
    set<poly> PPunion;
    for (set<poly>::iterator piter = compatiblePPs.begin(); piter != compatiblePPs.end(); ++piter)
      PPunion.insert(*piter);
    for (set<poly>::iterator piter = boundaryPPs.begin(); piter != boundaryPPs.end(); ++piter)
      PPunion.insert(*piter);
    for (
         list<PPWithIdeal>::iterator siter = possibleIdealsBasic.begin();
         searching && siter != possibleIdealsBasic.end();
         ++siter
        )
    {
      skeleton newSkeleton(currSkel);
      vector<constraint> newvecs;
      //cout << "testing "; p_Write((*siter).getPP(),Rx);
      //ConstraintsForNewPP(*siter, compatiblePPs, newvecs);
      ConstraintsForNewPP(*siter, PPunion, newvecs);
      if (newSkeleton.ddm(newvecs))
      {
        //cout << "consistent\n";
        //if (verifyAndModifyIfNecessary(newSkeleton, CurrentPolys))
        if (verifyAndModifyIfNecessary(Tx, newSkeleton, CurrentPolys, numPolys))
        {
          searching = false;
          currSkel = newSkeleton;
          winner = siter;
        }
      }
      else
      {
        // cout << "inconsistent\n";
        // cout << newSkeleton;
        // this monomial is not, in fact, compatible
        compatiblePPs.erase(siter->getPP());
      }
    }
  }
  else if (possibleIdealsBasic.size() == 1 && compatiblePPs.size() != 1)
  {
    vector<constraint> newvecs;
    ConstraintsForNewPP(*(possibleIdealsBasic.begin()), compatiblePPs, newvecs);
    currSkel.ddm(newvecs);
    verifyAndModifyIfNecessary(Tx, currSkel, CurrentPolys, numPolys);
  }
    
  // set marked lpp
  //r->setMarkedLPP(winner->getPP());
  poly t = p_Copy_noCheck(winner->getPP(), Rx);
  t->next = NULL;
  CurrentLPPs.push_back(t);
  // TODO: delete elements of allPPs (not clear how: elements are in a set)
  ray new_weight = ray_sum(currSkel.get_rays());
  new_weight.simplify_ray();
  for (int i = 0; i < (int )(new_weight.get_dimension()); ++i)
    ordering_changed = ordering_changed or ((int )new_weight[i]) != Rx->wvhdl[0][i];
  cout << "ordering changed? " << ordering_changed << endl;
  //cout << "finished with "; pWrite(t);
  /*for (unsigned long i = 0; i < CurrentLPPs.size(); ++i) pWrite(CurrentLPPs[i]);
  cout << endl;
  cout << "skeleton after:\n";
  cout << currSkel; */
  cout << "returning from selmon\n";
}

#endif