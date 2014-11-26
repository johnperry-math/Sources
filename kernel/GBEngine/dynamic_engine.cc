#ifndef __DYNAMIC_ENGINE_C__
#define __DYNAMIC_ENGINE_C__

#include "dynamic_engine.h"

#include <coeffs/coeffs.h>
#include <polys/coeffrings.h>
#include <polys/monomials/p_polys.h>

#include <iostream>
using std::cout; using std::endl;

ring QQt = NULL;

struct LessByHilbert {
  bool operator () (const PPWithIdeal &a, const PPWithIdeal &b)
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
  return result;
  }
};

/*poly monomial(int a, int *e, ring R)
{ poly p = p_One(R); p_SetExpV(p, e, R); p_SetCoeff(p, n_Init(a, R), R); p_Setm(p, R); return p; }

poly monomial_univariate(number a, ring R)
{ poly p = p_One(R); p_SetExp(p, 1, 0, R); p_SetCoeff(p, n_Init(a, R), R); p_Setm(p, R); return p;}*/

poly var_univariate(ring R)
{ poly p = p_One(R); p_SetExp(p, 1, 1, R); p_Setm(p, R); return p; }

poly constant(number a, ring R)
{ poly p = p_One(R); p_SetCoeff(p, a, R); p_Setm(p, R); return p; }

char * tstring = "t";

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
  const set<ray> &bndrys,             // known boundary vectors
  set<poly> &result           // returned as PPs for Hilbert function
                                      // ("easy" (& efficient?) to extract exps
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
  result.insert(currentLPP);
  // compare other monomials with LPP
  for (
        set<poly>::iterator b_ptr = allPPs.begin();
        b_ptr != allPPs.end();
        ++b_ptr
      )
  {
    // take the dot product of each monomial's exponents with pp,
    // add the ones that pass muster to result
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
      if ((*w_ptr) * bray > (*w_ptr) * aray)
      {
        // only need one success
        result.insert(*b_ptr);
        searching = false;
      }
    }
  }
  delete [] a; delete [] b; delete [] along;
}

bool verifyAndModifyIfNecessary(
  skeleton &skel,
  const vector<poly> &currentPolys
)
{
  bool consistent = true; // innocent until proven guilty
  ray w = ray_sum(skel.get_rays()); // our tentative ordering
  unsigned long n = w.get_dimension();
  ring Rx = currRing;
  int * pExp = new int[n + 1];
  int * qExp = new int[n + 1];
  unsigned long long *entries = new unsigned long long [n]; // used for entries for new rays
  long *coefficients = new long [n]; // used for coefficients for new constraints
  // loop through all polynomials; verify leading power product is unchanged
  for (
       vector<const poly>::iterator piter = currentPolys.begin();
       piter != currentPolys.end();
       ++piter
      )
  {
    // create a ray for the current LPP's exponents
    poly pp = p_Head(*piter, Rx);
    p_GetExpV(pp, pExp, Rx);
    for (unsigned long i = 0; i < n; ++i)
      entries[i] = pExp[i + 1];
    ray a(n, entries);
    // loop through the polynomial's remaining monomials
    for (poly titer = (*piter)->next; consistent and titer != NULL; titer = titer->next)
    {
      if (pp != titer) // don't compare with LPP; that would be Bad (TM)
      {
        // create a ray for the PP's exponents
        p_GetExpV(titer, qExp, Rx);
        for (unsigned long i = 0; i < n; ++i)
          entries[i] = qExp[i + 1];
        ray b(n, entries);
        // compare weights between a and b; if this fails,
        // recompute the skeleton with a new constraint
        if (a*w < b*w)
        {
          if (coefficients == NULL) // ensure we have space
            coefficients = new long[n];
          for (unsigned long i = 0; i < n; ++i)
            coefficients[i] = a[i] - b[i];
          constraint new_constraint(n, coefficients);
          skeleton newskel(skel);
          consistent = newskel.ddm(new_constraint);
          // if we're consistent, we need to recompute the ordering
          if (consistent)
          {
            w = ray_sum(skel.get_rays());
            skel = newskel;
          } // if consistent
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
    poly &r,                          // changes
    vector<poly> &CurrentLPPs,      // changes
    const vector<poly> &CurrentPolys, // changes   
    skeleton & currSkel,                // possibly changes
    DynamicHeuristic method
)
{
  cout << "current ring: "; rWrite(currRing);
  cout << "entering selmon with "; pWrite(r); cout << endl;
  cout << "skeleton before: " << currSkel << endl;
  ring Rx = currRing;
  ray w = ray_sum(currSkel.get_rays());
  vector<unsigned long long> ord;
  for (unsigned long i = 0; i < w.get_dimension(); ++i) { ord.push_back(w[i]); }
  poly currentLPP = p_Head(r, Rx);
  set<poly> compatiblePPs, allPPs;
  // transform monomials into exponent vectors
  for (
       poly titer = r;
       titer != NULL;
       titer = titer->next
      )
  {
    poly t = pHead(titer);
    allPPs.insert(t);
    //allPPs.insert(p_Head(titer, Rx));
  }
  // loop through all exponent vectors
  compatiblePP(currentLPP, allPPs, currSkel.get_rays(), compatiblePPs);
  // list possible future ideals, sort by Hilbert Function
  // using a set sorts the ideals automagically by appropriate DynamicHeuristics
  set<PPWithIdeal, LessByHilbert> possibleIdealsBasic;
  //set<set<PPWithIdeal>, LessByGradedHilbert> possibleIdealsGraded;
  //set<set<PPWithIdeal>, LessByPseudoSig> possibleIdealsPseudoSig;
  //set<set<PPWithIdeal>, LessByVolume> possibleIdealsVolume;
  for (
       set<poly>::iterator piter = compatiblePPs.begin();
       piter != compatiblePPs.end();
       ++piter
      )
  {
    poly t = *piter;
    PPWithIdeal newIdeal(*piter, CurrentLPPs, w, r);
    possibleIdealsBasic.insert(newIdeal);
  }
  set<PPWithIdeal>::iterator winner = possibleIdealsBasic.begin();
  if (possibleIdealsBasic.size() != 1)
  {
    // test each combination of LPPs for consistency
    // one of them must work (current LPP, if nothing else -- see previous case) 
    bool searching = true;
    for (
         set<PPWithIdeal>::iterator siter = possibleIdealsBasic.begin();
         searching && siter != possibleIdealsBasic.end();
         ++siter
        )
    {
      skeleton newSkeleton(currSkel);
      vector<constraint> newvecs;
      ConstraintsForNewPP(*siter, compatiblePPs, newvecs);
      if (newSkeleton.ddm(newvecs))
      {
        if (verifyAndModifyIfNecessary(newSkeleton, CurrentPolys))
        {
          currSkel = newSkeleton;
          searching = false;
          winner = siter;
        }
      }
      else
        // this monomial is not, in fact, compatible
        compatiblePPs.erase(siter->getPP());
    }
  }
  // set marked lpp
  //r->setMarkedLPP(winner->getPP());
  poly t = p_Copy_noCheck(winner->getPP(), Rx);
  t->next = NULL;
  CurrentLPPs.push_back(t);
  // TODO: delete elements of allPPs (not clear how: elements are in a set
  cout << "finished with ";
  for (unsigned long i = 0; i < CurrentLPPs.size(); ++i) pWrite(CurrentLPPs[i]);
  cout << endl;
  cout << "skeleton after:\n";
  cout << currSkel;
  cout << "returning from selmon\n";
}

#endif