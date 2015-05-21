#ifndef __DYNAMIC_ENGINE_H__
#define __DYNAMIC_ENGINE_H__

//#define NULL 0
#include <cstddef>

#include <kernel/combinatorics/hilb.h>
#include <kernel/combinatorics/stairc.h>
#include <kernel/ideals.h>
#include <kernel/GBEngine/kutil.h>
#include <polys/monomials/p_polys.h>

#include "skeleton.h"
#include <vector>
#include <iostream>

using namespace std;

extern ring currRing;

/**
  Pass this value to SelectMonomial. The methods are:
    -# `ORD_HILBERT_THEN_DEG` considers the standard Hilbert polynomial,
       breaking ties by the standard Hilbert series,
       breaking remaining ties by the (weighted) degree.
    -# `DEG_THEN_ORD_HILBERT` considers the degree of the leading power products,
       breaking ties by the standard Hilbert polynomial,
       breaking remaining ties by the standard Hilbert series.
    -# `GRAD_HILB_THEN_DEG` considers the standard Hilbert polynomial,
       breaking ties by the graded Hilbert series,
       breaking remaining ties by the (weighted) degree.
       (The graded Hilbert polynomial is actuall a quasi-polynomial,
       and we do not yet have a criterion to apply to this.)
    -# `DEG_THEN_GRAD_HILB` considers the degree of the leading power products,
       breaking ties by the standard Hilbert polynomial,
       breaking remaining ties by the graded Hilbert series.
    -# `SIGNATURE_HILBERT` considers the standard Hilbert polynomial,
       then the standard Hilbert series, to the embedded ideal
       of the products of the leading power products and their signatures.
    -# `SIGNATURE_GRAD_HILB` considers the standard Hilbert polynomial,
       then the graded Hilbert series, to the embedded ideal
       of the products of the leading power products and their signatures.
    -# `SIGNATURE_DEG` considers the (weighted) degree of the product
       of the leading power product and its signature.
    -# `MIN_CRIT_PAIRS` tries to minimize the number of new critical pairs.

   Ties that persist beyond these are broken lexicographically.
   \warning Most are not yet implemented!
*/
enum DynamicHeuristic
{
  ORD_HILBERT_THEN_LEX = 1, ORD_HILBERT_THEN_DEG, DEG_THEN_ORD_HILBERT,
  GRAD_HILB_THEN_DEG, DEG_THEN_GRAD_HILB,
  SIGNATURE_HILBERT, SIGNATURE_GRAD_HILB, SIGNATURE_DEG, MIN_CRIT_PAIRS
};

poly HilbertPolynomial(intvec *, int);

/**
  Used to associate a potential leading power product
  with the resulting monomial ideal if it were chosen
  as the actual leading power product.
*/
class PPWithIdeal
{
  public:
    PPWithIdeal(poly u, vector<poly> F, ray w)
    : t(p_Copy_noCheck(u,currRing)), ordering(w)
    {
      I = idInit(F.size() + 1);
      for (unsigned long i = 0; i < F.size(); ++i) { I->m[i] = F[i]; }
      I->m[F.size()] = u;
      hNum = NULL;
      hPol = NULL;
      // for series comparison, need the same denom, so first series
      //hNum = hFirstSeries(I, NULL, NULL); // third one to be modified for graded
      // for polynomial comparison, need the second series
      //hPol = HilbertPolynomial(hSecondSeries(hNum), scDimInt(I));
    };
    ~PPWithIdeal() { /*idDelete(&I, currRing);*/ }
    inline poly getPP() const { return t; };
    inline ideal getIdeal() const { return I; };
    inline ray getOrdering() const { return ordering; };
    inline intvec * getHilbertNumerator()
    {
      if (hNum != NULL) return hNum; else return hNum = hFirstSeries(I, NULL, NULL);
    };
    inline poly getHilbertPolynomial()
    {
      if (hPol != NULL) return hPol;
      else
        return hPol = HilbertPolynomial(hSecondSeries(getHilbertNumerator()),
                                                      scDimInt(I));
    };
  protected:
    poly t; ideal I; ray ordering; intvec * hNum; poly hPol;
};

/**
  Compute the compatible leading monomials of a polynomial.
  \param currentLPP the current leading power product of the polynomial
  \param allPPs set of all power products of the polynomial
  \param bndrys boundary (or, &ldquo;corner&rdquo;) vectors of current cone
  \param result set of power products of the polynomial compatible with `bndrys`
  \todo compare each potential PP with all other potential PPs,
    reducing the number of false positives
*/
void compatiblePP(
  poly currentLPP,            // the current LPP
  const set<poly> &allPPs,    // the monomials to consider;
                              // some will be removed
  const set<ray> &bndrys,     // known boundary vectors
  set<poly> &result           // returned as PPs for Hilbert function
                              // ("easy" (& efficient?) to extract exps
);

/**
  Verifies that the leading power products of the current basis
  remain compatible with the proposed refinement of ordering.
  \param skel the skeleton corresponding to the choice of \f$w\f$
  \param currentPolys the current basis
*/
bool verifyAndModifyIfNecessary(
  skeleton &skel,
  const vector<poly> &currentPolys
);

/**
  Create constraints for a candidate LPP.
  \param I pair of PP with the ideal it would have.
  \param monomialsForComparison monomials used to generate constraints with LPP
  \param result the new constraints
*/
void ConstraintsForNewPP(
  const PPWithIdeal &I,
  const set<poly> &monomialsForComparison,
  vector<constraint> &result
);

/**
  Selects a leading power product for a polynomial,
  applying a particular DynamicHeuristic (default is `ORD_HILBERT_THEN_DEG`)
  and ensuring compatibility with previous choices of LPP for other monomials.
  \param r the polynomial in need of a new choice of LPP
  \param CurrentLPPs the current choices of LPPs for `CurrentPolys`
  \param CurrentPolys the current basis of the ideal
  \param currSkel the current skeleton, corresponding to the choices `CurrentLPPs`
  \param method the method to apply; see `DynamicHeuristic`
*/
void SelectMonomial(
    poly &r,                          // changes
    vector<poly> &CurrentLPPs,      // changes
    //const vector<poly> &CurrentPolys,
    polyset CurrentPolys,
    int numPolys,
    skeleton & currSkel,                // possibly changes
    bool &ordering_changed,
    DynamicHeuristic method = ORD_HILBERT_THEN_DEG
);

#endif