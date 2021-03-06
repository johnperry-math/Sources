// -*- Mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*-
/*****************************************************************************\
 * Computer Algebra System SINGULAR
\*****************************************************************************/
/** @file syzextra.cc
 *
 * New implementations for the computation of syzygies and resolutions
 *
 * ABSTRACT: Computation of Syzygies due to Schreyer
 *
 * @author Oleksandr Motsak
 *
 **/
/*****************************************************************************/

// include header file
#include <kernel/mod2.h>

#include "syzextra.h"

#include "DebugPrint.h"

#include <omalloc/omalloc.h>

#include <misc/intvec.h>
#include <misc/options.h>

#include <coeffs/coeffs.h>

#include <polys/monomials/p_polys.h>
#include <polys/monomials/ring.h>
#include <polys/simpleideals.h>

#include <polys/kbuckets.h> // for kBucket*
#include <polys/sbuckets.h> // for sBucket*
//#include <polys/nc/summator.h> // for CPolynomialSummator
#include <polys/operations/p_Mult_q.h> // for MIN_LENGTH_BUCKET

#include <kernel/GBEngine/kstd1.h>
#include <kernel/polys.h>
#include <kernel/GBEngine/syz.h>
#include <kernel/ideals.h>

#include <kernel/oswrapper/timer.h>


#include <Singular/tok.h>
#include <Singular/ipid.h>
#include <Singular/lists.h>
#include <Singular/attrib.h>

#include <Singular/ipid.h>
#include <Singular/ipshell.h> // For iiAddCproc

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define NOPRODUCT 1

// USING_NAMESPACE_SINGULARXX;
USING_NAMESPACE( SINGULARXXNAME :: DEBUG )


BEGIN_NAMESPACE_SINGULARXX     BEGIN_NAMESPACE(SYZEXTRA)


BEGIN_NAMESPACE_NONAME

static inline poly pp_Add_qq( const poly a, const poly b, const ring R)
{
  return p_Add_q( p_Copy(a, R), p_Copy(b, R), R );
}

static inline poly p_VectorProductLT( poly s,  const ideal& L, const ideal& T, const ring& R)
{
  assume( IDELEMS(L) == IDELEMS(T) );
  poly vp = NULL; // resulting vector product

  while( s != NULL )
  {
    const poly nxt = pNext(s);
    pNext(s) = NULL;

    if( !n_IsZero( p_GetCoeff(s, R), R) )
    {
      const int i = p_GetComp(s, R) - 1;
      assume( i >= 0 ); assume( i < IDELEMS(L) );
      p_SetComp(s, 0, R); p_SetmComp(s, R);

      vp = p_Add_q( vp, pp_Mult_qq( s, L->m[i], R ), R); 
      vp = p_Add_q( vp, pp_Mult_qq( s, T->m[i], R ), R); 
    }

    p_Delete(&s, R);

    s = nxt;
  };

  assume( s == NULL );

  return vp;
}

static inline int atGetInt(idhdl rootRingHdl, const char* attribute, long def)
{
  return ((int)(long)(atGet(rootRingHdl, attribute, INT_CMD, (void*)def)));
}

END_NAMESPACE

BEGIN_NAMESPACE(SORT_c_ds)

#if (defined(HAVE_QSORT_R) && (defined __APPLE__ || defined __MACH__ || defined __DARWIN__ || defined __FreeBSD__ || defined __BSD__ || defined OpenBSD3_1 || defined OpenBSD3_9))
static int cmp_c_ds(void *R, const void *p1, const void *p2){
#elif (defined(HAVE_QSORT_R) && (defined _GNU_SOURCE || defined __GNU__ || defined __linux__))
static int cmp_c_ds(const void *p1, const void *p2, void *R){
#else
static int cmp_c_ds(const void *p1, const void *p2){ void *R = currRing;
#endif
  assume(R != NULL);
  const int YES = 1;
  const int NO = -1;

  const ring r =  (const ring) R; // TODO/NOTE: the structure is known: C, lp!!!

  assume( r == currRing ); // for now...

  const poly a = *(const poly*)p1;
  const poly b = *(const poly*)p2;

  assume( a != NULL );
  assume( b != NULL );

  p_LmTest(a, r);
  p_LmTest(b, r);


  const signed long iCompDiff = p_GetComp(a, r) - p_GetComp(b, r);

  // TODO: test this!!!!!!!!!!!!!!!!

  //return -( compare (c, qsorts) )

#ifndef SING_NDEBUG
  const int __DEBUG__ = 0;
  if( __DEBUG__ )
  {
    PrintS("cmp_c_ds: a, b: \np1: "); dPrint(a, r, r, 0);
    PrintS("b: "); dPrint(b, r, r, 0);
    PrintLn();
  }
#endif


  if( iCompDiff > 0 )
    return YES;

  if( iCompDiff < 0 )
    return  NO;

  assume( iCompDiff == 0 );

  const signed long iDegDiff = p_Totaldegree(a, r) - p_Totaldegree(b, r);

  if( iDegDiff > 0 )
    return YES;

  if( iDegDiff < 0 )
    return  NO;

  assume( iDegDiff == 0 );

#ifndef SING_NDEBUG
  if( __DEBUG__ )
  {
    PrintS("cmp_c_ds: a & b have the same comp & deg! "); PrintLn();
  }
#endif

  for (int v = rVar(r); v > 0; v--)
  {
    assume( v > 0 );
    assume( v <= rVar(r) );

    const signed int d = p_GetExp(a, v, r) - p_GetExp(b, v, r);

    if( d > 0 )
      return YES;

    if( d < 0 )
      return NO;

    assume( d == 0 );
  }

  return 0;
}

/*
static int cmp_poly(const poly &a, const poly &b)
{
  const int YES = 1;
  const int NO = -1;

  const ring r =  (const ring) currRing; // TODO/NOTE: the structure is known: C, lp!!!

  assume( r == currRing );

  assume( a != NULL );
  assume( b != NULL );

  p_LmTest(a, r);
  p_LmTest(b, r);
  assume( p_GetComp(a, r) == 0 );
  assume( p_GetComp(b, r) == 0 );

#ifndef SING_NDEBUG
  const int __DEBUG__ = 0;
  if( __DEBUG__ )
  {
    PrintS("cmp_lex: a, b: \np1: "); dPrint(a, r, r, 0);
    PrintS("b: "); dPrint(b, r, r, 0);
    PrintLn();
  }
#endif

  for (int v = rVar(r); v > 0; v--)
  {
    assume( v > 0 );
    assume( v <= rVar(r) );

    const signed int d = p_GetExp(a, v, r) - p_GetExp(b, v, r);

    if( d > 0 )
      return YES;

    if( d < 0 )
      return NO;

    assume( d == 0 );
  }

  return 0;
}
*/

END_NAMESPACE
/* namespace SORT_c_ds */

/// writes a monomial (p),
/// uses form x*gen(.) if ko != coloumn number of p
static void writeLatexTerm(const poly t, const ring r, const bool bCurrSyz = true, const bool bLTonly = true)
{
  if( t == NULL )
  {
    PrintS(" 0 ");
    return;
  }

  assume( r != NULL );
  const coeffs C = r->cf; assume(C != NULL);

  poly p = t;
  BOOLEAN writePlus = FALSE;

  do {
  assume( p != NULL );

  // write coef...
  number& n = p_GetCoeff(p, r);

  n_Normalize(n, C);

  BOOLEAN writeMult = FALSE; ///< needs * before pp or module generator

  BOOLEAN writeOne = FALSE; ///< need to write something after '-'!

  if( n_IsZero(n, C) )
  {
    PrintS( writePlus? " + 0" : " 0 " ); writePlus = TRUE; writeMult = TRUE;
//    return; // yes?
  }

  if (n_IsMOne(n, C))
  {
    PrintS(" - "); writeOne = TRUE; writePlus = FALSE;
  }
  else if (!n_IsOne(n, C))
  {
    if( writePlus && n_GreaterZero(n, C) )
      PrintS(" + ");

    StringSetS(""); n_WriteLong(n, C);
    if (true)
    {
      char *s = StringEndS(); PrintS(s); omFree(s);
    }


    writeMult = TRUE;
    writePlus = TRUE;
  } else
     writeOne = TRUE;

  // missing '1' if no PP and gen...!?
  // write monom...
  const short N = rVar(r);

  BOOLEAN wrotePP = FALSE; ///< needs * before module generator?

  for (short i = 0; i < N; i++)
  {
    const long ee = p_GetExp(p, i+1, r);

    if (ee!=0L)
    {
      if (writeMult)
      {
        PrintS(" ");
        writeMult = FALSE;
      } else
      if( writePlus )
        PrintS(" + ");

      writePlus = FALSE;

      if (ee != 1L)
        Print(" %s^{%ld} ", rRingVar(i, r), ee);
      else
        Print(" %s ", rRingVar(i, r));

      writeOne = FALSE;
      wrotePP = TRUE;
    }
  }

  writePlus = writePlus || wrotePP;
  writeMult = writeMult || wrotePP;

  // write module gen...
  const long comp = p_GetComp(p, r);

  if (comp > 0 )
  {
    if (writeMult)
      PrintS("  ");
     else
      if( writePlus )
        PrintS(" + ");

    if (bCurrSyz)
      Print(" \\\\GEN{%ld} ", comp);
    else
      Print(" \\\\GENP{%ld} ", comp);

      writeOne = FALSE;
  }

  if ( writeOne )
    PrintS( writePlus? " + 1 " : " 1 " );


  pIter(p);

  writePlus = TRUE;
  } while( (!bLTonly) && (p != NULL) );

}





/// return a new term: leading coeff * leading monomial of p
/// with 0 leading component!
poly leadmonom(const poly p, const ring r, const bool bSetZeroComp)
{
  poly m = NULL;

  if( p != NULL )
  {
    assume( p != NULL );
    p_LmTest(p, r);

    m = p_LmInit(p, r);
    p_SetCoeff0(m, n_Copy(p_GetCoeff(p, r), r), r);

    if( bSetZeroComp )
      p_SetComp(m, 0, r);
    p_Setm(m, r);


    assume( m != NULL );
    assume( pNext(m) == NULL );
    p_LmTest(m, r);

    if( bSetZeroComp )
      assume( p_GetComp(m, r) == 0 );
  }

  return m;
}



poly p_Tail(const poly p, const ring r)
{
  if( p == NULL)
    return NULL;
  else
    return p_Copy( pNext(p), r );
}


ideal id_Tail(const ideal id, const ring r)
{
  if( id == NULL)
    return NULL;

  const ideal newid = idInit(IDELEMS(id),id->rank);

  for (int i=IDELEMS(id) - 1; i >= 0; i--)
    newid->m[i] = p_Tail( id->m[i], r );

  newid->rank = id_RankFreeModule(newid, currRing);

  return newid;
}



void Sort_c_ds(const ideal id, const ring r)
{
  const int sizeNew = IDELEMS(id);

#if ( (defined(HAVE_QSORT_R)) && (defined __APPLE__ || defined __MACH__ || defined __DARWIN__ || defined __FreeBSD__ || defined __BSD__ || defined OpenBSD3_1 || defined OpenBSD3_9) )
#define qsort_my(m, s, ss, r, cmp) qsort_r(m, s, ss, r, cmp)
#elif ( (defined(HAVE_QSORT_R)) && (defined _GNU_SOURCE || defined __GNU__ || defined __linux__))
#define qsort_my(m, s, ss, r, cmp) qsort_r(m, s, ss, cmp, r)
#else
#define qsort_my(m, s, ss, r, cmp) qsort(m, s, ss, cmp)
#endif

  if( sizeNew >= 2 )
    qsort_my(id->m, sizeNew, sizeof(poly), r, FROM_NAMESPACE(SORT_c_ds, cmp_c_ds));

#undef qsort_my

  id->rank = id_RankFreeModule(id, r);
}

/// Clean up all the accumulated data
void SchreyerSyzygyComputation::CleanUp()
{
  extern void id_Delete (ideal*, const ring);

  id_Delete(const_cast<ideal*>(&m_idTails), m_rBaseRing); // TODO!!!

  if( m_sum_bucket != NULL )
  {
    assume ( sIsEmpty(m_sum_bucket) );
    sBucketDestroy(&m_sum_bucket);
    m_sum_bucket = NULL;
  }

  if( m_spoly_bucket != NULL )
  {
    kBucketDestroy(&m_spoly_bucket);
    m_spoly_bucket = NULL;
  }

  for( TCache::iterator it = m_cache.begin(); it != m_cache.end(); it++ )
  {
    TP2PCache& T = it->second;

    for(TP2PCache::iterator vit = T.begin(); vit != T.end(); vit++ )
    {
      p_Delete( (&(vit->second)), m_rBaseRing);
      p_Delete( const_cast<poly*>(&(vit->first)), m_rBaseRing);
    }
  }
}
  /*
  for( TTailTerms::const_iterator it = m_idTailTerms.begin(); it != m_idTailTerms.end(); it++ )
  {
    const TTail& v = *it;
    for(TTail::const_iterator vit = v.begin(); vit != v.end(); vit++ )
      delete const_cast<CTailTerm*>(*vit);
  }
  */



int CReducerFinder::PreProcessTerm(const poly t, CReducerFinder& syzChecker) const
{
  assume( t != NULL );

  if( __DEBUG__ & __TAILREDSYZ__ )
    assume( !IsDivisible(t) ); // each input term should NOT be in <L>

  const ring r = m_rBaseRing;


  if( __TAILREDSYZ__ )
    if( p_LmIsConstant(t, r) ) // most basic case of baing coprime with L, whatever that is...
      return 1; // TODO: prove this...?

  //   return false; // appears to be fine

  const long comp = p_GetComp(t, r);

  CReducersHash::const_iterator itr = m_hash.find(comp);

  if ( itr == m_hash.end() )
    return 2; // no such leading component!!!

  assume( itr->first == comp );

  const bool bIdealCase = (comp == 0);
  const bool bSyzCheck = syzChecker.IsNonempty(); // need to check even in ideal case????? proof?  "&& !bIdealCase"

  if( __TAILREDSYZ__ && (bIdealCase || bSyzCheck) )
  {
    const TReducers& v = itr->second;
    const int N = rVar(r);
    // TODO: extract exps of t beforehand?!
    bool coprime = true;
    for(TReducers::const_iterator vit = v.begin(); (vit != v.end()) && coprime; ++vit )
    {
      assume( (*vit)->CheckLT( m_L ) );

      const poly p = (*vit)->lt();

      assume( p_GetComp(p, r) == comp );

      // TODO: check if coprime with Leads... if __TAILREDSYZ__ !
      for( int var = N; var > 0; --var )
        if( (p_GetExp(p, var, r) != 0) && (p_GetExp(t, var, r) != 0) )
        {
#ifndef SING_NDEBUG
          if( __DEBUG__ | 0)
          {
            PrintS("CReducerFinder::PreProcessTerm, 't' is NOT co-prime with the following leading term: \n");
            dPrint(p, r, r, 0);
          }
#endif
          coprime = false; // t not coprime with p!
          break;
        }

      if( bSyzCheck && coprime )
      {
        poly ss = p_LmInit(t, r);
        p_SetCoeff0(ss, n_Init(1, r), r); // for delete & printout only!...
        p_SetComp(ss, (*vit)->label() + 1, r); // coeff?
        p_Setm(ss, r);

        coprime = ( syzChecker.IsDivisible(ss) );

#ifndef SING_NDEBUG
        if( __DEBUG__ && !coprime)
        {
          PrintS("CReducerFinder::PreProcessTerm, 't' is co-prime with p but may lead to NOT divisible syz.term: \n");
          dPrint(ss, r, r, 0);
        }
#endif

        p_LmDelete(&ss, r); // deletes coeff as well???
      }

      assume( p == (*vit)->lt() );
      assume( (*vit)->CheckLT( m_L ) );
    }

#ifndef SING_NDEBUG
    if( __DEBUG__ && coprime )
      PrintS("CReducerFinder::PreProcessTerm, the following 't' is 'co-prime' with all of leading terms! \n");
#endif

    return coprime? 3: 0; // t was coprime with all of leading terms!!!

  }
  //   return true; // delete the term

  return 0;
}


void SchreyerSyzygyComputation::SetUpTailTerms()
{
  const ideal idTails = m_idTails;
  assume( idTails != NULL );
  assume( idTails->m != NULL );
  const ring r = m_rBaseRing;

  unsigned long pp[4] = {0,0,0,0}; // count preprocessed terms...
  
#ifndef SING_NDEBUG
  if( __DEBUG__ | 0)
  {
    PrintS("SchreyerSyzygyComputation::SetUpTailTerms(): Tails: \n");
    dPrint(idTails, r, r, 0);
  }
#endif

  for( int p = IDELEMS(idTails) - 1; p >= 0; --p )
    for( poly* tt = &(idTails->m[p]); (*tt) != NULL;  )
    {
      const poly t = *tt;
      const int k = m_div.PreProcessTerm(t, m_checker); // 0..3
      assume( 0 <= k && k <= 3 );

      if( __PROT__ )
        pp[k]++;

      if( k )
      {
#ifndef SING_NDEBUG
        if( __DEBUG__)
        {
          Print("SchreyerSyzygyComputation::SetUpTailTerms(): PP (%d) the following TT: \n", k);
          dPrint(t, r, r, 0);
        }
#endif

        (*tt) = p_LmDeleteAndNext(t, r); // delete the lead and next...
      }
      else
        tt = &pNext(t); // go next?

    }

#ifndef SING_NDEBUG
  if( !__TREEOUTPUT__ )
  if( __DEBUG__ | 0)
  {
    PrintS("SchreyerSyzygyComputation::SetUpTailTerms(): Preprocessed Tails: \n");
    dPrint(idTails, r, r, 0);
  }
#endif

  if( __PROT__ )
    Print("(PP/ST: {c: %lu, C: %lu, P: %lu} + %lu)", pp[1], pp[2], pp[3], pp[0]);

  
}
/*
  m_idTailTerms.resize( IDELEMS(idTails) );

  for( unsigned int p = IDELEMS(idTails) - 1; p >= 0; p -- )
  {
    TTail& v = m_idTailTerms[p];
    poly t = idTails->m[p];
    v.resize( pLength(t) );

    unsigned int pp = 0;

    while( t != NULL )
    {
      assume( t != NULL );
      // TODO: compute L:t!
//      ideal reducers;
//      CReducerFinder m_reducers

      CTailTerm* d = v[pp] = new CTailTerm();

      ++pp; pIter(t);
    }
  }
*/



ideal SchreyerSyzygyComputation::Compute1LeadingSyzygyTerms()
{
  const ideal& id = m_idLeads;
  const ring& r = m_rBaseRing;
//  const SchreyerSyzygyComputationFlags& attributes = m_atttributes;

  assume(!__LEAD2SYZ__);

  // 1. set of components S?
  // 2. for each component c from S: set of indices of leading terms
  // with this component?
  // 3. short exp. vectors for each leading term?

  const int size = IDELEMS(id);

  if( size < 2 )
  {
    const ideal newid = idInit(1, 0); newid->m[0] = NULL; // zero ideal...
    return newid;
  }

  // TODO/NOTE: input is supposed to be (reverse-) sorted wrt "(c,ds)"!??

  // components should come in groups: count elements in each group
  // && estimate the real size!!!


  // use just a vector instead???
  const ideal newid = idInit( (size * (size-1))/2, size); // maximal size: ideal case!

  int k = 0;

  for (int j = 0; j < size; j++)
  {
    const poly p = id->m[j];
    assume( p != NULL );
    const int  c = p_GetComp(p, r);

    for (int i = j - 1; i >= 0; i--)
    {
      const poly pp = id->m[i];
      assume( pp != NULL );
      const int  cc = p_GetComp(pp, r);

      if( c != cc )
        continue;

      const poly m = p_Init(r); // p_New???

      // m = LCM(p, pp) / p! // TODO: optimize: knowing the ring structure: (C/lp)!
      for (int v = rVar(r); v > 0; v--)
      {
        assume( v > 0 );
        assume( v <= rVar(r) );

        const short e1 = p_GetExp(p , v, r);
        const short e2 = p_GetExp(pp, v, r);

        if( e1 >= e2 )
          p_SetExp(m, v, 0, r);
        else
          p_SetExp(m, v, e2 - e1, r);

      }

      assume( (j > i) && (i >= 0) );

      p_SetComp(m, j + 1, r);
      pNext(m) = NULL;
      p_SetCoeff0(m, n_Init(1, r->cf), r); // for later...

      p_Setm(m, r); // should not do anything!!!

      newid->m[k++] = m;
    }
  }

//   if( __DEBUG__ && FALSE )
//   {
//     PrintS("ComputeLeadingSyzygyTerms::Temp0: \n");
//     dPrint(newid, r, r, 0);
//   }

  // the rest of newid is assumed to be zeroes...

  // simplify(newid, 2 + 32)??
  // sort(newid, "C,ds")[1]???
  id_DelDiv(newid, r); // #define SIMPL_LMDIV 32

//   if( __DEBUG__ && FALSE )
//   {
//     PrintS("ComputeLeadingSyzygyTerms::Temp1: \n");
//     dPrint(newid, r, r, 0);
//   }

  idSkipZeroes(newid); // #define SIMPL_NULL 2

//   if( __DEBUG__ )
//   {
//     PrintS("ComputeLeadingSyzygyTerms::Output: \n");
//     dPrint(newid, r, r, 0);
//   }

  Sort_c_ds(newid, r);

  return newid;
}

ideal SchreyerSyzygyComputation::Compute2LeadingSyzygyTerms()
{
  const ideal& id = m_idLeads;
  const ring& r = m_rBaseRing;
//  const SchreyerSyzygyComputationFlags& attributes = m_atttributes;

  // 1. set of components S?
  // 2. for each component c from S: set of indices of leading terms
  // with this component?
  // 3. short exp. vectors for each leading term?

  const int size = IDELEMS(id);

  if( size < 2 )
  {
    const ideal newid = idInit(1, 1); newid->m[0] = NULL; // zero module...
    return newid;
  }


  // TODO/NOTE: input is supposed to be sorted wrt "C,ds"!??

  // components should come in groups: count elements in each group
  // && estimate the real size!!!


  // use just a vector instead???
  ideal newid = idInit( (size * (size-1))/2, size); // maximal size: ideal case!

  int k = 0;

  for (int j = 0; j < size; j++)
  {
    const poly p = id->m[j];
    assume( p != NULL );
    const int  c = p_GetComp(p, r);

    for (int i = j - 1; i >= 0; i--)
    {
      const poly pp = id->m[i];
      assume( pp != NULL );
      const int  cc = p_GetComp(pp, r);

      if( c != cc )
        continue;

        // allocate memory & zero it out!
      const poly m = p_Init(r); const poly mm = p_Init(r);


        // m = LCM(p, pp) / p! mm = LCM(p, pp) / pp!
        // TODO: optimize: knowing the ring structure: (C/lp)!

      for (int v = rVar(r); v > 0; v--)
      {
        assume( v > 0 );
        assume( v <= rVar(r) );

        const short e1 = p_GetExp(p , v, r);
        const short e2 = p_GetExp(pp, v, r);

        if( e1 >= e2 )
          p_SetExp(mm, v, e1 - e2, r); //            p_SetExp(m, v, 0, r);
        else
          p_SetExp(m, v, e2 - e1, r); //            p_SetExp(mm, v, 0, r);

      }

      assume( (j > i) && (i >= 0) );

      p_SetComp(m, j + 1, r);
      p_SetComp(mm, i + 1, r);

      const number& lc1 = p_GetCoeff(p , r);
      const number& lc2 = p_GetCoeff(pp, r);

      number g = n_Lcm( lc1, lc2, r->cf );

      p_SetCoeff0(m ,       n_Div(g, lc1, r), r);
      p_SetCoeff0(mm, n_InpNeg(n_Div(g, lc2, r), r), r);

      n_Delete(&g, r);

      p_Setm(m, r); // should not do anything!!!
      p_Setm(mm, r); // should not do anything!!!

      pNext(m) = mm; //        pNext(mm) = NULL;

      newid->m[k++] = m;
    }
  }

//   if( __DEBUG__ && FALSE )
//   {
//     PrintS("Compute2LeadingSyzygyTerms::Temp0: \n");
//     dPrint(newid, r, r, 0);
//   }

  if( !__TAILREDSYZ__ )
  {
      // simplify(newid, 2 + 32)??
      // sort(newid, "C,ds")[1]???
    id_DelDiv(newid, r); // #define SIMPL_LMDIV 32

//     if( __DEBUG__ && FALSE )
//     {
//       PrintS("Compute2LeadingSyzygyTerms::Temp1 (deldiv): \n");
//       dPrint(newid, r, r, 0);
//     }
  }
  else
  {
      //      option(redSB); option(redTail);
      //      TEST_OPT_REDSB
      //      TEST_OPT_REDTAIL
    assume( r == currRing );

    BITSET _save_test; SI_SAVE_OPT1(_save_test);
    SI_RESTORE_OPT1(Sy_bit(OPT_REDTAIL) | Sy_bit(OPT_REDSB) | _save_test);

    intvec* w=new intvec(IDELEMS(newid));
    ideal tmp = kStd(newid, currRing->qideal, isHomog, &w);
    delete w;

    SI_RESTORE_OPT1(_save_test)

    id_Delete(&newid, r);
    newid = tmp;

//     if( __DEBUG__ && FALSE )
//     {
//       PrintS("Compute2LeadingSyzygyTerms::Temp1 (std): \n");
//       dPrint(newid, r, r, 0);
//     }

  }

  idSkipZeroes(newid);

  Sort_c_ds(newid, r);

  return newid;
}

poly SchreyerSyzygyComputation::TraverseNF(const poly a, const poly a2) const
{
#ifndef SING_NDEBUG
  if( __DEBUG__ )
  {
    m_div.Verify();
    m_checker.Verify();
  }
#endif

  const ideal& L = m_idLeads;
  const ideal& T = m_idTails;

  const ring& R = m_rBaseRing;

  const int r = p_GetComp(a, R) - 1;

  assume( r >= 0 && r < IDELEMS(T) );
  assume( r >= 0 && r < IDELEMS(L) );

  assume( a != NULL );

#ifndef SING_NDEBUG
  if( __DEBUG__ )
  {
    PrintS("SchreyerSyzygyComputation::TraverseNF(syz_lead, poly syz_2), \n");
    PrintS("syz_lead: \n");
    dPrint(a, R, R, 0);
    PrintS("syz_2: \n");
    dPrint(a2, R, R, 0);
    PrintLn();
  }
#endif
  
  if( __TREEOUTPUT__ )
  {
     PrintS("{ \"proc\": \"TraverseNF\", \"nodelabel\": \""); 
     writeLatexTerm(a, R); 
     PrintS("\", \"children\": [");
  }   

  poly aa = leadmonom(a, R); assume( aa != NULL); // :(

  poly t = TraverseTail(aa, r);

#ifndef SING_NDEBUG
  if( __DEBUG__ )
  {
    m_div.Verify();
    m_checker.Verify();
  }
#endif
  
  if( a2 != NULL )
  {
    assume( __LEAD2SYZ__ );

    if( __TREEOUTPUT__ )
    {
  
       PrintS("{ \"proc\": \"TraverseNF2\", \"nodelabel\": \""); 
       writeLatexTerm(a2, R); 
       PrintS("\", \"children\": [");
    }

    // replace the following... ?
    const int r2 = p_GetComp(a2, R) - 1; poly aa2 = leadmonom(a2, R); // :(

    assume( r2 >= 0 && r2 < IDELEMS(T) );

    t = p_Add_q(a2, p_Add_q(t, TraverseTail(aa2, r2), R), R);

    p_Delete(&aa2, R);


#ifndef SING_NDEBUG
    if( __DEBUG__ )
    {
      m_div.Verify();
      m_checker.Verify();
    }
#endif

    if( __TREEOUTPUT__ )
    {
       PrintS("], \"noderesult\": \"");
       writeLatexTerm(t, R, true, false);
       PrintS("\" },");
    }
     
    
  } else
    t = p_Add_q(t, ReduceTerm(aa, L->m[r], a), R);

  p_Delete(&aa, R);

  if( __TREEOUTPUT__ )
  {
     poly tt = pp_Add_qq( a, t, R);    
     PrintS("], \"noderesult\": \"");
     writeLatexTerm(tt, R, true, false);
     PrintS("\" },");
     p_Delete(&tt, R);
  }
#ifndef SING_NDEBUG
  if( __DEBUG__ )
  {
    PrintS("SchreyerSyzygyComputation::TraverseNF(syz_lead, poly syz_2), ==>>> \n");
    dPrint(t, R, R, 0);
    PrintLn();
  }
#endif

#ifndef SING_NDEBUG
  if( __DEBUG__ )
  {
    m_div.Verify();
    m_checker.Verify();
  }
#endif
  
  return t;
}

void SchreyerSyzygyComputation::ComputeSyzygy()
{
#ifndef SING_NDEBUG
  if( __DEBUG__ )
  {
    m_div.Verify();
    m_checker.Verify();
  }
#endif

  assume( m_idLeads != NULL );
  assume( m_idTails != NULL );

  const ideal& L = m_idLeads;
  const ideal& T = m_idTails;

  ideal& TT = m_syzTails;
  const ring& R = m_rBaseRing;

  if( m_sum_bucket == NULL )
    m_sum_bucket = sBucketCreate(R);

  assume ( sIsEmpty(m_sum_bucket) );

  if( m_spoly_bucket == NULL )
    m_spoly_bucket = kBucketCreate(R);


  assume( IDELEMS(L) == IDELEMS(T) );
#ifndef SING_NDEBUG
  int t, r;
#endif

  if( __TREEOUTPUT__ )
    Print("\n{ \"syzygylayer\": \"%d\", \"hybridnf\": \"%d\", \"diagrams\": \n[", __SYZNUMBER__, __HYBRIDNF__ );
  
  if( __PROT__ )
    Print("SYZ{%d}:", __SYZNUMBER__ );

  if( m_syzLeads == NULL )
  {
#ifndef SING_NDEBUG
    if( __PROT__ )
    {
      t = getTimer(); r = getRTimer();
      Print("\n%% %5d **!TIME4!** SchreyerSyzygyComputation::ComputeSyzygy::ComputeLeadingSyzygyTerms: t: %d, r: %d\n", r, t, r);
    }
#endif
    ComputeLeadingSyzygyTerms( __LEAD2SYZ__ && !__IGNORETAILS__ ); // 2 terms OR 1 term!
#ifndef SING_NDEBUG
    if( __PROT__ )
    {
      t = getTimer() - t; r = getRTimer() - r;
      Print("\n%% %5d **!TIME4!** SchreyerSyzygyComputation::ComputeSyzygy::ComputeLeadingSyzygyTerms: dt: %d, dr: %d\n", getRTimer(), t, r);
    }
#endif

  }

  assume( m_syzLeads != NULL );
  ideal& LL = m_syzLeads;
  const int size = IDELEMS(LL);

  TT = idInit(size, 0);

  if( size == 1 && LL->m[0] == NULL )
  {
     if( __TREEOUTPUT__ )
       PrintS("]},");
     return;
  }


  // use hybrid (Schreyer NF) method?
  const bool method = (__HYBRIDNF__  == 1); //  || (__HYBRIDNF__ == 2 && __SYZNUMBER__ < 3);

  if( __PROT__ )
     {
	if (method)
	  Print("(NF:PR, %s)", (NOPRODUCT == 1)? "*_*": "-*-" );
	else
	  Print("(NF:TT, %s)", (NOPRODUCT == 1)? "*_*": "-*-" );
     }
   
   

  if(  !__IGNORETAILS__)
  {
    if( T != NULL )
    {
#ifndef SING_NDEBUG
      if( __PROT__ )
      {
        t = getTimer(); r = getRTimer();
        Print("\n%% %5d **!TIME4!** SchreyerSyzygyComputation::ComputeSyzygy::SetUpTailTerms(): t: %d, r: %d\n", r, t, r);
      }
#endif

      SetUpTailTerms();
       
#ifndef SING_NDEBUG
      if( __PROT__ )
      {
        t = getTimer() - t; r = getRTimer() - r;
        Print("\n%% %5d **!TIME4!** SchreyerSyzygyComputation::ComputeSyzygy::SetUpTailTerms(): dt: %d, dr: %d\n", getRTimer(), t, r);
      }
#endif
    }
  }

#ifndef SING_NDEBUG
  if( __PROT__ )
  {
    t = getTimer(); r = getRTimer();
    Print("\n%% %5d **!TIME4!** SchreyerSyzygyComputation::ComputeSyzygy::SyzygyLift: t: %d, r: %d\n", r, t, r);
  }
#endif

#ifndef SING_NDEBUG
  if( __DEBUG__ )
  {
    m_div.Verify();
    m_checker.Verify();
  }
#endif
  
//  for( int k = 0; k < size; ++k ) // TODO: should be fine now!
  for( int k = size - 1; k >= 0; --k ) 
  {
    const poly a = LL->m[k]; assume( a != NULL );

    poly a2 = pNext(a);

    // Splitting 2-terms Leading syzygy module
    if( a2 != NULL )
      pNext(a) = NULL;

    if( __IGNORETAILS__ )
    {
      TT->m[k] = NULL;

      assume( a2 != NULL );

      if( a2 != NULL )
        p_Delete(&a2, R);

      continue;
    }

    //    TT->m[k] = a2;

#ifndef SING_NDEBUG
    if( __DEBUG__ )
    {
      m_div.Verify();
      m_checker.Verify();
    }
#endif

    poly nf;
    
    if( method )
      nf = SchreyerSyzygyNF(a, a2);
    else
      nf = TraverseNF(a, a2); // TODO: WRONG

#ifndef SING_NDEBUG
    if( __DEBUG__ )
    {
      m_div.Verify();
      m_checker.Verify();
    }
#endif

    TT->m[k] = nf; // ???
    
    if( __SYZCHECK__ )      
    {
      // TODO: check the correctness (syzygy property): a + TT->m[k] should be a syzygy!!!

      poly s = pp_Add_qq( a, TT->m[k], R); // current syzygy

      poly vp = p_VectorProductLT(s, L, T, R); 

      if( __DEBUG__ && (vp != NULL) && ! __TREEOUTPUT__ )
      {
        Warn("SchreyerSyzygyComputation::ComputeSyzygy: failed syzygy property for syzygy [%d], non-zero image is as follows: ", k);        
        dPrint(vp, R, R, 0);       p_Delete(&vp, R);
        
        PrintS("SchreyerSyzygyComputation::ComputeSyzygy: Wrong syzygy is as follows: ");
        s = pp_Add_qq( a, TT->m[k], R);
        dPrint(s, R, R, 0); p_Delete(&s, R);

        PrintS("SchreyerSyzygyComputation::ComputeSyzygy: Testing with the other method");
        
        if( !method )
          s = SchreyerSyzygyNF(a, a2);
        else
          s = TraverseNF(a, a2);

        s = p_Add_q( p_Copy(a, R), s, R); // another syzygy // :((((
        PrintS("SchreyerSyzygyComputation::ComputeSyzygy: The other method gives the following  syzygy: ");
        dPrint(s, R, R, 0); 
        
        vp = p_VectorProductLT(s, L, T, R);
        
        if( vp == NULL )
        {
          PrintS("SchreyerSyzygyComputation::ComputeSyzygy: .... which is correct!!! ");
        } else
        {
          Warn("SchreyerSyzygyComputation::ComputeSyzygy: failed to compute syzygy tail[%d] with both methods!!! Non-zero image (2nd) is as follows: ", k);        
          dPrint(vp, R, R, 0);
        }

#ifndef SING_NDEBUG
        if( __DEBUG__ )
        {
          m_div.Verify();
          m_checker.Verify(); 
        }
#endif
        
      } else
        assume( vp == NULL );

      if( __PROT__ && (vp != NULL) )
        Print("!%d!", k); // check k'th syzygy failed
      
      p_Delete(&vp, R);
    }

#ifndef SING_NDEBUG
    if( __DEBUG__ )
    {
      m_div.Verify();
      m_checker.Verify();
    }
#endif
  }

#ifndef SING_NDEBUG
    if( __PROT__ )
    {
      t = getTimer() - t; r = getRTimer() - r;
      Print("\n%% %5d **!TIME4!** SchreyerSyzygyComputation::ComputeSyzygy::SyzygyLift: dt: %d, dr: %d\n", getRTimer(), t, r);
    }
#endif

  TT->rank = id_RankFreeModule(TT, R);

  if( __TREEOUTPUT__ )
    PrintS("\n]},");
  
  if( __PROT__ )
    PrintLn();
}

void SchreyerSyzygyComputation::ComputeLeadingSyzygyTerms(bool bComputeSecondTerms)
{
//  const SchreyerSyzygyComputationFlags& attributes = m_atttributes;

//  const BOOLEAN __LEAD2SYZ__   = attributes.__LEAD2SYZ__;
//  const BOOLEAN __TAILREDSYZ__ = attributes.__TAILREDSYZ__;

  assume( m_syzLeads == NULL );

  if( bComputeSecondTerms )
  {
    assume( __LEAD2SYZ__ );
//    m_syzLeads = FROM_NAMESPACE(INTERNAL, _Compute2LeadingSyzygyTerms(m_idLeads, m_rBaseRing, m_atttributes));
    m_syzLeads = Compute2LeadingSyzygyTerms();
  }
  else
  {
    assume( !__LEAD2SYZ__ );

    m_syzLeads = Compute1LeadingSyzygyTerms();
  }
//    m_syzLeads = FROM_NAMESPACE(INTERNAL, _ComputeLeadingSyzygyTerms(m_idLeads, m_rBaseRing, m_atttributes));

  // NOTE: set m_LS if tails are to be reduced!
  assume( m_syzLeads!= NULL );

  if (__TAILREDSYZ__ && !__IGNORETAILS__ && (IDELEMS(m_syzLeads) > 0) && !((IDELEMS(m_syzLeads) == 1) && (m_syzLeads->m[0] == NULL)))
  {
    m_LS = m_syzLeads;
    m_checker.Initialize(m_syzLeads);
#ifndef SING_NDEBUG
    if( __DEBUG__ )
    {
      const ring& r = m_rBaseRing;
      PrintS("SchreyerSyzygyComputation::ComputeLeadingSyzygyTerms: \n");
      PrintS("m_syzLeads: \n");
      dPrint(m_syzLeads, r, r, 0);
      PrintS("m_checker.Initialize(m_syzLeads) => \n");
      m_checker.DebugPrint();
    }
#endif
    assume( m_checker.IsNonempty() ); // TODO: this always fails... BUG????
  }

  if( __PROT__ )
    Print("[L%dS:%d]", bComputeSecondTerms ? 2 : 1, IDELEMS(m_syzLeads));
      
}

poly SchreyerSyzygyComputation::SchreyerSyzygyNF(const poly syz_lead, poly syz_2) const
{
  assume( !__IGNORETAILS__ );

  const ideal& L = m_idLeads;
  const ideal& T = m_idTails;
  const ring& r = m_rBaseRing;

  assume( syz_lead != NULL );


#ifndef SING_NDEBUG
  if( __DEBUG__ )
  {
    PrintS("SchreyerSyzygyComputation::SchreyerSyzygyNF(syz_lead, poly syz_2), \n");
    PrintS("syz_lead: \n");
    dPrint(syz_lead, r, r, 0);
    PrintS("syz_2: \n");
    dPrint(syz_2, r, r, 0);
    PrintLn();
  }
#endif

  if( __TREEOUTPUT__ )
  {
    PrintS("{   \"nodelabel\": \""); writeLatexTerm(syz_lead, r);
    PrintS("\", \"children\": [");
  }
  
  if( syz_2 == NULL )
  {
    const int rr = p_GetComp(syz_lead, r) - 1;

    assume( rr >= 0 && rr < IDELEMS(T) );
    assume( rr >= 0 && rr < IDELEMS(L) );

#if NOPRODUCT
    syz_2 = m_div.FindReducer(syz_lead, L->m[rr], syz_lead, m_checker);

    if( __TREEOUTPUT__ )
    {
      PrintS("{ \"nodelabel\": \""); writeLatexTerm(syz_2, r); PrintS("\" },");
    }
#else
    poly aa = leadmonom(syz_lead, r); assume( aa != NULL); // :(
    aa = p_Mult_mm(aa, L->m[rr], r);

    if( __TREEOUTPUT__ )
    {
      PrintS("{ \"nodelabel\": \""); writeLatexTerm(syz_2, r); PrintS("\", \"edgelabel\": \""); writeLatexTerm(aa, r, false); PrintS("\" },");
    }

    syz_2 = m_div.FindReducer(aa, syz_lead, m_checker);

    p_Delete(&aa, r);
#endif

  }

  assume( syz_2 != NULL ); // by construction of S-Polynomial

  assume( L != NULL );
  assume( T != NULL );

  assume( IDELEMS(L) == IDELEMS(T) );

  int  c = p_GetComp(syz_lead, r) - 1;

  assume( c >= 0 && c < IDELEMS(T) );

  if( m_sum_bucket == NULL )
    m_sum_bucket = sBucketCreate(r);

  assume( sIsEmpty(m_sum_bucket) );
  
  if( m_spoly_bucket == NULL )
    m_spoly_bucket = kBucketCreate(r);

  sBucket_pt tail   = m_sum_bucket; assume( tail != NULL ); m_sum_bucket = NULL;
  
  kBucket_pt bucket = m_spoly_bucket; assume( bucket != NULL ); kbTest(bucket); m_spoly_bucket = NULL;
  
//  kBucketInit(bucket, NULL, 0); // not needed!?

  poly p = leadmonom(syz_lead, r); // :(
//  poly spoly = pp_Mult_qq(p, T->m[c], r);
  kBucket_Plus_mm_Mult_pp(bucket, p, T->m[c], 0); // TODO: store pLength(T->m[c]) separately!?
  p_Delete(&p, r);

  kbTest(bucket);

  c = p_GetComp(syz_2, r) - 1;
  assume( c >= 0 && c < IDELEMS(T) );

  p = leadmonom(syz_2, r); // :(
//  spoly = p_Add_q(spoly, pp_Mult_qq(p, T->m[c], r), r);
  kBucket_Plus_mm_Mult_pp(bucket, p, T->m[c], 0); // pLength(T->m[c])?!
  kbTest(bucket);
  p_Delete(&p, r);

  // TODO: use bucket!?
//  const bool bUsePolynomial = TEST_OPT_NOT_BUCKETS; //  || (pLength(spoly) < MIN_LENGTH_BUCKET);
//  CPolynomialSummator tail(r, bUsePolynomial);
  sBucket_Add_p(tail, syz_2, 1); // tail.AddAndDelete(syz_2, 1);

  kbTest(bucket);
  for( poly spoly = kBucketExtractLm(bucket); spoly != NULL; p_LmDelete(&spoly, r), spoly = kBucketExtractLm(bucket))
  {
    kbTest(bucket);
    poly t = m_div.FindReducer(spoly, NULL, m_checker);

    if( t != NULL )
    {
      p = leadmonom(t, r); // :(
      c = p_GetComp(t, r) - 1;

      assume( c >= 0 && c < IDELEMS(T) );

      if( __TREEOUTPUT__ )
      {
        PrintS("{ \"nodelabel\": \""); writeLatexTerm(t, r); PrintS("\", \"edgelabel\": \""); writeLatexTerm(spoly, r, false); PrintS("\" },");
      }

      kBucket_Plus_mm_Mult_pp(bucket, p, T->m[c], 0); // pLength(T->m[c])?
//      spoly = p_Add_q(spoly, pp_Mult_qq(p, T->m[c], r), r);

      p_Delete(&p, r);

      sBucket_Add_p(tail, t, 1); // tail.AddAndDelete(t, 1);
    } // otherwise discard that leading term altogether!
    else
      if( __PROT__ )
        PrintS("$"); // LOT
    
    kbTest(bucket);
  }

  // now bucket must be empty!
  kbTest(bucket);
  
  assume( kBucketClear(bucket) == NULL );
  
  

  poly result; int len;
  sBucketClearAdd(tail, &result, &len); // TODO: use Merge with sBucket???
  assume( pLength(result) == len );

  assume( tail != NULL ); assume ( sIsEmpty(tail) );

  if( m_sum_bucket == NULL )
    m_sum_bucket = tail;
  else
    sBucketDestroy(&tail);
  

  if( m_spoly_bucket == NULL )
    m_spoly_bucket = bucket;
  else
    kBucketDestroy(&bucket);
  
  
  if( __TREEOUTPUT__ )
  {
    PrintS("]},");
  }

#ifndef SING_NDEBUG
  if( __DEBUG__ )
  {
    PrintS("SchreyerSyzygyComputation::SchreyerSyzygyNF(syz_lead, poly syz_2) =>>> \n");
    dPrint(result, r, r, 0);
    PrintLn();
    // TODO: Add SyzCheck!!!???
  }
#endif

  return result;
}

// namespace     {

// };

bool my_p_LmCmp (poly a, poly b, const ring r) { return p_LmCmp(a, b, r) == -1; } // TODO: change to simple lex. memory compare!

// NOTE: need p_Copy?????? for image + multiplier!!???
// NOTE: better store complete syz. terms!!?
poly SchreyerSyzygyComputation::TraverseTail(poly multiplier, const int tail) const
{
#ifndef SING_NDEBUG
  if( __DEBUG__ )
  {
    m_div.Verify();
    m_checker.Verify();
  }
#endif

  const ring& r = m_rBaseRing;

  assume(m_idTails !=  NULL && m_idTails->m != NULL);
  assume( tail >= 0 && tail < IDELEMS(m_idTails) );


  

/*  return ComputeImage(multiplier, tail); */

  // TODO: store (multiplier, tail) -.-^-.-^-.--> !
  TCache::iterator top_itr = m_cache.find(tail);

  if ( top_itr != m_cache.end() )
  {
     assume( top_itr->first == tail );

     TP2PCache& T = top_itr->second;

     TP2PCache::iterator itr = T.find(multiplier);

     if( itr != T.end() ) // Yey - Reuse!!!
     {
       assume( p_LmEqual(itr->first, multiplier, r) );

       if( itr->second == NULL )
         return (NULL);

       if( __TREEOUTPUT__ )
       {
//         PrintS("{ \"nodelabel\": \""); writeLatexTerm(multiplier, r, false);
//         Print("  \\\\GEN{%d}\", \"children\": [ ", tail + 1);
         PrintS("{ \"proc\": \"TTLookup\", \"nodelabel\": \"");
         writeLatexTerm(itr->first, r, false); Print(" \\\\GEN{%d}\", \"Lookup\": \"", tail + 1);
         writeLatexTerm(itr->second, r, true, false);
         PrintS("\", ");
       }
       
       poly p = p_Copy(itr->second, r); // COPY!!!


       if( !n_Equal( pGetCoeff(multiplier), pGetCoeff(itr->first), r) ) // normalize coeffs!?
       {
         number n = n_Div( pGetCoeff(multiplier), pGetCoeff(itr->first), r); // new number
         
         if( __TREEOUTPUT__ )
         {
           StringSetS("");
           n_Write(n, r);
           char* s = StringEndS();
           Print("\"recale\": \"%s\", ", s);
           omFree(s);
         }
         
         if( __PROT__ )
           PrintS("l*"); // lookup & rescale
           
         
         p = p_Mult_nn(p, n, r); // !
         n_Delete(&n, r);        
       } else
         if( __PROT__ )
           PrintS("l"); // lookup no rescale
       

       if( __TREEOUTPUT__ )
       {
         PrintS("\"noderesult\": \"");         writeLatexTerm(p, r, true, false);         PrintS("\" },");
       }

#ifndef SING_NDEBUG
       if( __DEBUG__ )       {         m_div.Verify();         m_checker.Verify();       }
#endif
       
       return p;
     }


     if( __TREEOUTPUT__ )
     {
       Print("{ \"proc\": \"TTStore%d\", \"nodelabel\": \"", tail + 1); writeLatexTerm(multiplier, r, false); Print(" \\\\GEN{%d}\", \"children\": [", tail + 1);
     }
     
     const poly p = ComputeImage(multiplier, tail);

     if( __TREEOUTPUT__ )
     {
       PrintS("], \"noderesult\": \""); writeLatexTerm(p, r, true, false); PrintS("\" },");
     }

     if( __PROT__ )
       PrintS("S"); // store

     T.insert( TP2PCache::value_type(p_Copy(multiplier, r), p) ); //     T[ multiplier ] = p;

//     if( p == NULL )
//        return (NULL);

#ifndef SING_NDEBUG
     if( __DEBUG__ )
     {
       m_div.Verify();
       m_checker.Verify();
     }
#endif
     
     return p_Copy(p, r);
  }

  CCacheCompare o(r); TP2PCache T(o);

  if( __TREEOUTPUT__ )
  {
    Print("{ \"proc\": \"TTStore%d\", \"nodelabel\": \"", 0); writeLatexTerm(multiplier, r, false); Print(" \\\\GEN{%d}\", \"children\": [", tail + 1);
  }

  const poly p = ComputeImage(multiplier, tail);

  if( __TREEOUTPUT__ )
  {
    PrintS("], \"noderesult\": \""); writeLatexTerm(p, r, true, false); PrintS("\" },");
  }

  if( __PROT__ )
    PrintS("S"); // %d", tail + 1);
  
  T.insert( TP2PCache::value_type(p_Copy(multiplier, r), p) );

  m_cache.insert( TCache::value_type(tail, T) );

//  if( p == NULL )
//    return (NULL);

#ifndef SING_NDEBUG
  if( __DEBUG__ )
  {
    m_div.Verify();
    m_checker.Verify();
  }
#endif

  return p_Copy(p, r);
}

poly SchreyerSyzygyComputation::ComputeImage(poly multiplier, const int tail) const
{
  const ring& r = m_rBaseRing;

  assume(m_idTails !=  NULL && m_idTails->m != NULL);
  assume( tail >= 0 && tail < IDELEMS(m_idTails) );
  
  const poly t = m_idTails->m[tail]; // !!!

  if(t != NULL)
  {
    if( __TREEOUTPUT__ )
    {
      PrintS("{ \"proc\": \"ComputeImage\", \"nodelabel\": \"");
      writeLatexTerm(multiplier, r, false);
      Print(" \\\\GEN{%d}\", \"edgelabel\": \"", tail + 1);
      writeLatexTerm(t, r, false);
      PrintS("\", \"children\": [");
    }
    
    const poly p = TraverseTail(multiplier, t);

    if( __TREEOUTPUT__ )
    {
      PrintS("], \"noderesult\": \""); writeLatexTerm(p, r, true, false); PrintS("\" },");
    }
    
    return p;
    
  }

  return NULL;
}


poly SchreyerSyzygyComputation::TraverseTail(poly multiplier, poly tail) const
{
  assume( !__IGNORETAILS__ );

  const ideal& L = m_idLeads;
  const ideal& T = m_idTails;
  const ring& r = m_rBaseRing;

  assume( multiplier != NULL );

  assume( L != NULL );
  assume( T != NULL );

#ifndef SING_NDEBUG
  if( __DEBUG__ )
  {
    m_div.Verify();
    m_checker.Verify();
  }
#endif

  if(!(  (!__TAILREDSYZ__)   ||   m_lcm.Check(multiplier)     ))
  {
    if( __TAILREDSYZ__ && __PROT__ )
      PrintS("%"); // check LCM !
    
    return NULL;
  }
    
  //    const bool bUsePolynomial = TEST_OPT_NOT_BUCKETS; //  || (pLength(tail) < MIN_LENGTH_BUCKET);

  sBucket_pt sum;

  if( m_sum_bucket == NULL )
    sum = sBucketCreate(r);
  else
  {
    if( !sIsEmpty(m_sum_bucket) )
      sum = sBucketCreate(r);
    else
    {
      sum = m_sum_bucket;
      m_sum_bucket = NULL;
    }
  }


  assume( sum != NULL ); assume ( sIsEmpty(sum) );
  assume( r == sBucketGetRing(sum) );

  poly s; int len;

  //    CPolynomialSummator sum(r, bUsePolynomial);
  //    poly s = NULL;

  if( __TREEOUTPUT__ )
  {
    Print("{ \"proc\": \"TTPoly\", \"nodelabel\": \""); writeLatexTerm(multiplier, r, false); Print(" * \\\\ldots \", \"children\": [");
  }

  for(poly p = tail; p != NULL; p = pNext(p))   // iterate over the tail
  {
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    const poly rt = ReduceTerm(multiplier, p, NULL); // TODO: also return/store length?
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    const int lp = pLength(rt);

    if( rt != NULL && lp != 0 )
      sBucket_Add_p(sum, rt, lp);
  }

  sBucketClearAdd(sum, &s, &len); // Will Not Clear?!?

  assume( sum != NULL ); assume ( sIsEmpty(sum) );

  if( m_sum_bucket == NULL )
    m_sum_bucket = sum;
  else
    sBucketDestroy(&sum);   

  assume( pLength(s) == len );
#ifndef SING_NDEBUG
  if( __DEBUG__ )
  {
    m_div.Verify();
    m_checker.Verify();
  }
#endif

  if( __TREEOUTPUT__ )
  {
    PrintS("], \"noderesult\": \""); writeLatexTerm(s, r, true, false); PrintS("\" },");
  }

  return s;

}




poly SchreyerSyzygyComputation::ReduceTerm(poly multiplier, poly term4reduction, poly syztermCheck) const
{
#ifndef SING_NDEBUG
  if( __DEBUG__ )
  {
    m_div.Verify();
    m_checker.Verify();
  }
#endif

  assume( !__IGNORETAILS__ );

  const ideal& L = m_idLeads;
  const ideal& T = m_idTails;
  const ring& r = m_rBaseRing;

  assume( multiplier != NULL );
  assume( term4reduction != NULL );


  assume( L != NULL );
  assume( T != NULL );

  // simple implementation with FindReducer:
  poly s = NULL;

  if( (!__TAILREDSYZ__) || m_lcm.Check(multiplier) )
  {
#if NOPRODUCT
    s = m_div.FindReducer(multiplier, term4reduction, syztermCheck, m_checker);
    
    if( s == NULL ) // No Reducer?
    {      
      if( __PROT__ )
        PrintS("$" ); // LOT?!
      
      return s;
    }

    if( __TREEOUTPUT__ )
    {
      poly product = pp_Mult_mm(multiplier, term4reduction, r);
      PrintS("{ \"proc\": \"RdTrmNoP\", \"nodelabel\": \""); writeLatexTerm(s, r); PrintS("\", \"edgelabel\": \""); writeLatexTerm(product, r, false);
      p_Delete(&product, r);
    }

#else
    // NOTE: only LT(term4reduction) should be used in the following:
    poly product = pp_Mult_mm(multiplier, term4reduction, r);
    s = m_div.FindReducer(product, syztermCheck, m_checker); // ??

    if( s == NULL ) // No Reducer?
    {
      if( __PROT__ )
        PrintS("$" ); // LOT?!

      return s;
    }

    if( __TREEOUTPUT__ )
    {
      PrintS("{ \"proc\": \"RdTrmP\", \"nodelabel\": \""); writeLatexTerm(s, r); PrintS("\", \"edgelabel\": \""); writeLatexTerm(product, r, false);
    }

    p_Delete(&product, r);
#endif
  }

#ifndef SING_NDEBUG
  if( __DEBUG__ )
  {
    m_div.Verify();
    m_checker.Verify();
  }
#endif

  if( s == NULL ) // No Reducer?
  {
    if( __TAILREDSYZ__&& __PROT__ )
      PrintS("%" ); // LCM check
    
    return s;
  }
  
  poly b = leadmonom(s, r);

  const int c = p_GetComp(s, r) - 1;
  assume( c >= 0 && c < IDELEMS(T) );


  if( __TREEOUTPUT__ )
     PrintS("\", \"children\": [");

  const poly t = TraverseTail(b, c); // T->m[c];

  if( t != NULL )
    s = p_Add_q(s, t, r);

  if( __TREEOUTPUT__ )
  {

    PrintS("], \"noderesult\": \"");
    writeLatexTerm(s, r, true, false);
    PrintS("\"");

    if( syztermCheck != NULL )
    {
      PrintS(", \"syztermCheck\":\"" );
      writeLatexTerm(syztermCheck, r, true, false);
      PrintS("\" },");      
    } else
      PrintS(" },");
  }

  
#ifndef SING_NDEBUG
  if( __DEBUG__ )
  {
    m_div.Verify();
    m_checker.Verify();
  }
#endif

  
  return s;
}

SchreyerSyzygyComputationFlags::SchreyerSyzygyComputationFlags(idhdl rootRingHdl):
    __DEBUG__( atGetInt(rootRingHdl,"DEBUG", 0) ),
    __LEAD2SYZ__( atGetInt(rootRingHdl, "LEAD2SYZ", 0) ),
    __TAILREDSYZ__( atGetInt(rootRingHdl, "TAILREDSYZ", 1) ),
    __HYBRIDNF__( atGetInt(rootRingHdl, "HYBRIDNF", 0) ),
    __IGNORETAILS__( atGetInt(rootRingHdl, "IGNORETAILS", 0) ),
    __SYZNUMBER__( atGetInt(rootRingHdl, "SYZNUMBER", 0) ),
    __TREEOUTPUT__( atGetInt(rootRingHdl, "TREEOUTPUT", 0) ),
    __SYZCHECK__( atGetInt(rootRingHdl, "SYZCHECK", 0) ),
    __PROT__(TEST_OPT_PROT),
    m_rBaseRing( rootRingHdl->data.uring )
{
#ifndef SING_NDEBUG
  if( __DEBUG__ && 0 )
  {
    PrintS("SchreyerSyzygyComputationFlags: \n");
    Print("        DEBUG: \t%d\n", __DEBUG__);
//    Print("   SYZCHECK  : \t%d\n", __SYZCHECK__);
    Print("     LEAD2SYZ: \t%d\n", __LEAD2SYZ__);
    Print("   TAILREDSYZ: \t%d\n", __TAILREDSYZ__);
    Print("  IGNORETAILS: \t%d\n", __IGNORETAILS__);
    Print("   TREEOUTPUT: \t%d\n", __TREEOUTPUT__);
    Print("     SYZCHECK: \t%d\n", __SYZCHECK__);
  }
#endif

  // TODO: just current setting!
  assume( rootRingHdl == currRingHdl );
  assume( rootRingHdl->typ == RING_CMD );
  assume( m_rBaseRing == currRing );
  // move the global ring here inside???
}



CLeadingTerm::CLeadingTerm(unsigned int _label,  const poly _lt, const ring R):
    m_sev( p_GetShortExpVector(_lt, R) ),  m_label( _label ),  m_lt( _lt )
#ifndef SING_NDEBUG
    , _R(R), m_lt_copy( p_Head(_lt, R) )
#endif  
{
#ifndef SING_NDEBUG
  assume( pNext(m_lt_copy) == NULL );
#endif  
  assume( sev() == p_GetShortExpVector(lt(), R) );
}

#ifndef SING_NDEBUG
CLeadingTerm::~CLeadingTerm()
{
  assume( p_LmEqual(m_lt, m_lt_copy, _R) );
  assume( m_sev == p_GetShortExpVector(m_lt, _R) );
  
  poly p = const_cast<poly>(m_lt_copy);
  p_Delete(&p, _R);
}
poly CLeadingTerm::lt() const
{
  assume( p_LmEqual(m_lt, m_lt_copy, _R) );
  assume( m_sev == p_GetShortExpVector(m_lt, _R) );  
  return m_lt;
}

unsigned long CLeadingTerm::sev() const
{
  assume( p_LmEqual(m_lt, m_lt_copy, _R) );
  assume( m_sev == p_GetShortExpVector(m_lt, _R) );  
  return m_sev;
}

unsigned int CLeadingTerm::label() const
{
  assume( p_LmEqual(m_lt, m_lt_copy, _R) );
  assume( m_sev == p_GetShortExpVector(m_lt, _R) );  
  return m_label;
}
#endif  



CReducerFinder::~CReducerFinder()
{
  for( CReducersHash::const_iterator it = m_hash.begin(); it != m_hash.end(); it++ )
  {
    const TReducers& v = it->second;
    for(TReducers::const_iterator vit = v.begin(); vit != v.end(); vit++ )
      delete const_cast<CLeadingTerm*>(*vit);
  }
}


void CReducerFinder::Initialize(const ideal L)
{
  assume( m_L == NULL || m_L == L );
  if( m_L == NULL )
    m_L = L;

  assume( m_L == L );

  if( L != NULL )
  {
    const ring& R = m_rBaseRing;
    assume( R != NULL );

    for( int k = IDELEMS(L) - 1; k >= 0; k-- )
    {
      const poly a = L->m[k]; // assume( a != NULL );

      // NOTE: label is k \in 0 ... |L|-1!!!
      if( a != NULL )
        m_hash[p_GetComp(a, R)].push_back( new CLeadingTerm(k, a, R) );
    }
  }
}

CReducerFinder::CReducerFinder(const ideal L, const SchreyerSyzygyComputationFlags& flags):
    SchreyerSyzygyComputationFlags(flags),
    m_L(const_cast<ideal>(L)), // for debug anyway
    m_hash()
{
  assume( flags.m_rBaseRing == m_rBaseRing );
  if( L != NULL )
    Initialize(L);
}

/// _p_LmDivisibleByNoComp for a | b*c
static inline BOOLEAN _p_LmDivisibleByNoComp(const poly a, const poly b, const poly c, const ring r)
{
  int i=r->VarL_Size - 1;
  unsigned long divmask = r->divmask;
  unsigned long la, lb;

  if (r->VarL_LowIndex >= 0)
  {
    i += r->VarL_LowIndex;
    do
    {
      la = a->exp[i];
      lb = b->exp[i] + c->exp[i];
      if ((la > lb) ||
          (((la & divmask) ^ (lb & divmask)) != ((lb - la) & divmask)))
      {
        pDivAssume(p_DebugLmDivisibleByNoComp(a, b, r) == FALSE);
        return FALSE;
      }
      i--;
    }
    while (i>=r->VarL_LowIndex);
  }
  else
  {
    do
    {
      la = a->exp[r->VarL_Offset[i]];
      lb = b->exp[r->VarL_Offset[i]] + c->exp[r->VarL_Offset[i]];
      if ((la > lb) ||
          (((la & divmask) ^ (lb & divmask)) != ((lb - la) & divmask)))
      {
        pDivAssume(p_DebugLmDivisibleByNoComp(a, b, r) == FALSE);
        return FALSE;
      }
      i--;
    }
    while (i>=0);
  }
#ifdef HAVE_RINGS
  assume( !rField_is_Ring(r) ); // not implemented for rings...!
#endif
  return TRUE;
}


bool CLeadingTerm::CheckLT( const ideal & L ) const
{
//  for( int i = IDELEMS(L); i >= 0; --i) assume( pNext(L->m[i]) == NULL ); // ???
  return ( L->m[label()] == lt() );
} 

bool CLeadingTerm::DivisibilityCheck(const poly product, const unsigned long not_sev, const ring r) const
{
  // may have no coeff yet
//  assume ( !n_IsZero( p_GetCoeff(product, r), r ) );

  assume ( !n_IsZero( p_GetCoeff(lt(), r), r ) );
  assume( sev() == p_GetShortExpVector(lt(), r) );

  assume( product != NULL );
  assume( (p_GetComp(lt(), r) == p_GetComp(product, r)) || (p_GetComp(lt(), r) == 0) ); 

#ifndef SING_NDEBUG
  assume( r == _R );
#endif   
  
//  const int k = label();
//  assume( m_L->m[k] == p );

  return p_LmShortDivisibleByNoComp(lt(), sev(), product, not_sev, r);

}

/// as DivisibilityCheck(multiplier * t, ...) for monomial 'm'
/// and a module term 't'
bool CLeadingTerm::DivisibilityCheck(const poly m, const poly t, const unsigned long not_sev, const ring r) const
{
  assume ( !n_IsZero( p_GetCoeff(lt(), r), r ) );
  assume( sev() == p_GetShortExpVector(lt(), r) );

  assume( m != NULL );
  assume( t != NULL );  
  assume ( !n_IsZero( p_GetCoeff(m, r), r ) );
  assume ( !n_IsZero( p_GetCoeff(t, r), r ) );

// assume( p_GetComp(m, r) == 0 );
  assume( (p_GetComp(lt(), r) == p_GetComp(t, r))  || (p_GetComp(lt(), r) == 0)  );

//  const int k = label();
//  assume( m_L->m[k] == p );

#ifndef SING_NDEBUG
  assume( r == _R );
#endif   
  
  if (sev() & not_sev)
    return false;

  return _p_LmDivisibleByNoComp(lt(), m, t, r);
//  return p_LmShortDivisibleByNoComp(p, p_sev, product, not_sev, r);
}



/// TODO:
class CDivisorEnumerator: public SchreyerSyzygyComputationFlags
{
  private:
    const CReducerFinder& m_reds;
    const poly m_product;
    const unsigned long m_not_sev;
    const long m_comp;

    CReducerFinder::CReducersHash::const_iterator m_itr;
    CReducerFinder::TReducers::const_iterator m_current, m_finish;

    bool m_active;

  public:
    CDivisorEnumerator(const CReducerFinder& self, const poly product):
        SchreyerSyzygyComputationFlags(self),
        m_reds(self),
        m_product(product),
        m_not_sev(~p_GetShortExpVector(product, m_rBaseRing)),
        m_comp(p_GetComp(product, m_rBaseRing)),
        m_itr(), m_current(), m_finish(),
        m_active(false)
    {
      assume( m_comp >= 0 );
      assume( m_reds.m_L != NULL ); /// TODO: m_L should stay the same!!!
      
      assume( product != NULL ); // may have no coeff yet :(
//      assume ( !n_IsZero( p_GetCoeff(product, m_rBaseRing), m_rBaseRing ) );
#ifndef SING_NDEBUG
      if( __DEBUG__ )
        m_reds.Verify();
#endif
    }

    inline bool Reset()
    {
      m_active = false;

      m_itr = m_reds.m_hash.find(m_comp);

      if( m_itr == m_reds.m_hash.end() )
        return false;

      assume( m_itr->first == m_comp );

      m_current = (m_itr->second).begin();
      m_finish = (m_itr->second).end();

      if (m_current == m_finish)
        return false;

//      m_active = true;
      return true;
    }

    const CLeadingTerm& Current() const
    {
      assume( m_active );
      assume( m_current != m_finish );

      return *(*m_current);
    }

    inline bool MoveNext()
    {
      assume( m_current != m_finish );

      if( m_active )
        ++m_current;
      else
        m_active = true; // for Current()

      // looking for the next good entry
      for( ; m_current != m_finish; ++m_current )
      {
        assume( Current().CheckLT( m_reds.m_L ) );

        if( Current().DivisibilityCheck(m_product, m_not_sev, m_rBaseRing) )
        {
#ifndef SING_NDEBUG
          if( __DEBUG__ )
          {
            Print("CDivisorEnumerator::MoveNext::est LS: q is divisible by LS[%d] !:((, diviser is: ", 1 + Current().label());
            dPrint(Current().lt(), m_rBaseRing, m_rBaseRing, 0);
          }
#endif
//          m_active = true;
          assume( Current().CheckLT( m_reds.m_L ) );
          return true;
        }
        assume( Current().CheckLT( m_reds.m_L ) );        
      }

      // the end... :(
      assume( m_current == m_finish );

      m_active = false;
      return false;
    }
};



bool CReducerFinder::IsDivisible(const poly product) const
{
#ifndef SING_NDEBUG
  if( __DEBUG__ )
    Verify();
#endif
  
  assume( product != NULL );

  // NOTE: q may have no coeff!!!  

  CDivisorEnumerator itr(*this, product);
  if( !itr.Reset() )
    return false;

  return itr.MoveNext();

/*
  const ring& r = m_rBaseRing;

  const long comp = p_GetComp(product, r);
  const unsigned long not_sev = ~p_GetShortExpVector(product, r);

  assume( comp >= 0 );

  CReducersHash::const_iterator it = m_hash.find(comp); // same module component

  assume( m_L != NULL );

  if( it == m_hash.end() )
    return false;
  // assume comp!

  const TReducers& reducers = it->second;

  for(TReducers::const_iterator vit = reducers.begin(); vit != reducers.end(); vit++ )
  {
    assume( (*vit)->CheckLT( m_L ) );

    if( (*vit)->DivisibilityCheck(product, not_sev, r) )
    {
      if( __DEBUG__ )
      {
        Print("_FindReducer::Test LS: q is divisible by LS[%d] !:((, diviser is: ", 1 + (*vit)->label());
        dPrint((*vit)->lt(), r, r, 0);
      }

      return true;
    }
  }

  return false;
*/
}


#ifndef SING_NDEBUG
void CReducerFinder::Verify() const
{
  const ring& r = m_rBaseRing;

  for( CReducersHash::const_iterator it = m_hash.begin(); it != m_hash.end(); it++)
  {
    const TReducers& reducers = it->second;

    for(TReducers::const_iterator vit = reducers.begin(); vit != reducers.end(); vit++ )
    {
      assume( (*vit)->CheckLT( m_L ) );
      
      const poly p = (*vit)->lt();
      
      const unsigned long p_sev = (*vit)->sev();
      assume( p_sev == p_GetShortExpVector(p, r) );
      
      assume( p_GetComp(p, r) == it->first );

      const int k = (*vit)->label();        
      assume( m_L->m[k] == p );
       
      pp_Test(p, r, r);
    }
  }
}
   
   

void CReducerFinder::DebugPrint() const
{
  const ring& r = m_rBaseRing;

  for( CReducersHash::const_iterator it = m_hash.begin(); it != m_hash.end(); it++)
  {
    Print("Hash Key: %ld, Values: \n", it->first);
    const TReducers& reducers = it->second;

    for(TReducers::const_iterator vit = reducers.begin(); vit != reducers.end(); vit++ )
    {
      assume( (*vit)->CheckLT( m_L ) );
      
      const int k = (*vit)->label();
      const poly p = (*vit)->lt();
       
      pp_Test(p, r, r);
      
      assume( m_L->m[k] == p );

      const unsigned long p_sev = (*vit)->sev();
      assume( p_sev == p_GetShortExpVector(p, r) );
      
      assume( p_GetComp(p, r) == it->first );

      Print("L[%d]: ", k); dPrint(p, r, r, 0); Print("SEV: %ld\n", p_sev);
      
      assume( m_L->m[k] == p );
    }
  }
}
#endif

/// TODO:
class CDivisorEnumerator2: public SchreyerSyzygyComputationFlags
{
  private:
    const CReducerFinder& m_reds;
    const poly m_multiplier, m_term;
    const unsigned long m_not_sev;
    const long m_comp;

    CReducerFinder::CReducersHash::const_iterator m_itr;
    CReducerFinder::TReducers::const_iterator m_current, m_finish;

    bool m_active;

  public:
    CDivisorEnumerator2(const CReducerFinder& self, const poly m, const poly t):
        SchreyerSyzygyComputationFlags(self),
        m_reds(self),
        m_multiplier(m), m_term(t),
        m_not_sev(~p_GetShortExpVector(m, t, m_rBaseRing)),
        m_comp(p_GetComp(t, m_rBaseRing)),
        m_itr(), m_current(), m_finish(),
        m_active(false)
    {
      assume( m_comp >= 0 );
      assume( m_reds.m_L != NULL );
      assume( m_multiplier != NULL );
      assume( m_term != NULL );
      
      assume( m != NULL );
      assume( t != NULL );  
      assume ( !n_IsZero( p_GetCoeff(m, m_rBaseRing), m_rBaseRing ) );
      assume ( !n_IsZero( p_GetCoeff(t, m_rBaseRing), m_rBaseRing ) );
      
//      assume( p_GetComp(m_multiplier, m_rBaseRing) == 0 );
#ifndef SING_NDEBUG
      if( __DEBUG__ )
        m_reds.Verify();
#endif
    }

    inline bool Reset()
    {
      m_active = false;

      m_itr = m_reds.m_hash.find(m_comp);

      if( m_itr == m_reds.m_hash.end() )
        return false;

      assume( m_itr->first == m_comp );

      m_current = (m_itr->second).begin();
      m_finish = (m_itr->second).end();

      if (m_current == m_finish)
        return false;

      return true;
    }

    const CLeadingTerm& Current() const
    {
      assume( m_active );
      assume( m_current != m_finish );

      return *(*m_current);
    }

    inline bool MoveNext()
    {
      assume( m_current != m_finish );

      if( m_active )
        ++m_current;
      else
        m_active = true;


      // looking for the next good entry
      for( ; m_current != m_finish; ++m_current )
      {
        assume( Current().CheckLT( m_reds.m_L ) );

        if( Current().DivisibilityCheck(m_multiplier, m_term, m_not_sev, m_rBaseRing) )
        {
#ifndef SING_NDEBUG
          if( __DEBUG__ )
          {
            Print("CDivisorEnumerator::MoveNext::est LS: q is divisible by LS[%d] !:((, diviser is: ", 1 + Current().label());
            dPrint(Current().lt(), m_rBaseRing, m_rBaseRing, 0);
          }
#endif
//          m_active = true;
          assume( Current().CheckLT( m_reds.m_L ) );        
          return true;

        }
        assume( Current().CheckLT( m_reds.m_L ) );        
      }

      // the end... :(
      assume( m_current == m_finish );

      m_active = false;
      return false;
    }
};

poly CReducerFinder::FindReducer(const poly multiplier, const poly t,
                                 const poly syzterm,
                                 const CReducerFinder& syz_checker) const
{
#ifndef SING_NDEBUG
  if( __DEBUG__ )
  {
    Verify();
    syz_checker.Verify();
  }
#endif

  CDivisorEnumerator2 itr(*this, multiplier, t);
  if( !itr.Reset() )
    return NULL;

  // don't care about the module component of multiplier (as it may be the syzygy term)
  // product = multiplier * t?
  const ring& r = m_rBaseRing;

  assume( multiplier != NULL ); assume( t != NULL );

  const ideal& L = m_L; assume( L != NULL ); // for debug/testing only!

  long c = 0;

  if (syzterm != NULL)
    c = p_GetComp(syzterm, r) - 1;

  assume( c >= 0 && c < IDELEMS(L) );

  if (__DEBUG__ && (syzterm != NULL))
  {
    const poly m = L->m[c]; assume( m != NULL ); assume( pNext(m) == NULL );

    //  def @@c = leadcomp(syzterm); int @@r = int(@@c);
    //  def @@product = leadmonomial(syzterm) * L[@@r];
    poly lm = p_Mult_mm( leadmonom(syzterm, r, true), m, r); // !NC :(
    poly pr = p_Mult_q( leadmonom(multiplier, r, true), leadmonom(t, r, false), r); // !NC :(

    assume( p_EqualPolys(lm, pr, r) );

    p_Delete(&lm, r);
    p_Delete(&pr, r);
  }

#ifndef SING_NDEBUG
  if( __DEBUG__ )
  {
    Verify();
    syz_checker.Verify();
  }
#endif
  
  const BOOLEAN to_check = (syz_checker.IsNonempty()); // __TAILREDSYZ__ &&

  const poly q = p_New(r); pNext(q) = NULL;

  if( __DEBUG__ )
    p_SetCoeff0(q, 0, r); // for printing q

  while( itr.MoveNext() )
  {
    assume( itr.Current().CheckLT( L ) ); // ???
    
    const poly p = itr.Current().lt(); // ???
    const int k  = itr.Current().label();

    p_ExpVectorSum(q, multiplier, t, r); // q == product == multiplier * t // TODO: do it once?
    p_ExpVectorDiff(q, q, p, r); // (LM(product) / LM(L[k]))

    p_SetComp(q, k + 1, r);
    p_Setm(q, r);

    // cannot allow something like: a*gen(i) - a*gen(i)
    if (syzterm != NULL && (k == c))
      if (p_ExpVectorEqual(syzterm, q, r))
      {
#ifndef SING_NDEBUG
        if( __DEBUG__ )
        {
          Print("_FindReducer::Test SYZTERM: q == syzterm !:((, syzterm is: ");
          dPrint(syzterm, r, r, 0);
        }
#endif
        assume( itr.Current().CheckLT( L ) ); // ???
        continue;
      }

    // while the complement (the fraction) is not reducible by leading syzygies
    if( to_check && syz_checker.IsDivisible(q) )
    {
#ifndef SING_NDEBUG
      if( __DEBUG__ )
      {
        PrintS("_FindReducer::Test LS: q is divisible by LS[?] !:((: ");
      }
#endif
      assume( itr.Current().CheckLT( L ) ); // ???
      continue;
    }

    number n = n_Mult( p_GetCoeff(multiplier, r), p_GetCoeff(t, r), r);
    p_SetCoeff0(q, n_InpNeg( n_Div(n, p_GetCoeff(p, r), r), r), r);
    n_Delete(&n, r);

#ifndef SING_NDEBUG
    if( __DEBUG__ )
    {
      Verify();
      syz_checker.Verify();
    }
#endif
    
    assume( itr.Current().CheckLT( L ) ); // ???
    return q;
  }

/*
  const long comp = p_GetComp(t, r); assume( comp >= 0 );
  const unsigned long not_sev = ~p_GetShortExpVector(multiplier, t, r); // !

//   for( int k = IDELEMS(L)-1; k>= 0; k-- )
//   {
//     const poly p = L->m[k];
//
//     if ( p_GetComp(p, r) != comp )
//       continue;
//
//     const unsigned long p_sev = p_GetShortExpVector(p, r); // to be stored in m_hash!!!

   // looking for an appropriate diviser p = L[k]...
  CReducersHash::const_iterator it = m_hash.find(comp); // same module component

  if( it == m_hash.end() )
    return NULL;

  // assume comp!

  assume( m_L != NULL );

  const TReducers& reducers = it->second;

  for(TReducers::const_iterator vit = reducers.begin(); vit != reducers.end(); vit++ )
  {

    const poly p = (*vit)->lt(); // ??
    const int k = (*vit)->label();

    assume( L->m[k] == p ); // CheckLT

//    const unsigned long p_sev = (*vit)->sev();
//    assume( p_sev == p_GetShortExpVector(p, r) );

//    if( !p_LmShortDivisibleByNoComp(p, p_sev, product, not_sev, r) )
//      continue;

    if( !(*vit)->DivisibilityCheck(multiplier, t, not_sev, r) )
      continue;


//    if (p_sev & not_sev) continue;
//    if( !_p_LmDivisibleByNoComp(p, multiplier, t, r) ) continue;


    p_ExpVectorSum(q, multiplier, t, r); // q == product == multiplier * t
    p_ExpVectorDiff(q, q, p, r); // (LM(product) / LM(L[k]))

    p_SetComp(q, k + 1, r);
    p_Setm(q, r);

    // cannot allow something like: a*gen(i) - a*gen(i)
    if (syzterm != NULL && (k == c))
      if (p_ExpVectorEqual(syzterm, q, r))
      {
        if( __DEBUG__ )
        {
          Print("_FindReducer::Test SYZTERM: q == syzterm !:((, syzterm is: ");
          dPrint(syzterm, r, r, 0);
        }

        continue;
      }

    // while the complement (the fraction) is not reducible by leading syzygies
    if( to_check && syz_checker.IsDivisible(q) )
    {
      if( __DEBUG__ )
      {
        PrintS("_FindReducer::Test LS: q is divisible by LS[?] !:((: ");
      }

      continue;
    }

    number n = n_Mult( p_GetCoeff(multiplier, r), p_GetCoeff(t, r), r);
    p_SetCoeff0(q, n_InpNeg( n_Div(n, p_GetCoeff(p, r), r), r), r);
    n_Delete(&n, r);

    return q;
  }
*/

  p_LmFree(q, r);
  
#ifndef SING_NDEBUG
  if( __DEBUG__ )
  {
    Verify();
    syz_checker.Verify();
  }
#endif

  return NULL;

}


poly CReducerFinder::FindReducer(const poly product, const poly syzterm, const CReducerFinder& syz_checker) const
{

#ifndef SING_NDEBUG
  if( __DEBUG__ )
  {
    Verify();
    syz_checker.Verify();
  }
#endif

  CDivisorEnumerator itr(*this, product);
  if( !itr.Reset() )
    return NULL;


  
  const ring& r = m_rBaseRing;

  assume( product != NULL );

  const ideal& L = m_L; assume( L != NULL ); // for debug/testing only!

  long c = 0;

  if (syzterm != NULL)
    c = p_GetComp(syzterm, r) - 1;

  assume( c >= 0 && c < IDELEMS(L) );

  if (__DEBUG__ && (syzterm != NULL))
  {
    const poly m = L->m[c];

    assume( m != NULL ); assume( pNext(m) == NULL );

    poly lm = p_Mult_mm(leadmonom(syzterm, r), m, r);
    assume( p_EqualPolys(lm, product, r) );

    //  def @@c = leadcomp(syzterm); int @@r = int(@@c);
    //  def @@product = leadmonomial(syzterm) * L[@@r];

    p_Delete(&lm, r);
  }

#ifndef SING_NDEBUG
  if( __DEBUG__ )
  {
    Verify();
    syz_checker.Verify();
  }
#endif

  const BOOLEAN to_check = (syz_checker.IsNonempty()); // __TAILREDSYZ__ &&

  const poly q = p_New(r); pNext(q) = NULL;

  if( __DEBUG__ )
    p_SetCoeff0(q, 0, r); // ONLY for printing q

  while( itr.MoveNext() )
  {
    assume( itr.Current().CheckLT( L ) ); // ???
    
    const poly p = itr.Current().lt(); // ??
    const int k  = itr.Current().label();

    p_ExpVectorDiff(q, product, p, r); // (LM(product) / LM(L[k]))
    p_SetComp(q, k + 1, r);
    p_Setm(q, r);

    // cannot allow something like: a*gen(i) - a*gen(i)
    if (syzterm != NULL && (k == c))
      if (p_ExpVectorEqual(syzterm, q, r))
      {
#ifndef SING_NDEBUG
        if( __DEBUG__ )
        {
          Print("_FindReducer::Test SYZTERM: q == syzterm !:((, syzterm is: ");
          dPrint(syzterm, r, r, 0);
        }
#endif
        assume( itr.Current().CheckLT( L ) ); // ???
        continue;
      }

    // while the complement (the fraction) is not reducible by leading syzygies
    if( to_check && syz_checker.IsDivisible(q) ) // ????? 
    {
#ifndef SING_NDEBUG
      if( __DEBUG__ )
      {
        PrintS("_FindReducer::Test LS: q is divisible by LS[?] !:((: ");
      }
#endif
      assume( itr.Current().CheckLT( L ) ); // ???
      continue;
    }

    p_SetCoeff0(q, n_InpNeg( n_Div( p_GetCoeff(product, r), p_GetCoeff(p, r), r), r), r);

    assume( itr.Current().CheckLT( L ) ); // ???
    
#ifndef SING_NDEBUG
    if( __DEBUG__ )
    {
      Verify();
      syz_checker.Verify();
    }
#endif
    
    return q;
  }



/*
  const long comp = p_GetComp(product, r);
  const unsigned long not_sev = ~p_GetShortExpVector(product, r);

  assume( comp >= 0 );

//   for( int k = IDELEMS(L)-1; k>= 0; k-- )
//   {
//     const poly p = L->m[k];
//
//     if ( p_GetComp(p, r) != comp )
//       continue;
//
//     const unsigned long p_sev = p_GetShortExpVector(p, r); // to be stored in m_hash!!!

   // looking for an appropriate diviser p = L[k]...
  CReducersHash::const_iterator it = m_hash.find(comp); // same module component

  if( it == m_hash.end() )
    return NULL;

  assume( m_L != NULL );

  const TReducers& reducers = it->second;

  const BOOLEAN to_check = (syz_checker.IsNonempty()); // __TAILREDSYZ__ &&

  const poly q = p_New(r); pNext(q) = NULL;

  if( __DEBUG__ )
    p_SetCoeff0(q, 0, r); // for printing q

  for(TReducers::const_iterator vit = reducers.begin(); vit != reducers.end(); vit++ )
  {
    const poly p = (*vit)->lt(); // ???

    assume( p_GetComp(p, r) == comp );

    const int k = (*vit)->label();

    assume( L->m[k] == p ); // CheckLT

    const unsigned long p_sev = (*vit)->sev();

    assume( p_sev == p_GetShortExpVector(p, r) );

    if( !p_LmShortDivisibleByNoComp(p, p_sev, product, not_sev, r) )
      continue;

//     // ... which divides the product, looking for the _1st_ appropriate one!
//     if( !p_LmDivisibleByNoComp(p, product, r) ) // included inside  p_LmShortDivisibleBy!
//       continue;

    p_ExpVectorDiff(q, product, p, r); // (LM(product) / LM(L[k]))
    p_SetComp(q, k + 1, r);
    p_Setm(q, r);

    // cannot allow something like: a*gen(i) - a*gen(i)
    if (syzterm != NULL && (k == c))
      if (p_ExpVectorEqual(syzterm, q, r))
      {
        if( __DEBUG__ )
        {
          Print("_FindReducer::Test SYZTERM: q == syzterm !:((, syzterm is: ");
          dPrint(syzterm, r, r, 0);
        }

        continue;
      }

    // while the complement (the fraction) is not reducible by leading syzygies
    if( to_check && syz_checker.IsDivisible(q) )
    {
      if( __DEBUG__ )
      {
        PrintS("_FindReducer::Test LS: q is divisible by LS[?] !:((: ");
      }

      continue;
    }

    p_SetCoeff0(q, n_InpNeg( n_Div( p_GetCoeff(product, r), p_GetCoeff(p, r), r), r), r);
    return q;
  }
*/

  p_LmFree(q, r);

#ifndef SING_NDEBUG
  if( __DEBUG__ )
  {
    Verify();
    syz_checker.Verify();
  }
#endif
  
  return NULL;
}



CLCM::CLCM(const ideal& L, const SchreyerSyzygyComputationFlags& flags):
    SchreyerSyzygyComputationFlags(flags), std::vector<bool>(),
    m_compute(false), m_N(rVar(flags.m_rBaseRing))
{
  const ring& R = m_rBaseRing;
  assume( flags.m_rBaseRing == R );
  assume( R != NULL );

  assume( L != NULL );

  if( __TAILREDSYZ__ && !__HYBRIDNF__ && (L != NULL))
  {
    const int l = IDELEMS(L);

    assume( l > 0 );

    resize(l, false);

    for( int k = l - 1; k >= 0; k-- )
    {
      const poly a = L->m[k]; assume( a != NULL );

      for (unsigned int j = m_N; j > 0; j--)
        if ( !(*this)[j] )
          (*this)[j] = (p_GetExp(a, j, R) > 0);
    }

    m_compute = true;
  }
}


bool CLCM::Check(const poly m) const
{
  assume( m != NULL );
  if( m_compute && (m != NULL))
  {
    const ring& R = m_rBaseRing;

    assume( __TAILREDSYZ__ && !__HYBRIDNF__ );

    for (unsigned int j = m_N; j > 0; j--)
      if ( (*this)[j] )
        if(p_GetExp(m, j, R) > 0)
          return true;

    return false;

  } else return true;
}

CCacheCompare::CCacheCompare(): m_ring(currRing) {}


template class std::vector<bool>;
template class std::vector<CLeadingTerm const*>;
template class std::map< CReducerFinder::TComponentKey, CReducerFinder::TReducers >;

template class std::map<TCacheKey, TCacheValue, struct CCacheCompare>;
template class std::map<int, TP2PCache>;



END_NAMESPACE               END_NAMESPACE_SINGULARXX


// Vi-modeline: vim: filetype=c:syntax:shiftwidth=2:tabstop=8:textwidth=0:expandtab
