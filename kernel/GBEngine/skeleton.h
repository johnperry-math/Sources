/**
  \brief implementation of classes for double description method
  \author John Perry
  \version 1.0
  \date October 2014
  \copyright The University of Southern Mississippi
*/

#ifndef SKELETON_H
#define SKELETON_H

#include <set>
#include <vector>
#include <iostream>
#include "../polys.h"
using namespace std;

#define ulong     unsigned long
#define ulonglong unsigned long long

/**
  \brief a constraint \f$ c_1 x_1 + \ldots + c_n x_n \geq 0 \f$
  \author John Perry
  \version 1.0
  \date October 2014
  \copyright The University of Southern Mississippi
  \details This class encapsulates a simple constraint for a skeleton; that is,
  an inequality of the form \f$ c_1 x_1 + \ldots + c_n x_n \geq 0 \f$.
  Constraints can be ordered lexicographically using the less-than operator,
  allowing for their inclusion in ordered collections, such as sets.
  \todo Add a hash mechanism.
  \todo Have `ddm` save state, and roll back if one of the constraints does not work.
*/
class constraint
{

public:

  // construction

  /**
    Initialize constraint to the given coefficients.
    The resulting constraint is \f$ c_1x_1 + \cdots + c_nx_n \geq 0, \f$
    where \f$ c_i \f$ is the coefficient of \f$ x_i \f$.
    \param num_variables length of coeffs
    \param coeffs copies this array of coefficients
    \pre the size of the array needs to be at least as long as the dimension!
  */
  constraint(ulong, long []);

  /**
    Initialize constraint to the given coefficients.
    The resulting constraint is \f$ c_1x_1 + \cdots + c_nx_n \geq 0, \f$
    where \f$ c_i \f$ is the coefficient of \f$ x_i \f$.
    \param coeffs copies thiis vector of coefficients
    \post \c nvars will have the value equal to `coeffs.size()`
  */
  constraint(vector<long> &);

  /**
    Copies the coefficients of the other constraint,
    including the allocation of new memory.
  */
  constraint(const constraint &);

  // destruction

  /**
    Deletes memory allocated by the constructor.
    Currently, that means it deletes an array created by the constructors.
  */
  ~constraint();

  //access

  /**
    Returns the number of variables in the constraint.
  */
  inline ulong get_number_of_variables() const { return nvars; }

  /**
    Returns the coefficient indicated. Numbering starts at 0.
  */
  inline long operator[](ulong index) const { return coefficients[index]; };

  /**
    Lexicographic comparison of constraints.
    \warning This is unsafe when number of variables is not the same.
      It does not check, since the assumption is that you know what you're doing.
  */
  friend bool operator<(const constraint &a, const constraint &b);

  /**
    \warning This is unsafe when number of variables is not the same.
      It does not check, since the assumption is that you know what you're doing.
  */
  friend bool operator==(const constraint &a, const constraint &b);
  friend bool operator!=(const constraint &a, const constraint &b);

  /**
    \warning This is unsafe when number of variables is not the same.
      It does not check, since the assumption is that you know what you're doing.
  */
  friend bool operator!=(constraint &a, constraint &b);

  // I/O

  /**
    Output is of the form \f$ c_1 x_1 + \ldots + c_n x_n \f$ , where \f$ c_i \f$
    is the coefficient of \f$ x_i \f$.
  */
  friend ostream & operator<<(ostream &, const constraint &);

private:

  /** number of variables/coefficients */
  ulong nvars;
  /** coefficients */
  long * coefficients;

};

/**
  \brief a ray defined by nonnegative coordinates \f$(a_1,\ldots,a_n)\f$
  \author John Perry
  \version 1.0
  \date October 2014
  \copyright The University of Southern Mississippi
  \details This class encapsulates a ray, one major part of the definition of a skeleton.
  Rays can be initialized to a particular set of coefficients, or to a particular
  axis (which is then translated into the corresponding coefficients).

  A special feature is that a rays can track the constraints known to be
  active at the ray, allowing for more efficient computation in the double
  description method. Adding known constraints can be done with or without
  checking whether the constraint actually is active, so this should be done
  with care.
*/
class ray
{

public:

  // construction

  /**
    Creates a ray with the given number of variables, all set to 0.
    The optional second argument specifies a direction, and sets that coordinate
    to 1. In this case, there is no need to set the ray's known active constraints,
    as this is known and populated automatically.
    \pre The dimension should be greater than zero. While the direction need not
      be specified&hellip; (see postcondition)
    \post &hellip;the result when the direction is zero is a zero ray.
      If the direction is \f$ i \f$, then the result is the \f$i\f$th canonical vector.
  */
  ray(ulong, long = -1);

  /**
    Creates a ray with the given number of variables,
    with coordinates set to the value of the array.
    \pre the size of the array needs to be at least as long
      as the number of variables!
  */
  ray(ulong, ulonglong []);

  /**
    Creates a ray whose coordinates are given by the vector.
    \post The dimension of this ray will equal the number of entries in the vector,
      and the values of their entries will be equal.
  */
  ray(vector<ulonglong> &);

  /**
    Copies the coordinates of the other ray.
    Allocates new memory, and copies the active constraints.
  */
  ray(const ray &);

  // destruction

  /**
    Deletes memory allocated by the constructor.
    Currently, that means it deletes `coords`.
  */
  ~ray();

  // access
  
  /** Returns the dimension of this ray. */
  inline ulong get_dimension() const { return dim; };

  /** Returns the entry indicated. Numbering starts at 0. */
  inline ulonglong operator[](ulong index) const { return coords[index]; };

  /**
    Indicates whether the two rays are equal.
    \warning This is unsafe when number of variables is not the same.
      It does not check, since the assumption is that you know what you're doing.
  */
  friend bool operator==(const ray &, const ray &);
  friend bool operator!=(const ray &, const ray &);

  /**
    Indicates whether the two rays are unequal.
    \warning This is unsafe when number of variables is not the same.
      It does not check, since the assumption is that you know what you're doing.
  */
  inline friend bool operator!=(ray &r, ray &s) { return !(r==s); }

  /** Synonym for []. I have no idea why I added this. */
  inline ulonglong coordinate(ulong index) { return coords[index]; };

  /**
    Returns `true` if and only if the hyperplane is active at this ray.
    Practically speaking, if the hyperplane is defined by the vector
    \f$ \mathbf c \f$ and the ray is defined by \f$ \mathbf r \f$ ,
    this function returns true if and only if \f$ c\cdot r = 0 \f$.
  */
  inline bool is_active_at(const constraint &hyperplane)
  {
    return 0 == obtain_dot_product(hyperplane);
  };

  /**
    Returns `true` if and only if this ray is above the hyperplane.
    Practically speaking, if the hyperplane is defined by the vector
    \f$ \mathbf c \f$ and the ray is defined by \f$ \mathbf r \f$ ,
    this function returns true if and only if \f$ c\cdot r > 0 \f$.
  */
  inline bool is_above(constraint &hyperplane)
  {
    return 0 < obtain_dot_product(hyperplane);
  };

  /**
    Returns `true` if and only if this ray is below the hyperplane.
    Practically speaking, if the hyperplane is defined by the vector
    \f$ \mathbf c \f$ and the ray is defined by \f$ \mathbf r \f$ ,
    this function returns true if and only if \f$ c\cdot r < 0 \f$.
  */
  inline bool is_below(constraint &hyperplane)
  {
    return 0 > obtain_dot_product(hyperplane);
  };

  /** Returns \c true if and only if the coordinates of the two rays are equal. */
  friend bool operator == (const ray &, const ray &);

  /**
    Lexicographic comparison of rays.
    \warning This is unsafe when dimension is not the same.
      It does not check, since the assumption is that you know what you're doing.
  */
  friend bool operator < (const ray &, const ray &);

  /**
    Convenience function to compute dot product between ray and constraint.
    \warning This is unsafe when dimension is not the same.
      It does not check, since the assumption is that you know what you're doing.
  */
  long long obtain_dot_product(const constraint &) const;

  // modify

  /**
    Simplifies the ray by dividing its components by the least common denominator.
  */
  void simplify_ray();

  /**
    Adds a constraint that is claimed to be active, and returns `true`.
    If the optional boolean argument is set to `true` (default is `false`),
    then it first checks whether the constraint really is active.
  */
  bool add_active_constraint(const constraint &, bool = false);

  /**
    Finds and adds the constraints in this vector that are active at this ray.
    \warning This is unsafe when dimension is not the same
      as the number of variables.
      It does not check, since the assumption is that you know what you're doing.
  */
  void find_and_add_active_constraints(vector<constraint> &);

  /** Returns the list of known active constraints. */
  inline set<constraint> * get_known_active_constraints()
  { return & known_active_constraints; };

  /**
    Assignment operator; assigns the value of `other` to `this`.
    \warning This is unsafe when dimension is not the same.
      It does not check, since the assumption is that you know what you're doing.
  */
  ray & operator=(const ray &);

  /**
    Swap two rays of equal dimension by swapping their data,
      avoiding memory reallocation.
    \warning This is unsafe when dimension is not the same.
      It does not check, since the assumption is that you know what you're doing.
  */
  void swap(ray &);

  // I/O

  /** Output is of the form \f$(r_1, \ldots, r_n)\f$. */
  friend ostream & operator<<(ostream &, const ray &);

private:

  ulong dim; /**< number of entries in \c coords */

  ulonglong * coords; /**< coordinates of the ray */

  /**
    list of constraints known to be active at this ray --
    or, if the client was lazy, believed to be active at this ray
  */
  set<constraint> known_active_constraints;

};

/**@{*/

/**
  \defgroup RayArithmetic arithmetic on rays
  \brief convenience functions for arithmetic on rays
  \warning These operations are unsafe when the dimensions of two rays differ,
    or when the dimension of a ray differs from the number of variables in a
    constraint.
    It does not check, since the assumption is that you know what you're doing.
*/

/**
  \ingroup RayArithmetic
  Multiply every coordinate in the given ray by the given scalar.
  \warning This is unsafe when dimension is not the same.
    It does not check, since the assumption is that you know what you're doing.
*/
ray operator*(ulonglong, ray &);

/**
  \ingroup RayArithmetic
  Add the two rays.
  \warning This is unsafe when dimension is not the same.
    It does not check, since the assumption is that you know what you're doing.
*/
ray operator+(ray &, ray &);

/**
  \ingroup RayArithmetic
  Subtract the two rays.
  \warning This is unsafe when dimension is not the same.
    It does not check, since the assumption is that you know what you're doing.
*/
ray operator-(const ray &, const ray &);

/**
  \ingroup RayArithmetic
  Add all the rays in a set.
  \warning This is unsafe when dimension is not the same.
    It does not check, since the assumption is that you know what you're doing.
*/
ray ray_sum(const set<ray> &);

/**
  \ingroup RayArithmetic
  Compute the dot product on the rays.
  \warning This is unsafe when dimension is not the same.
    It does not check, since the assumption is that you know what you're doing.
*/
ulonglong operator*(const ray &, const ray &);

/**
  \ingroup RayArithmetic
  Compute the dot product between the ray and the constraint.
  \warning This is unsafe when dimension is not the same.
    It does not check, since the assumption is that you know what you're doing.
*/
inline long long operator*(const ray &r, const constraint &c)
{ return r.obtain_dot_product(c); }

/**
  \ingroup RayArithmetic
  Compute the dot product between the ray and the constraint.
  \warning This is unsafe when dimension is not the same.
    It does not check, since the assumption is that you know what you're doing.
*/
inline long long operator*(constraint &c, ray &r)
{ return r.obtain_dot_product(c); }

/**@}*/

/**
  \brief an edge \f$(r_1,r_2)\f$ connecting the two rays \f$ r_1 \f$ and \f$ r_2 \f$
  \author John Perry
  \version 1.0
  \date October 2014
  \copyright The University of Southern Mississippi
  \details  This class encapsulates an edge, the other major part of a skeleton.
  Edges describe how the rays of the skeleton are connected.
  Edges are ordered, so that the smaller ray always comes first.

  \warning An edge's rays should have the same dimension.
    To start with, it doesn't make mathematical sense to &ldquo;join&rdquo;
    two rays of different dimension.
    Moreover, comparison of edges requires comparison of rays,which requires
    that the rays have the same dimension.
    (But you wouldn't be dumb enough to do this in the first place.)
*/
class edge

{

public:

  // construction

  /** Creates a new edge that joins the two rays. */
  edge(const ray &, const ray &);

  /** Copies the rays in \c other to two new rays. */
  edge(const edge &);

  // destruction: destroy nothing inside!

  ~edge() {} /**< Does nothing beyond what the compiler would do. */

  // access

  /** Returns the first ray listed in this edge. */
  inline ray get_first_ray() { return first; };

  /** Returns the second ray listed in this edge. */
  inline ray get_second_ray() { return second; };

  /**
    Compares two edges lexicographically.
    If the first ray in \c this edge is smaller, then \c this edge is smaller.
    Otherwise, if the first rays are equal, and the second ray in \c this edge
    is smaller, then \c this edge is smaller.
  */
  friend bool operator<(const edge &, const edge &);

  // I/O

  /**
    Output has the form \f$ \{ \mathbf{r}_1, \mathbf{r}_2 \} \f$
    where \f$ \mathbf{r}_1 \f$ is the first ray in this edge, etc.
  */
  friend ostream & operator<<(ostream &, const edge &);

  /** Assignment operator */
  edge & operator=(const edge &);

  /**
    Equal if and only if the first and second rays are true.
    (We can restrict ourselves to this because we have ordered the rays.)
  */
  inline friend bool operator==(const edge &e1, const edge &e2)
  { return e1.first == e2.first && e1.second == e2.second; }

private:

  ray first, second; /**< the rays defining this edge */

};

/**@{*/

/**
  \defgroup SetArithmetic arithmetic on sets
  \brief convenience functions for computing intersections, unions, and related
    properties
*/

/**
  \ingroup SetArithmetic
  Returns the number of constraints common to both sets.
*/
int number_of_common_constraints(
    set<constraint> &, set<constraint> &
);


/**
  \ingroup SetArithmetic
  Returns the intersection between the given sets of constraints.
*/
set<constraint> intersections_of_active_constraints(
      set<constraint> &, set<constraint> &
  );

/**
  \ingroup SetArithmetic
  Returns `true` if and only if the first set is a subset of the second.
*/
bool is_first_subset_of_second(
      set<constraint> &, set<constraint> &
  );

/**
  \ingroup SetArithmetic
  Returns the unions of two sets of edges.
*/
set<edge> union_of_edge_sets(const set<edge> &, const set<edge> &);

/**@}*/


/**
  \brief skeleton of a polyhedral cone, with methods allowing definition and refinement
  \author John Perry
  \version 1.0
  \date October 2014
  \copyright The University of Southern Mississippi
  \details  This class encapsulates the skeleton of a polyhedral cone,
  defined by a sequence of inequalities of the form
  \f$ c_1 x_1 + \cdots c_n x_n \geq 0 \f$.
  It also provides an implementation of the Double Description Method,
  an iterative algorithm for computing the skeleton of a cone.
  The iterative nature means that the cone can be updated with new constraints,
  passed to that algorithm, and the skeleton will be automatically recomputed.

  \warning Clients must ensure two things.
      -# First, the rays must
         have the same number \f$ m \f$ of dimensions, constraints must have
         the same number \f$ n \f$ of variables, and \f$ m=n \f$. Violating any
         of these three conditions will lead to unpredictable behavior.
         &mdash; well, no, that's too generous. Violating any of these three
         conditions will bring your program to a screeching halt!
      -# Second, if refining the cone, it is essential to check that
         the return value of `ddm` is `true`; for if it is not,
         then the cone is no longer be consistent.
         Please read the accompanying instructions.
*/
class skeleton
{

public:

  // construction

  /**
  */
  void common_initialization(ulong);

  /**
    Constructs a basic skeleton in the given number of dimensions,
    initialized to the axes, or (equivalently) to the set of constraints
    \f$ x_i \geq 0 \f$. The rays are informed of their active constraints.
    \pre the argument should be at least two
    \post the skeleton of the positive orthant
  */
  skeleton(ulong);

  /**
    Constructs a skeleton described by the given system of constraints.
    Practically speaking, it first generates a basic skeleton,
    then iterates on the given constraints.
    \pre `u.size() == v.size()` for all `u`, `v` in the vector
    \post unless the system supplied was inconsistent, a valid skeleton of the
        corresponding polyhedral cone
    \warning Your program will almost certainly fail if you do not respect the precondition.
  */
  skeleton(ulong, vector<constraint> &);

  /** Performs a deep copy of `other`. */
  skeleton(skeleton &);

  // destruction

  /**
    Currently does nothing the compiler wouldn't do.
  */
  ~skeleton();

  // access

  /** Returns the dimension of this skeleton. */
  inline ulong get_dimension() const { return dim; };

  /** Returns the number of rays defining the skeleton. */
  inline ulong get_number_of_rays() { return rays.size(); };

  /** Returns the rays that define the skeleton. */
  inline const set<ray> & get_rays() { return rays; };

  /** Returns the number of edges defining the skeleton. */
  inline ulong get_number_of_edges() { return edges.size(); };

  /** Returns the edges that define the skeleton. */
  inline set<edge> get_edges() { return edges; };

  /** Returns the number of constraints defining the skeleton. */
  inline ulong get_number_of_constraints() { return constraints.size(); };

  /** Returns the constraints that define the skeleton. */
  inline const vector<constraint> & get_constraints() { return constraints; };

  /** Returns the indicated constraint. Numbering starts at 0. */
  inline constraint get_constraint(int index) { return constraints[index]; };

  /** prints out the constraints, then the rays, then the edges. */
  friend ostream & operator<<(ostream &, const skeleton &);

  /** tests for consistency of a potentially new constraint. */
  inline bool is_consistent(const constraint & c) const
  {
    bool inconsistent = true;
    for (set<ray>::iterator riter = rays.begin(); inconsistent and riter != rays.end(); ++riter)
      if (((*riter) * c) > 0)
        inconsistent = false;
    return not inconsistent;
  }

  /** tests for consistency of a constraint generated by two monomials. */
  inline bool makes_consistent_constraint(poly t, poly u, bool show_data = false)
  {
    bool inconsistent = true;
    //if (show_data) pWrite(t);
    for (set<ray>::iterator riter = rays.begin(); inconsistent and riter != rays.end(); ++riter)
    {
      int d = 0;
      for (int i = 0; i < riter->get_dimension(); ++i)
        d += (*riter)[i] * (pGetExp(t,i+1) - pGetExp(u,i+1));
      if (d > 0)
        inconsistent = false;
      //if (show_data) { cout << d << ' '; if (!inconsistent) cout << *riter << endl; }
    }
    //if (show_data) cout << endl;
    return not inconsistent;
  }

  // modification

  /**
    Adds the indicated constraints (plural!) and re-computes the skeleton.

    \return `true` if and only if the new constraints are consistent with the
      current constraints

    \warning Checking the return value is crucial!
      If the function returns `false`, you have an inconsistent system!
      While the <i>present</i> cone will remain consistent,
      <b>the function will not roll back previous changes you have made</b>,
      so if you want to iterate again,
      your best bet is to copy the skeleton, and try that copy.
      Accept the new constraints only if that copy succeeds,
      in which case, you might as well discard the original, and keep the copy.
  */
  bool ddm(vector<constraint> &);

  /**
    Adds the indicated constraint (singular!) and re-computes the skeleton.

    \return `true` if and only if the new constraint is consistent with the
      current constraints

    \warning Checking the return value is crucial!
      If the function returns `false`, you have an inconsistent system!
      While the <i>present</i> cone will remain consistent,
      <b>the function will not roll back previous changes you have made</b>,
      so if you want to iterate again,
      your best bet is to copy the skeleton, and try that copy.
      Accept the new constraints only if that copy succeeds,
      in which case, you might as well discard the original, and keep the copy.
  */
  bool ddm(constraint &);

  /**
    Re-computes the edges in the skeleton using Zolotych's `GraphAdj` algorithm
    and returns the result.
  */
  set<edge> adjacencies_by_graphs(set<ray>);

  /**
    Assignment operator; empties current set & copies from other.
  */
  skeleton & operator=(const skeleton &);

private:

  int dim; /**< dimension of skeleton */

  set<ray> rays; /**< rays defining skeleton */

  set<edge> edges; /**< edges defining skeleton */

  vector<constraint> constraints; /**< constraints defining skeleton */

};

#endif
