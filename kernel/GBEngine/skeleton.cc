/**
  \brief implementation of classes for double description method
  \author John Perry
  \version 1.0
  \date October 2014
  \copyright The University of Southern Mississippi
*/

#ifndef SKELETON_C
#define SKELETON_C

#include <iostream>
#include <cstdlib>
#include "skeleton.h"

constraint::constraint(ulong num_variables, const long coeffs [])
{
  nvars = num_variables;
  coefficients = new long[nvars];
  for (ulong i = 0; i < nvars; ++i)
    coefficients[i] = coeffs[i];
}

constraint::constraint(const vector<long> &coeffs)
{
  nvars = coeffs.size();
  coefficients = new long[nvars];
  for (ulong i = 0; i < nvars; ++i)
    coefficients[i] = coeffs[i];
}

constraint::constraint(const constraint &old_constraint)
{
  nvars = old_constraint.nvars;
  coefficients = new long[nvars];
  for (ulong i = 0; i < nvars; ++i)
    coefficients[i] = old_constraint.coefficients[i];
}

// ordering is lexicographic
bool operator < (const constraint &first, const constraint &second)
{
  bool result = !(first == second);
  bool checking = result;
  for (ulong i = 0; checking and i < first.nvars; ++i)
    if (first[i] != second[i])
    {
      checking = false;
      result = first[i] < second[i];
    }
  return result;
}

bool operator == (const constraint &first, const constraint &second)
{
  bool result = true;
  for (ulong i = 0; result and i < first.get_number_of_variables(); ++i)
    if (first[i] != second[i])
      result = false;
  return result;
}

// formatted as sum of products of coefficients and variables
ostream & operator<<(ostream & ostr, const constraint &c)
{
  bool first = true;
  ostr << "0 ≤ ";
  for (ulong i = 0; i < c.nvars; ++i)
  {
    if (c[i] < 0)
    {
      if (first) { ostr << '-'; first = false; } else ostr << "- ";
      if (c[i] != -1)
        ostr << -c[i] << 'x' << i << ' ';
      else
        ostr << 'x' << i << ' ';
    }
    else if (c[i] > 0)
    {
      if (first) { first = false; } else ostr << "+ ";
      if (c[i] != 1)
        ostr << c[i] << 'x' << i << ' ';
      else
        ostr << 'x' << i << ' ';
    }
    else
    {
      // do nothing when c[i] == 0
    }
  }
  return ostr;
}

constraint::~constraint()
{
  delete [] coefficients;
}

ray::ray(ulong dimension, long direction)
{
  dim = dimension;
  coords = new ulonglong[dim];
  for (ulong i = 0; i < dim; ++i) coords[i] = 0;
  coords[direction] = 1;
  // we always know the active constraints that define this ray
  // so we populate known_active_constraints accordingly
  long *coeffs = new long[dim];
  for (ulong i = 0; i < (ulong )dim; ++i)
    if (i != direction)
    {
      coeffs[i] = 1;
      add_active_constraint(constraint(dim, coeffs));
      coeffs[i] = 0;
    }
  delete [] coeffs;
}

ray::ray(ulong dimension, const ulonglong entries [])
{
  dim = dimension;
  coords = new ulonglong[dim];
  for (ulong i = 0; i < dim; ++i)
    coords[i] = entries[i];
}

ray::ray(const vector<ulonglong> &entries)
{
  dim = entries.size();
  coords = new ulonglong[dim];
  for (ulong i = 0; i < dim; ++i)
    coords[i] = entries[i];
}

ray::ray(const ray &old_ray)
{
  dim = old_ray.dim;
  coords = new ulonglong[dim];
  for (ulong i = 0; i < dim; ++i)
    coords[i] = old_ray.coords[i];
  // need to copy each constraint, to avoid errors in the destructor
  set<constraint> old_constraints = old_ray.known_active_constraints;
  for (set<constraint>::iterator citer = old_constraints.begin();
        citer != old_constraints.end();
        ++citer)
    known_active_constraints.insert(*citer);
}

ray::~ray()
{
  delete [] coords;
}

void ray::simplify_ray()
{
  ulonglong * w = coords;
  ulonglong gcd = 0;
  for (int i = 0; i < get_dimension(); ++i)
  {
    if (gcd == 0)
      gcd = w[i];
    else if (w[i] != 0)
    {
      ulonglong r = (gcd < w[i] ? gcd : w[i]);
      ulonglong s = (gcd < w[i] ? w[i] : gcd);
      while (r != 0)
      {
        ulonglong t = s % r;
        s = r;
        r = t;
      }
      gcd = s;
    }
  }
  if (gcd > 1)
    for (int i = 0; i < get_dimension(); ++i)
      w[i] /= gcd;
}

long long ray::obtain_dot_product(const constraint &hyperplane) const
{
  long long result = 0;
  for (ulong i = 0; i < dim; ++i)
    result += coords[i] * hyperplane[i];
  return result;
}

bool ray::add_active_constraint(const constraint & hyperplane, bool verify)
{
  bool can_add = true;
  if (verify)
    can_add = is_active_at(hyperplane);
  if (can_add)
    known_active_constraints.insert(hyperplane);
  return can_add;
}

void ray::find_and_add_active_constraints(const vector<constraint> &old_constraints)
{
  for (vector<const constraint>::iterator citer = old_constraints.begin();
       citer != old_constraints.end();
       ++citer
      )
    if (this->is_active_at(*citer))
      add_active_constraint(*citer);
}

ray operator*(ulonglong a, const ray &r)
{
  ulong d = r.get_dimension();
  ulonglong *coords = new ulonglong[d];
  for (ulong i = 0; i < d; ++i)
    coords[i] = a*r[i];
  ray result(d, coords);
  delete [] coords;
  return result;
}

ulonglong operator*(const ray &r1, const ray &r2)
{
  ulonglong result = 0;
  for (ulong i = 0; i < r1.get_dimension(); ++i)
    result += r1[i]*r2[i];
  return result;
}

ulonglong operator*(const ray &r1, vector<long> &r2)
{
  ulonglong result = 0;
  for (ulong i = 0; i < r1.get_dimension(); ++i)
    result += r1[i]*r2[i];
  return result;
}

ulonglong operator*( vector<long> &r1, const ray &r2)
{
  ulonglong result = 0;
  for (ulong i = 0; i < r1.size(); ++i)
    result += r1[i]*r2[i];
  return result;
}

ray operator+(const ray &r1, const ray &r2)
{
  ulong d = r1.get_dimension();
  ulonglong *coords = new ulonglong[d];
  for (ulonglong i = 0; i < d; ++i)
    coords[i] = r1[i] + r2[i];
  ray result(d, coords);
  delete [] coords;
  return result;
}

ray operator-(const ray &r1, const ray &r2)
{
  ulong d = r1.get_dimension();
  ulonglong *coords = new ulonglong[d];
  for (ulonglong i = 0; i < d; ++i)
    coords[i] = r1[i] - r2[i];
  ray result(d, coords);
  delete [] coords;
  return result;
}

ray ray_sum(const set<ray> &rs)
{
  ulong d = 0;
  ulonglong *coords = NULL;
  for (set<ray>::iterator riter = rs.begin(); riter != rs.end(); ++riter)
  {
    ray r = *riter;
    if (coords == NULL)
    {
      d = r.get_dimension();
      coords = new ulonglong[d];
      for (ulong i = 0; i < d; ++i) coords[i] = 0;
    }
    for (ulong i = 0; i < d; ++i)
    {
      coords[i] += r[i];
    }
  }
  ray result(d, coords);
  delete [] coords;
  return result;
}

bool operator==(const ray &r1, const ray &r2)
{
  bool result = true;
  for (int ulong i = 0; result and i < r1.get_dimension(); ++i)
    result = r1[i] == r2[i];
  return result;
}

ostream & operator<<(ostream & ostr, const ray &r)
{
  ostr << "( ";
  ulong i;
  for (i = 0; i < r.dim - 1; ++i)
    ostr << r[i] << ", ";
  ostr << r[i] << " )";
  return ostr;
}

ray & ray::operator=(const ray &other)
{
  if (!(*this == other))
  {
    // resize coords if need be
    if (dim != other.dim)
    {
      dim = other.dim;
      delete [] coords;
      coords = new ulonglong[dim];
    }
    // copy coords
    for (ulong i = 0; i < dim; ++i)
      coords[i] = other.coords[i];
    set<constraint> old_constraints = other.known_active_constraints;
    // need to copy each constraint, to avoid errors in the destructor
    for (set<constraint>::iterator citer = old_constraints.begin();
          citer != old_constraints.end();
          ++citer)
      known_active_constraints.insert(*citer);
  }
  return *this;
}

void ray::swap(ray &other)
{
  ulonglong tmpval;
  for (ulong i = 0; i < dim; ++i)
  {
    tmpval = coords[i];
    coords[i] = other.coords[i];
    other.coords[i] = tmpval;
  }
}

bool operator<(const ray &first_ray, const ray &second_ray)
{
  bool result = true;
  ulong i = 0;
  bool equal = true;
  while (equal and i < first_ray.dim)
  {
    if (first_ray[i] != second_ray[i])
    {
      equal = false;
      result = first_ray[i] < second_ray[i];
    }
    ++i;
  }
  return result and (not equal);
}

edge::edge(const ray &first_ray, const ray &second_ray)
    : first(first_ray), second(second_ray)
{
  if (first < second)
    first.swap(second);
}

edge::edge(const edge &old_edge)
    : first(old_edge.first), second(old_edge.second)
{
  // nothing to do
}

ostream & operator<<(ostream & ostr, const edge &e)
{
  ostr << "{ " << e.first << " , " << e.second << " }";
  return ostr;
}

bool operator < (const edge &first, const edge &second)
{
  bool result = true;
  if (first.first < second.first)
  {
    // do nothing
  }
  else if (second.first < first.first)
    result = false;
  else // first entries equal; look at second
    if (first.second < second.second)
    {
      // do nothing
    }
    else
      result = false;
  return result;
}

edge & edge::operator=(const edge &other)
{
  if (!(*this == other))
  {
    first = other.first;
    second = other.second;
  }
  return *this;
}

void skeleton::common_initialization(ulong dimension)
{
  //cout << "creating skeleton with dimension " << dimension << endl;
  dim = dimension;
  // initialize the constraints
  long * constr_coords = new long[dim];
  for (int i = 0; i < dim; ++i) constr_coords[i] = 0;
  // add the constraint xi >= 0 for each i = 0, ..., dim - 1
  for (int i = 0; i < dim; ++i)
  {
    constr_coords[i] = 1;
    constraints.push_back(constraint(dim, constr_coords));
    constr_coords[i] = 0;
  }
  delete [] constr_coords;
  // initialize the rays, one for each axis
  for (int i = 0; i < dim; ++i)
  {
    ray new_ray(dim, i);
    // add the constraints active at this ray: xj >= 0 for j =/= i
    for (int j = 0; j < dim; ++j)
      if (i != j)
        new_ray.add_active_constraint(constraints[j]);
    rays.insert(new_ray);
  }
  // initialize the edges
  // currently, all the rays are adjacent
  for (set<ray>::iterator riter = rays.begin(); riter != rays.end(); ++riter)
    for (set<ray>::iterator siter = rays.begin(); siter != rays.end(); ++siter)
      if (*riter != *siter)
      {
        edge new_edge(*riter, *siter);
        edges.insert(new_edge);
      }
}

skeleton::skeleton(ulong dimension)
{
  common_initialization(dimension);
}

skeleton::skeleton(ulong dimension, const vector<constraint> &constraints)
        //: skeleton(dimension)
{
  common_initialization(dimension);
  ddm(constraints);
}

skeleton::skeleton(const skeleton &old_skeleton)
        : constraints(old_skeleton.constraints), rays(old_skeleton.rays),
          edges(old_skeleton.edges), dim(old_skeleton.dim)
          
{
  // nothing more to do
  // WARNING: this may not be correct;
  // I may need to copy each constraint, ray, and edge individually
}

skeleton::~skeleton()
{
}

bool skeleton::ddm(const constraint &constraint)
{
  cout << "processing constraint " << constraint << endl;
  // innocent until proven guilty
  bool consistent = true;
  // sort the rays into the ones above, below, or on the constraint
  set<ray> rays_above, rays_below, rays_on;
  for (set<ray>::iterator riter = rays.begin(); riter != rays.end(); ++riter)
  {
    long long dp = (*riter) * constraint; // overloaded * as dot product :-)
    if (dp > 0)
    {
      rays_above.insert(*riter);
      cout << *riter << " is above constraint\n";
    }
    else if (dp < 0)
    {
      rays_below.insert(*riter);
      cout << *riter << " is below constraint\n";
    }
    else
    {
      ray old_ray = *riter;
      old_ray.add_active_constraint(constraint);
      rays_on.insert(old_ray);
      cout << *riter << " is on constraint\n";
    }
  }
  cout << rays_above.size() << " rays above; " << rays_below.size() << " rays below; " << rays_on.size() << " rays on\n";
  // check for constitency
  if (rays_above.size() == 0)
  {
    consistent = false;
    //cout << "inconsistent\n";
  }
  // proceed only if constraint is consistent, and *not* redundant;
  // redundancy can be checked by making sure
  // that at least one ray is below the constraint
  if (consistent and rays_below.size() != 0)
  {
    set<edge> edges_above, edges_on;
    for (set<edge>::iterator eiter = edges.begin(); eiter != edges.end(); ++eiter)
    {
      edge e = *eiter;
      ray u = e.get_first_ray();
      ray v = e.get_second_ray();
      //cout << "edge " << u << ',' << v << endl;
      // identify the edges that lie above and on this constraint
      if (
          (rays_above.find(u) != rays_above.end() or rays_on.find(u) != rays_on.end())
          and
          (rays_above.find(v) != rays_above.end() or rays_on.find(v) != rays_on.end())
         )
      {
        cout << "old edge preserved: " << u << ',' << v << "\n";
        edges_above.insert(e);
      }
    }
    for (set<edge>::iterator eiter = edges.begin(); eiter != edges.end(); ++eiter)
    {
      edge e = *eiter;
      ray u = e.get_first_ray();
      ray v = e.get_second_ray();
      // identify edges that pass through the constraint
      // (one ray above, one ray below)
      if (rays_above.find(u) != rays_above.end() and rays_below.find(v) != rays_below.end())
      {
        ray w = (u * constraint)*v - (v * constraint)*u; // overloaded * :-)
        w.simplify_ray();
        w.add_active_constraint(constraint);
        w.find_and_add_active_constraints(constraints);
        rays_on.insert(w);
        edges_on.insert(edge(u,w));
        cout << "new ray (u,v) is " << w << " with constraints " << endl;
      }
      else if (rays_above.find(v) != rays_above.end() and rays_below.find(u) != rays_below.end())
      {
        ray w = (v * constraint)*u - (u * constraint)*v;
        w.simplify_ray();
        w.add_active_constraint(constraint);
        w.find_and_add_active_constraints(constraints);
        rays_on.insert(w);
        edges_on.insert(edge(v,w));
        cout << "new ray (v,u) is " << w << " with constraints " << endl;
      }
    }
    // clear the old rays, add the new ones (above and on the constraint)
    rays.clear();
    for (set<ray>::iterator riter = rays_above.begin(); riter != rays_above.end(); ++riter)
    {
      rays.insert(*riter);
    }
    cout << "inserted rays above; rays on is\n";
    for (set<ray>::iterator riter = rays_on.begin(); riter != rays_on.end(); ++riter) { cout << '\t' << *riter << endl; }
    for (set<ray>::iterator riter = rays_on.begin(); riter != rays_on.end(); ++riter)
    {
      cout << "inserting " << *riter << endl;
      cout << "return value: " << *(get<0>(rays.insert(*riter)));
      for (set<ray>::iterator siter = rays.begin(); siter != rays.end(); ++siter) { cout << '\t' << *siter << endl; }
    }
    cout << rays.size() << " rays\n";
    for (set<ray>::iterator riter = rays.begin(); riter != rays.end(); ++riter) { cout << '\t' << *riter << endl; }
    // add the good constraint
    constraints.push_back(constraint);
    // determine new edges
    set<edge> edges_new = adjacencies_by_graphs(rays_on);
    // combine new edges with old ones that are known to be valid
    edges = union_of_edge_sets(union_of_edge_sets(edges_above, edges_on), edges_new);
    cout << edges.size() << " edges\n";
    for (set<edge>::iterator eiter = edges.begin(); eiter != edges.end(); ++eiter) { cout << *eiter << ' '; } cout << '\n';
  }
  return consistent;
}

bool skeleton::ddm(const vector<constraint> &new_constraints)
{
  // innocent until proven guilty
  bool consistent = true;
  //cout << "adding " << new_constraints.size() << "constraints\n";
  // process each constraint sequentially
  for (
        vector<const constraint>::iterator nciter = new_constraints.begin();
        consistent and nciter != new_constraints.end();
        ++nciter
      )
  {
    //cout << "adding constraint " << *nciter << endl;
    consistent = ddm(*nciter);
    //cout << "\tconsistent? " << consistent << endl;
  }
  return consistent;
}

int number_of_common_constraints(
    const set<constraint> &a, const set<constraint> &b
)
{
  int result = 0;
  for (set<constraint>::iterator aiter = a.begin(); aiter != a.end(); ++aiter)
  {
    //cout << "checking " << *aiter << " in other: " << (b.find(*aiter) != b.end()) << endl;
    if (b.find(*aiter) != b.end())
      ++result;
  }
  return result;
}

set<constraint> intersections_of_active_constraints(
    const set<constraint> &a, const set<constraint> &b
)
{
  // highly unoptimized, but off the top of my head i don't know how to do better
  set<constraint> result;
  for (set<constraint>::iterator aiter = a.begin(); aiter != a.end(); ++aiter)
    if (b.find(*aiter) != b.end())
      result.insert(*aiter);
  return result;
}

bool is_first_subset_of_second(
    const set<constraint>a, const set<constraint>b
)
{
  // highly unoptimized, but off the top of my head i don't know how to do better
  bool result = true;
  set<constraint>::iterator aiter = a.begin();
  while (result and aiter != a.end())
  {
    result = b.find(*aiter) != b.end();
    ++aiter;
  }
  return result;
}

set<edge> union_of_edge_sets(const set<edge> a, const set<edge> b)
{
  // highly unoptimized, but off the top of my head i don't know how to do better
  set<edge> result;
  for (set<edge>::iterator eiter = a.begin(); eiter != a.end(); ++eiter)
    result.insert(*eiter);
  for (set<edge>::iterator eiter = b.begin(); eiter != b.end(); ++eiter)
    result.insert(*eiter);
  return result;
}

set<edge> skeleton::adjacencies_by_graphs(set<ray> new_rays)
{
  set<edge> new_edges;
  set<ray> tested_rays;
  // loop through each new ray, examining active constraints shared with other rays
  for (set<ray>::iterator riter = new_rays.begin(); riter != new_rays.end(); ++riter)
  {
    ray u = *riter;
    tested_rays.insert(u);
    set<constraint> zero_u = u.get_known_active_constraints();
    // D's rays have at least dim - 2 active constraints in common with u
    // (see Proposition 3 in Zolotych's paper)
    set<ray> D;
    for (set<ray>::iterator siter = new_rays.begin(); siter != new_rays.end(); ++siter)
      if (*riter != *siter)
      {
        ray v = *siter;
        set<constraint> zero_v = v.get_known_active_constraints();
        //cout << "checking constraints of " << u << " against " << v  << " for " << dim << endl;
        if (number_of_common_constraints(zero_u, zero_v) >= dim - 2)
        {
          //cout << "accept " << u << ',' << v << " from active constraints\n";
          D.insert(v);
        } else {
          //cout << "reject " << u << ',' << v << " from active constraints\n";
        }
      }
    // check u with each v in D, making sure their active constraints
    // are not a subset of the active constraints of any w in D
    // (see Proposition 4 (graph test) in Zolotych's paper
    for (set<ray>::iterator diter = D.begin(); diter != D.end(); ++diter)
    {
      ray v = *diter;
      if (tested_rays.find(v) == tested_rays.end()) // avoid doubling edges
      {
        set<constraint> zero_v = v.get_known_active_constraints();
        // WARNING: I have commented out the following line, because it seems
        // unnecessary: v is in D iff the size of the intersection is at least
        // dim - 2. If there are unexpected bugs, this commenting should be
        // reconsidered.
        // if (intersections_of_active_constraints(zero_u, zero_v).size() >= dim - 2)
        {
          bool can_be_added = true;
          set<constraint> Zuv = intersections_of_active_constraints(zero_u, zero_v);
          for (
               set<ray>::iterator dditer = D.begin();
               dditer != D.end() and can_be_added;
               ++dditer
              )
          {
            ray w = *dditer;
            if (!(w == v))
              if (is_first_subset_of_second(Zuv,w.get_known_active_constraints()))
              {
                //cout << "rejecting " << u << ',' << v << " because of " << w << endl;
                can_be_added = false;
              }
          }
          if (can_be_added)
          {
            edge new_edge(u, v);
            //cout << "edge " << new_edge << " passes all criteria\n";
            new_edges.insert(new_edge);
          }
        }
      }
    }
  }
  return new_edges;
}

ostream & operator << (ostream & ostr, const skeleton &skel)
{
  // header, start constraints
  ostr << "Skeleton defined by constraints" << endl;
  for (
       vector<const constraint>::iterator citer=skel.constraints.begin();
       citer != skel.constraints.end();
       ++citer
      )
    ostr << '\t' << *citer << endl;
  // rays
  ostr << "has " << skel.rays.size() << " rays" << endl;
  for (set<ray>::iterator riter=skel.rays.begin(); riter != skel.rays.end(); ++riter)
    ostr << '\t' << *riter << endl;
  //edges
  ostr << "connected in " << skel.edges.size() << " edges" << endl;
  for (set<edge>::iterator eiter=skel.edges.begin(); eiter != skel.edges.end(); ++eiter)
    ostr << '\t' << *eiter << endl;
  // footer
  ostr << "End of skeleton" << endl;
  return ostr;
}

skeleton & skeleton::operator=(const skeleton & other)
{
  rays.clear();
  edges.clear();
  constraints.clear();
  dim = other.dim;
  for (
       set<ray>::iterator siter = other.rays.begin();
       siter != other.rays.end();
       ++siter
      )
    rays.insert(*siter);
  for (
       set<edge>::iterator eiter = other.edges.begin();
       eiter != other.edges.end();
       ++eiter
      )
    edges.insert(*eiter);
  for (
        vector<const constraint>::iterator citer = other.constraints.begin();
        citer != other.constraints.end();
        ++citer
      )
  constraints.push_back(*citer);
  return *this;
}

#endif