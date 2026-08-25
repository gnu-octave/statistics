/*
Copyright (C) 2026 Andreas Bertsatos <abertsatos@biol.uoa.gr>

This file is part of the statistics package for GNU Octave.

This program is free software; you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation; either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
details.

You should have received a copy of the GNU General Public License along with
this program; if not, see <http://www.gnu.org/licenses/>.
*/

#include <octave/oct.h>

#include <algorithm>
#include <cmath>
#include <utility>
#include <vector>

// The K smallest entries of each row, by partial selection.
//
// A nearest-neighbour search spends most of its time here rather than in the
// distances: taking the five smallest of four thousand by sorting all four
// thousand costs O(N log N) per row where O(N) will do, and on a 4000-by-4000
// matrix the sort is roughly nine times the cost of forming the distances.
//
// The order returned is Octave's own, so this is a drop-in replacement for
// sort rather than an approximation of it: non-NaN before NaN, then by value,
// then by column index, which is what a stable sort gives for ties.

template <typename T>
struct ValIdxLess
{
  typedef std::pair<T, octave_idx_type> ValIdx;

  bool operator() (const ValIdx& a, const ValIdx& b) const
  {
    bool na = std::isnan (a.first);
    bool nb = std::isnan (b.first);
    if (na != nb)
    {
      // A NaN sorts after every number, as it does in sort.
      return nb;
    }
    if (! na && a.first != b.first)
    {
      return a.first < b.first;
    }
    // Ties, NaN against NaN included, keep the lower column first.
    return a.second < b.second;
  }
};

// The class of the distances is carried through, single in and single out, as
// the searches promise.  The indices are always double.
template <typename MT, typename T>
static octave_value_list
knnselect (const MT& D, octave_idx_type K)
{
  typedef std::pair<T, octave_idx_type> ValIdx;

  const octave_idx_type m = D.rows ();
  const octave_idx_type n = D.columns ();

  Matrix OI (m, K);
  MT OD (m, K);

  const T *dp = D.data ();
  std::vector<ValIdx> row (n);
  ValIdxLess<T> cmp;

  for (octave_idx_type i = 0; i < m; i++)
  {
    for (octave_idx_type j = 0; j < n; j++)
    {
      row[j] = ValIdx (dp[i + j * m], j);
    }

    // nth_element leaves the K smallest below the pivot in no order, so the
    // kept head is sorted afterwards; both together stay linear in N.
    if (K < n)
    {
      std::nth_element (row.begin (), row.begin () + K, row.end (), cmp);
    }
    std::sort (row.begin (), row.begin () + K, cmp);

    for (octave_idx_type j = 0; j < K; j++)
    {
      OI(i, j) = (double) (row[j].second + 1);
      OD(i, j) = row[j].first;
    }
  }

  return ovl (OI, OD);
}

DEFUN_DLD(__knnselect__, args, ,
"-*- texinfo -*-\n\
@deftypefn {statistics} {[@var{idx}, @var{D}] =} __knnselect__ (@var{dist}, @var{K})\n\
\n\
The @var{K} smallest entries of each row of @var{dist}.  Internal; called by\n\
the nearest-neighbour searches and not meant to be used directly.\n\
\n\
@var{dist} is an @math{M*N} matrix of distances whose rows are query points,\n\
and @var{K} a positive integer no greater than @math{N}.  @var{idx} is the\n\
@math{M*K} matrix of column indices of the @var{K} smallest entries of each\n\
row, and @var{D} their values, in the class @var{dist} carries.\n\
\n\
The result is identical to sorting each row and taking the first @var{K}\n\
columns, ties and @qcode{NaN} included, and is obtained by partial selection\n\
rather than by a full sort.\n\
\n\
@end deftypefn")
{
  if (args.length () != 2)
  {
    print_usage ();
  }

  if (! args(0).isnumeric () || args(0).iscomplex () || args(0).isempty ())
  {
    error ("__knnselect__: DIST must be a real numeric matrix.");
  }
  if (! args(1).is_scalar_type () || ! args(1).isnumeric ()
      || args(1).iscomplex ())
  {
    error ("__knnselect__: K must be a real scalar.");
  }

  const octave_idx_type n = args(0).columns ();
  double kd = args(1).scalar_value ();

  if (kd != std::floor (kd) || kd < 1.0 || kd > (double) n)
  {
    error ("__knnselect__: K must be an integer between 1 and columns (DIST).");
  }
  octave_idx_type K = (octave_idx_type) kd;

  if (args(0).is_single_type ())
  {
    return knnselect<FloatMatrix, float> (args(0).float_matrix_value (), K);
  }

  return knnselect<Matrix, double> (args(0).matrix_value (), K);
}

/*
%!test
%! ## The result is sort's, not an approximation of it.
%! D = [3, 1, 2, 1; 5, 5, 5, 5; 9, 8, 1, 2; 0, -1, -1, 4];
%! [sv, so] = sort (D, 2);
%! for k = 1:4
%!   [idx, dst] = __knnselect__ (D, k);
%!   assert_equal (idx, so(:,1:k));
%!   assert_equal (dst, sv(:,1:k));
%! endfor

%!test
%! ## Ties keep the lower column first, as a stable sort does.
%! [idx, dst] = __knnselect__ ([5, 5, 5, 5], 3);
%! assert_equal (idx, [1, 2, 3]);
%! assert_equal (dst, [5, 5, 5]);

%!test
%! ## A NaN sorts after every number and never displaces one.
%! [idx, dst] = __knnselect__ ([3, NaN, 1, NaN, 2], 3);
%! assert_equal (idx, [3, 5, 1]);
%! assert_equal (dst, [1, 2, 3]);

%!test
%! ## A row of nothing but NaN keeps them in column order.
%! [idx, dst] = __knnselect__ ([NaN, NaN, NaN], 2);
%! assert_equal (idx, [1, 2]);
%! assert_equal (isnan (dst), [true, true]);

%!test
%! ## Inf is a value like any other and comes before NaN.
%! [idx, dst] = __knnselect__ ([Inf, NaN, 2, -Inf], 3);
%! assert_equal (idx, [4, 3, 1]);
%! assert_equal (dst, [-Inf, 2, Inf]);

%!test
%! ## K equal to the width returns the whole row, sorted.
%! [idx, dst] = __knnselect__ ([4, 2, 9, 1], 4);
%! assert_equal (idx, [4, 2, 1, 3]);
%! assert_equal (dst, [1, 2, 4, 9]);

%!test
%! ## Single distances come back single; the indices stay double.
%! [idx, dst] = __knnselect__ (single ([4, 2, 9, 1]), 2);
%! assert_equal (class (dst), 'single');
%! assert_equal (class (idx), 'double');
%! assert_equal (dst, single ([1, 2]));

%!test
%! ## It agrees with sort over random matrices at every width of K.
%! rand ("seed", 42);
%! for t = 1:50
%!   A = round (rand (6, 9) * 4);
%!   k = 1 + mod (t, 9);
%!   [idx, dst] = __knnselect__ (A, k);
%!   [sv, so] = sort (A, 2);
%!   assert_equal (idx, so(:,1:k));
%!   assert_equal (dst, sv(:,1:k));
%! endfor

%!error<Invalid call to __knnselect__> __knnselect__ (ones (2, 2))
%!error<__knnselect__: DIST must be a real numeric matrix.> ...
%! __knnselect__ ({1, 2}, 1)
%!error<__knnselect__: DIST must be a real numeric matrix.> ...
%! __knnselect__ (ones (2, 2) * i, 1)
%!error<__knnselect__: K must be a real scalar.> ...
%! __knnselect__ (ones (2, 2), [1, 2])
%!error<__knnselect__: K must be an integer between 1 and columns \(DIST\).> ...
%! __knnselect__ (ones (2, 3), 0)
%!error<__knnselect__: K must be an integer between 1 and columns \(DIST\).> ...
%! __knnselect__ (ones (2, 3), 4)
%!error<__knnselect__: K must be an integer between 1 and columns \(DIST\).> ...
%! __knnselect__ (ones (2, 3), 1.5)
*/
