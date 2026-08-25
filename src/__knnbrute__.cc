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

#include <cmath>
#include <string>
#include <vector>

// Exhaustive K nearest neighbours, without forming the distance matrix.
//
// The straightforward way to search exhaustively is to compute every distance,
// keep the whole M-by-N matrix, and then take the K smallest of each row.  The
// matrix is the problem: at four thousand points against four thousand it is
// 128 MB that is written once, read once and thrown away, and it puts a
// ceiling on N long before patience runs out.
//
// Here each distance is used the moment it is computed, against a running list
// of the K best for that query, and then discarded.  Memory is the size of the
// answer.  The list is kept sorted, so the output needs no sorting afterwards
// and the rejection test is a single compare against its last element, which
// fails for almost every point once the list has filled.
//
// Cosine is deliberately absent.  Its distance is 1 minus a quantity that is
// exactly 1 for parallel rows, so on collinear or one-dimensional data every
// true distance is 0 and the ordering is decided entirely by which way that
// subtraction rounds.  pdist2 does not agree with its own written formula
// there, returning 1.11e-16 where 1 - dot*sx*sy gives 0, so there is no fixed
// result to reproduce.  Cosine keeps the pdist2 route, where whatever it
// answers is what it has always answered.
//
// Ties resolve as a stable sort resolves them, the lower index first, because
// a candidate is admitted only when it is strictly better than the worst held
// and is inserted after every entry it ties with.

enum MetricCode
{
  MET_EUCLIDEAN = 0,
  MET_CITYBLOCK,
  MET_CHEBYCHEV,
  MET_MINKOWSKI
};

// Ordering of a candidate against a held entry, matching sort: a NaN comes
// after every number, and equal values keep the lower index.
template <typename T>
static inline bool
better (T da, octave_idx_type ia, T db, octave_idx_type ib)
{
  bool na = std::isnan (da);
  bool nb = std::isnan (db);
  if (na != nb)
  {
    return nb;
  }
  if (! na && da != db)
  {
    return da < db;
  }
  return ia < ib;
}

template <typename MT, typename T>
static octave_value_list
knnbrute (const MT& X, const MT& Y, octave_idx_type K, int metric, double q)
{
  const octave_idx_type n = X.rows ();
  const octave_idx_type p = X.columns ();
  const octave_idx_type m = Y.rows ();

  // Pack both sets so that a point's coordinates are contiguous.  Octave holds
  // them column-major, where the inner loop would stride by the row count and
  // miss the cache on every coordinate.
  std::vector<T> xr ((std::size_t) n * p);
  std::vector<T> yr ((std::size_t) m * p);
  for (octave_idx_type c = 0; c < p; c++)
  {
    for (octave_idx_type i = 0; i < n; i++)
    {
      xr[(std::size_t) i * p + c] = X(i, c);
    }
    for (octave_idx_type j = 0; j < m; j++)
    {
      yr[(std::size_t) j * p + c] = Y(j, c);
    }
  }

  Matrix OI (m, K);
  MT OD (m, K);

  std::vector<T> bd (K);
  std::vector<octave_idx_type> bi (K);

  for (octave_idx_type j = 0; j < m; j++)
  {
    const T *yp = &yr[(std::size_t) j * p];
    octave_idx_type held = 0;

    for (octave_idx_type i = 0; i < n; i++)
    {
      const T *xp = &xr[(std::size_t) i * p];

      // Euclidean and Minkowski are ranked on the sum and rooted at the end,
      // which saves a root per pair and cannot reorder anything, both roots
      // being monotone.
      T d;
      switch (metric)
      {
        case MET_CITYBLOCK:
        {
          d = 0;
          for (octave_idx_type c = 0; c < p; c++)
          {
            d += std::abs (xp[c] - yp[c]);
          }
          break;
        }
        case MET_CHEBYCHEV:
        {
          d = 0;
          for (octave_idx_type c = 0; c < p; c++)
          {
            T t = std::abs (xp[c] - yp[c]);
            if (t > d || std::isnan (t))
            {
              d = t;
            }
          }
          break;
        }
        case MET_MINKOWSKI:
        {
          d = 0;
          for (octave_idx_type c = 0; c < p; c++)
          {
            d += std::pow (std::abs (xp[c] - yp[c]), (T) q);
          }
          break;
        }
        default:
        {
          d = 0;
          for (octave_idx_type c = 0; c < p; c++)
          {
            T t = xp[c] - yp[c];
            d += t * t;
          }
          break;
        }
      }

      // The list is sorted, so its last entry is the one to beat.
      if (held == K && ! better<T> (d, i, bd[K-1], bi[K-1]))
      {
        continue;
      }

      octave_idx_type at = (held < K) ? held : K - 1;
      while (at > 0 && better<T> (d, i, bd[at-1], bi[at-1]))
      {
        bd[at] = bd[at-1];
        bi[at] = bi[at-1];
        at--;
      }
      bd[at] = d;
      bi[at] = i;
      if (held < K)
      {
        held++;
      }
    }

    for (octave_idx_type c = 0; c < K; c++)
    {
      OI(j, c) = (double) (bi[c] + 1);
      T v = bd[c];
      if (metric == MET_EUCLIDEAN)
      {
        v = std::sqrt (v);
      }
      else if (metric == MET_MINKOWSKI)
      {
        v = std::pow (v, (T) (1.0 / q));
      }
      OD(j, c) = v;
    }
  }

  return ovl (OI, OD);
}

DEFUN_DLD(__knnbrute__, args, ,
"-*- texinfo -*-\n\
@deftypefn {statistics} {[@var{idx}, @var{D}] =} __knnbrute__ (@var{X}, @var{Y}, @var{K}, @var{metric}, @var{param})\n\
\n\
Exhaustive nearest-neighbour search without forming the distance matrix.\n\
Internal; called by the nearest-neighbour searches and not meant to be used\n\
directly.\n\
\n\
@var{X} holds the reference points and @var{Y} the queries, one per row.\n\
@var{metric} is one of @qcode{'euclidean'}, @qcode{'cityblock'},\n\
@qcode{'chebychev'} or @qcode{'minkowski'}, and @var{param}\n\
the exponent for @qcode{'minkowski'} and empty otherwise.\n\
\n\
@var{idx} and @var{D} are the @math{M*K} indices and distances of the @var{K}\n\
nearest reference points to each query, in increasing distance, ordered as\n\
sorting the full row would order them.\n\
\n\
@end deftypefn")
{
  if (args.length () != 5)
  {
    print_usage ();
  }

  if (! args(0).isnumeric () || args(0).iscomplex () || args(0).isempty ())
  {
    error ("__knnbrute__: X must be a real numeric matrix.");
  }
  if (! args(1).isnumeric () || args(1).iscomplex () || args(1).isempty ())
  {
    error ("__knnbrute__: Y must be a real numeric matrix.");
  }
  if (args(0).columns () != args(1).columns ())
  {
    error ("__knnbrute__: X and Y must have the same number of columns.");
  }
  if (! args(2).is_scalar_type () || ! args(2).isnumeric ())
  {
    error ("__knnbrute__: K must be a real scalar.");
  }
  if (! args(3).is_string ())
  {
    error ("__knnbrute__: METRIC must be a character vector.");
  }

  const octave_idx_type n = args(0).rows ();
  double kd = args(2).scalar_value ();
  if (kd != std::floor (kd) || kd < 1.0 || kd > (double) n)
  {
    error ("__knnbrute__: K must be an integer between 1 and rows (X).");
  }
  octave_idx_type K = (octave_idx_type) kd;

  std::string name = args(3).string_value ();
  int metric;
  if (name == "euclidean")
  {
    metric = MET_EUCLIDEAN;
  }
  else if (name == "cityblock")
  {
    metric = MET_CITYBLOCK;
  }
  else if (name == "chebychev")
  {
    metric = MET_CHEBYCHEV;
  }
  else if (name == "minkowski")
  {
    metric = MET_MINKOWSKI;
  }
  else
  {
    error ("__knnbrute__: unsupported METRIC '%s'.", name.c_str ());
  }

  double q = 2.0;
  if (metric == MET_MINKOWSKI)
  {
    if (! args(4).is_scalar_type () || ! args(4).isnumeric ())
    {
      error ("__knnbrute__: PARAM must be the minkowski exponent.");
    }
    q = args(4).scalar_value ();
    if (! (q > 0) || ! octave::math::isfinite (q))
    {
      error ("__knnbrute__: PARAM must be a positive finite scalar.");
    }
  }

  if (args(0).is_single_type () || args(1).is_single_type ())
  {
    return knnbrute<FloatMatrix, float> (args(0).float_matrix_value (),
                                         args(1).float_matrix_value (),
                                         K, metric, q);
  }

  return knnbrute<Matrix, double> (args(0).matrix_value (),
                                   args(1).matrix_value (), K, metric, q);
}

/*
%!test
%! ## It answers what pdist2 followed by a sort answers, exactly.
%! X = [1, 2; 3, 4; 5, 6; 7, 8; 2, 2; -1, 0];
%! Y = [2, 3; 6, 7; 0, 0];
%! for m = {"euclidean", "cityblock", "chebychev"}
%!   D = pdist2 (X, Y, m{1})';
%!   [sv, so] = sort (D, 2);
%!   for k = 1:6
%!     [idx, dst] = __knnbrute__ (X, Y, k, m{1}, []);
%!     assert_equal (idx, so(:,1:k));
%!     assert_equal (dst, sv(:,1:k));
%!   endfor
%! endfor

%!test
%! ## Minkowski takes its exponent and agrees at each of them.
%! X = [1, 2; 3, 4; 5, 6; 7, 8];
%! Y = [2, 3; 6, 7];
%! for q = [1, 1.5, 2, 3]
%!   D = pdist2 (X, Y, "minkowski", q)';
%!   [sv, so] = sort (D, 2);
%!   [idx, dst] = __knnbrute__ (X, Y, 3, "minkowski", q);
%!   assert_equal (idx, so(:,1:3));
%!   assert_equal (dst, sv(:,1:3), 1e-12);
%! endfor

%!test
%! ## Ties keep the lower index, as a stable sort does.
%! X = [0, 0; 0, 0; 0, 0];
%! [idx, dst] = __knnbrute__ (X, [1, 1], 2, "euclidean", []);
%! assert_equal (idx, [1, 2]);
%! assert_equal (dst, [sqrt(2), sqrt(2)], 1e-12);

%!test
%! ## A NaN in the data puts that neighbour last, never first.
%! X = [0, 0; NaN, 0; 1, 1];
%! [idx, dst] = __knnbrute__ (X, [0, 0], 3, "euclidean", []);
%! assert_equal (idx(1:2), [1, 3]);
%! assert_equal (idx(3), 2);
%! assert_equal (isnan (dst(3)), true);

%!test
%! ## Single data gives single distances; the indices stay double.
%! [idx, dst] = __knnbrute__ (single ([1, 1; 4, 4]), single ([0, 0]), 1, ...
%!                            "euclidean", []);
%! assert_equal (class (dst), 'single');
%! assert_equal (class (idx), 'double');

%!test
%! ## K may be the whole reference set.
%! X = [3, 0; 1, 0; 2, 0];
%! [idx, dst] = __knnbrute__ (X, [0, 0], 3, "euclidean", []);
%! assert_equal (idx, [2, 3, 1]);
%! assert_equal (dst, [1, 2, 3], 1e-12);

%!error<Invalid call to __knnbrute__> __knnbrute__ (ones (2, 2), ones (1, 2), 1)
%!error<__knnbrute__: X must be a real numeric matrix.> ...
%! __knnbrute__ ({1}, ones (1, 2), 1, "euclidean", [])
%!error<__knnbrute__: X and Y must have the same number of columns.> ...
%! __knnbrute__ (ones (2, 3), ones (1, 2), 1, "euclidean", [])
%!error<__knnbrute__: K must be an integer between 1 and rows \(X\).> ...
%! __knnbrute__ (ones (2, 2), ones (1, 2), 3, "euclidean", [])
%!error<__knnbrute__: unsupported METRIC 'cosine'.> ...
%! __knnbrute__ (ones (2, 2), ones (1, 2), 1, "cosine", [])
%!error<__knnbrute__: PARAM must be a positive finite scalar.> ...
%! __knnbrute__ (ones (2, 2), ones (1, 2), 1, "minkowski", -1)
*/
