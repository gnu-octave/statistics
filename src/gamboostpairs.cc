/*
Copyright (C) 2026 Andreas Bertsatos <abertsatos@biol.uoa.gr>

This file is part of the statistics package for GNU Octave.

This program is free software; you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation; either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details.

You should have received a copy of the GNU General Public License along with
this program; if not, see <http://www.gnu.org/licenses/>.
*/

#include <cmath>
#include <vector>
#include <string>
#include <octave/oct.h>
#include <octave/ov-struct.h>
#include <octave/Cell.h>
#include "gamboost.cpp"

DEFUN_DLD(gamboostpairs, args, ,
          "-*- texinfo -*-\n\
@deftypefn {statistics} {@var{S} =} gamboostpairs (@var{X}, @var{R})\n\
\n\
Score every pair of predictors for an interaction.\n\
\n\
@code{@var{S} = gamboostpairs (@var{X}, @var{R})} lays the residuals @var{R}\n\
of a fitted additive model on the coarse grid of each pair of columns of\n\
@var{X} and returns the two-way analysis of variance @math{F} ratio testing\n\
what only the cells explain.  It is used to rank candidate interactions for\n\
@code{ClassificationGAM} and @code{RegressionGAM}, and it is not meant to be\n\
called directly.\n\
\n\
The @math{p}-values are deliberately not computed here.  Turning @math{F}\n\
into a probability needs @code{fcdf}, which the package already ships and\n\
which is verified against MATLAB, so the caller applies it rather than a\n\
second implementation being carried in the compiled engine.  Which pairs to\n\
keep is policy and belongs beside that.\n\
\n\
@itemize\n\
@item @var{X} is an @math{NxP} numeric matrix of predictors, and @var{P} must\n\
be at least 2 for any pair to exist.\n\
\n\
@item @var{R} is the @math{Nx1} residual vector of the additive fit.\n\
@end itemize\n\
\n\
@var{S} is a structure with the following fields, one row per pair, ordered\n\
as @code{nchoosek} orders them.\n\
\n\
@itemize\n\
@item @qcode{Pairs}, the @math{Mx2} matrix of predictor index pairs.\n\
@item @qcode{F}, the @math{Mx1} vector of @math{F} ratios.  A pair with too\n\
few observations, no spare degrees of freedom or no within-cell scatter\n\
scores @qcode{0}.\n\
@item @qcode{DF1} and @qcode{DF2}, the @math{Mx1} numerator and denominator\n\
degrees of freedom.\n\
@item @qcode{BinEdges}, a @math{1xP} cell of the coarse cut points each\n\
predictor was laid on.  The grid is fixed at eight equal-frequency bins,\n\
which is what MATLAB reports for pair detection at every sample size.\n\
@end itemize\n\
\n\
@seealso{gamboosttrain, gamboostpredict, fcdf, ClassificationGAM}\n\
@end deftypefn")
{
  if (args.length () != 2)
  {
    print_usage ();
  }

  if (! args(0).isnumeric () || args(0).iscomplex () || args(0).isempty ())
  {
    error ("gamboostpairs: X must be a numeric matrix.");
  }
  if (! args(1).isnumeric () || args(1).iscomplex () || args(1).isempty ()
      || args(1).columns () != 1)
  {
    error ("gamboostpairs: R must be a numeric column vector.");
  }
  if (args(0).rows () != args(1).rows ())
  {
    error ("gamboostpairs: X and R must have the same number of rows.");
  }

  Matrix X = args(0).matrix_value ();
  ColumnVector R = args(1).column_vector_value ();
  octave_idx_type n = X.rows ();
  octave_idx_type d = X.columns ();

  if (d < 2)
  {
    error ("gamboostpairs: X must have at least two columns.");
  }

  std::vector<BinnedPredictor> B ((std::size_t) d);
  Cell edges (1, d);
  for (octave_idx_type j = 0; j < d; j++)
  {
    ColumnVector xj (n);
    for (octave_idx_type i = 0; i < n; i++)
    {
      xj(i) = X(i, j);
    }
    B[(std::size_t) j] = gamb_bin (xj, GAMB_PAIR_EDGES);
    edges(j) = octave_value (B[(std::size_t) j].edges);
  }

  octave_idx_type m = d * (d - 1) / 2;
  Matrix pairs (m, 2);
  ColumnVector F (m), DF1 (m), DF2 (m);

  octave_idx_type row = 0;
  for (octave_idx_type j = 0; j < d - 1; j++)
  {
    for (octave_idx_type k = j + 1; k < d; k++)
    {
      GamPairStat S = gamb_pair_stat (B[(std::size_t) j], B[(std::size_t) k],
                                      R, j, k);
      pairs(row, 0) = (double) (j + 1);
      pairs(row, 1) = (double) (k + 1);
      F(row) = S.F;
      DF1(row) = S.df1;
      DF2(row) = S.df2;
      row++;
    }
  }

  octave_scalar_map out;
  out.assign ("Pairs", octave_value (pairs));
  out.assign ("F", octave_value (F));
  out.assign ("DF1", octave_value (DF1));
  out.assign ("DF2", octave_value (DF2));
  out.assign ("BinEdges", octave_value (edges));

  return ovl (out);
}

/*
%!test
%! ## Every pair is scored once, in nchoosek order, and the grid is reported.
%! x = randn (200, 4);
%! r = randn (200, 1);
%! S = gamboostpairs (x, r);
%! assert_equal (size (S.Pairs), [6, 2]);
%! assert_equal (S.Pairs, nchoosek (1:4, 2));
%! assert_equal (numel (S.F), 6);
%! assert_equal (numel (S.BinEdges), 4);

%!test
%! ## The detection grid is fixed at eight equal-frequency bins, so seven cut
%! ## points, whatever the sample size.  MATLAB reports seven at 60, 250, 1000
%! ## and 4000 observations alike.
%! for n = [60, 250, 1000]
%!   S = gamboostpairs (randn (n, 2), randn (n, 1));
%!   assert_equal (numel (S.BinEdges{1}), 7);
%! endfor

%!test
%! ## The cut points are the octiles, which is what MATLAB's coincide with.
%! x = [(1:800)', (1:800)'];
%! S = gamboostpairs (x, randn (800, 1));
%! q = linspace (0, 1, 9);
%! assert_equal (S.BinEdges{1}, quantile (x(:,1), q(2:end-1)), 1);

%!test
%! ## A planted interaction outscores every pair that carries none.  The
%! ## response depends on x1 * x2 alone, so that pair must rank first.
%! rand ('seed', 7);
%! randn ('seed', 7);
%! x = randn (400, 4);
%! r = x(:,1) .* x(:,2) + 0.1 * randn (400, 1);
%! S = gamboostpairs (x, r);
%! [~, best] = max (S.F);
%! assert_equal (S.Pairs(best,:), [1, 2]);

%!test
%! ## Structureless residuals score low: an F ratio near one is what a pair
%! ## with no interaction should give, and the planted pair above is orders
%! ## above it.
%! randn ('seed', 11);
%! S = gamboostpairs (randn (400, 3), randn (400, 1));
%! assert_equal (max (S.F) < 5, true);

%!test
%! ## A pair with nothing to explain scores zero rather than a NaN: a constant
%! ## residual leaves no within-cell scatter to divide by.
%! S = gamboostpairs (randn (50, 2), ones (50, 1));
%! assert_equal (S.F, 0);
%! assert_equal (S.DF1, 0);

%!test
%! ## Degrees of freedom follow the occupied grid, not its nominal size.
%! x = [(1:200)', (1:200)'];
%! S = gamboostpairs (x, randn (200, 1));
%! assert_equal (S.DF2 > 0, true);
%! assert_equal (S.DF1 > 0, true);

%!error<Invalid call> gamboostpairs (1)
%!error<gamboostpairs: X must be a numeric matrix.> ...
%! gamboostpairs ('a', [1;2])
%!error<gamboostpairs: R must be a numeric column vector.> ...
%! gamboostpairs ([1, 2; 3, 4], [1, 2])
%!error<gamboostpairs: X and R must have the same number of rows.> ...
%! gamboostpairs ([1, 2; 3, 4], [1; 2; 3])
%!error<gamboostpairs: X must have at least two columns.> ...
%! gamboostpairs ([1; 2; 3], [1; 2; 3])
*/
