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

DEFUN_DLD(gamboostinter, args, ,
          "-*- texinfo -*-\n\
@deftypefn {statistics} {@var{Mdl} =} gamboostinter (@var{X}, @var{Y}, @\n\
@var{F0}, @var{Method}, @var{Pairs}, @var{NumTrees}, @var{LearnRate}, @\n\
@var{MaxNumSplits})\n\
\n\
Boost trees over selected pairs of predictors.\n\
\n\
@code{@var{Mdl} = gamboostinter (@dots{})} fits the interaction phase of a\n\
generalized additive model, continuing from the additive prediction the\n\
predictor phase left rather than refitting it.  It is used by\n\
@code{ClassificationGAM} and @code{RegressionGAM}, and it is not meant to be\n\
called directly.\n\
\n\
@itemize\n\
@item @var{X} is an @math{NxP} numeric matrix of predictors and @var{Y} the\n\
@math{Nx1} response, as @code{gamboosttrain} takes them.\n\
\n\
@item @var{F0} is the @math{Nx1} additive prediction of the predictor phase.\n\
The interaction phase starts from it, so its deviance is where this phase\n\
begins.\n\
\n\
@item @var{Method} selects what is boosted, @qcode{1} the logistic deviance\n\
and @qcode{2} the squared error.\n\
\n\
@item @var{Pairs} is an @math{Mx2} matrix of predictor index pairs, one-based\n\
and within range.  Choosing them is the caller's business; see\n\
@code{gamboostpairs}.\n\
\n\
@item @var{NumTrees}, @var{LearnRate} and @var{MaxNumSplits} are the\n\
interaction phase's own budget, initial step and split limit.\n\
@end itemize\n\
\n\
@var{Mdl} is a structure with the following fields.\n\
\n\
@itemize\n\
@item @qcode{PairBinEdges}, a @math{1xP} cell of the coarse cut points the\n\
surfaces are held on.  Interactions are binned coarser than main effects: a\n\
tree limited to @var{MaxNumSplits} splits carves no more regions than that,\n\
so a finer grid buys nothing and costs memory in every pair.\n\
@item @qcode{PairValues}, a @math{1xM} cell of matrices, one value per cell of\n\
the pair's grid.\n\
@item @qcode{Intercept}, the constant the recentred surfaces gave up.  Add it\n\
to the intercept of the predictor phase.\n\
@item @qcode{NumTrees}, @qcode{ReasonForTermination}, @qcode{Deviance} and\n\
@qcode{Residuals}, as @code{gamboosttrain} reports them.\n\
@end itemize\n\
\n\
@seealso{gamboosttrain, gamboostpairs, gamboostpredict}\n\
@end deftypefn")
{
  if (args.length () != 8)
  {
    print_usage ();
  }

  if (! args(0).isnumeric () || args(0).iscomplex () || args(0).isempty ())
  {
    error ("gamboostinter: X must be a numeric matrix.");
  }
  if (! args(1).isnumeric () || args(1).iscomplex () || args(1).isempty ()
      || args(1).columns () != 1)
  {
    error ("gamboostinter: Y must be a numeric column vector.");
  }
  if (! args(2).isnumeric () || args(2).iscomplex () || args(2).isempty ()
      || args(2).columns () != 1)
  {
    error ("gamboostinter: F0 must be a numeric column vector.");
  }
  if (args(0).rows () != args(1).rows ()
      || args(0).rows () != args(2).rows ())
  {
    error ("gamboostinter: X, Y and F0 must have the same number of rows.");
  }

  Matrix X = args(0).matrix_value ();
  ColumnVector Y = args(1).column_vector_value ();
  ColumnVector F0 = args(2).column_vector_value ();
  octave_idx_type d = X.columns ();

  if (! args(3).is_scalar_type () || ! args(3).isnumeric ()
      || args(3).iscomplex ())
  {
    error ("gamboostinter: Method must be a numeric scalar.");
  }
  int method = args(3).int_value ();
  if (method != 1 && method != 2)
  {
    error ("gamboostinter: Method must be 1 or 2.");
  }

  if (! args(4).isnumeric () || args(4).iscomplex () || args(4).isempty ()
      || args(4).columns () != 2)
  {
    error ("gamboostinter: Pairs must be a numeric matrix with two columns.");
  }
  Matrix pairs = args(4).matrix_value ();

  for (octave_idx_type q = 0; q < pairs.rows (); q++)
  {
    for (octave_idx_type c = 0; c < 2; c++)
    {
      double v = pairs(q, c);
      if (v != std::floor (v) || v < 1.0 || v > (double) d)
      {
        error ("gamboostinter: Pairs must hold predictor indices "
               "between 1 and %d.", (int) d);
      }
    }
    if (pairs(q, 0) == pairs(q, 1))
    {
      error ("gamboostinter: a pair must name two different predictors.");
    }
  }

  if (! args(5).is_scalar_type () || ! args(5).isnumeric ()
      || args(5).iscomplex ())
  {
    error ("gamboostinter: NumTrees must be a numeric scalar.");
  }
  octave_idx_type maxtrees = (octave_idx_type) args(5).int_value ();
  if (maxtrees < 1)
  {
    error ("gamboostinter: NumTrees must be a positive integer.");
  }

  if (! args(6).is_scalar_type () || ! args(6).isnumeric ()
      || args(6).iscomplex ())
  {
    error ("gamboostinter: LearnRate must be a numeric scalar.");
  }
  double lrate = args(6).scalar_value ();
  if (! (lrate > 0.0) || lrate > 1.0)
  {
    error ("gamboostinter: LearnRate must be greater than 0 and at most 1.");
  }

  if (! args(7).is_scalar_type () || ! args(7).isnumeric ()
      || args(7).iscomplex ())
  {
    error ("gamboostinter: MaxNumSplits must be a numeric scalar.");
  }
  octave_idx_type maxsplits = (octave_idx_type) args(7).int_value ();
  if (maxsplits < 1)
  {
    error ("gamboostinter: MaxNumSplits must be a positive integer.");
  }

  if (method == 1)
  {
    for (octave_idx_type i = 0; i < Y.numel (); i++)
    {
      if (Y(i) != 0.0 && Y(i) != 1.0)
      {
        error ("gamboostinter: Y must hold zeros and ones for Method 1.");
      }
    }
  }

  GamInterFit F = gamb_boost_inter (X, Y, F0, method, pairs, maxtrees,
                                    lrate, maxsplits);

  Cell edges (1, d);
  for (octave_idx_type j = 0; j < d; j++)
  {
    edges(j) = octave_value (F.edges[(std::size_t) j]);
  }

  octave_idx_type np = (octave_idx_type) F.term.size ();
  Cell values (1, np);
  for (octave_idx_type q = 0; q < np; q++)
  {
    values(q) = octave_value (F.term[(std::size_t) q].value);
  }

  octave_scalar_map Mdl;
  Mdl.assign ("PairBinEdges", octave_value (edges));
  Mdl.assign ("PairValues", octave_value (values));
  Mdl.assign ("Intercept", octave_value (F.shift));
  Mdl.assign ("NumTrees", octave_value ((double) F.ntrees));
  Mdl.assign ("ReasonForTermination", octave_value (F.reason));
  Mdl.assign ("Deviance", octave_value (F.deviance));
  Mdl.assign ("Residuals", octave_value (F.residuals));

  return ovl (Mdl);
}

/*
%!test
%! ## Every selected pair gets a surface on the coarse grid of its two
%! ## predictors, and the phase reports the usual fields.
%! x = randn (200, 3);
%! y = double (x(:,1) .* x(:,2) + 0.1 * randn (200, 1) > 0);
%! M = gamboosttrain (x, y, 1, 20, 1, 1);
%! f0 = gamboostpredict (M.BinEdges, M.ShapeValues, x, M.Intercept);
%! I = gamboostinter (x, y, f0, 1, [1, 2], 20, 1, 4);
%! assert_equal (numel (I.PairValues), 1);
%! assert_equal (size (I.PairValues{1}), [8, 8]);
%! assert_equal (numel (I.PairBinEdges), 3);
%! assert_equal (numel (I.Residuals), 200);

%!test
%! ## The phase starts from the deviance the predictor phase left, so it can
%! ## only improve on it.
%! x = randn (200, 3);
%! y = double (x(:,1) .* x(:,2) + 0.1 * randn (200, 1) > 0);
%! M = gamboosttrain (x, y, 1, 20, 1, 1);
%! f0 = gamboostpredict (M.BinEdges, M.ShapeValues, x, M.Intercept);
%! I = gamboostinter (x, y, f0, 1, [1, 2], 30, 1, 4);
%! assert_equal (I.Deviance <= M.Deviance, true);

%!test
%! ## A pair term is recentred like a shape function, so what it gives up is
%! ## reported for the intercept rather than left inside the surface.
%! x = randn (150, 2);
%! y = double (x(:,1) .* x(:,2) > 0);
%! M = gamboosttrain (x, y, 1, 20, 1, 1);
%! f0 = gamboostpredict (M.BinEdges, M.ShapeValues, x, M.Intercept);
%! I = gamboostinter (x, y, f0, 1, [1, 2], 25, 1, 4);
%! assert_equal (isfinite (I.Intercept), true);

%!test
%! ## Several pairs are fitted in one call, each on its own grid.
%! x = randn (200, 4);
%! y = double (x(:,1) .* x(:,2) + x(:,3) .* x(:,4) > 0);
%! M = gamboosttrain (x, y, 1, 20, 1, 1);
%! f0 = gamboostpredict (M.BinEdges, M.ShapeValues, x, M.Intercept);
%! I = gamboostinter (x, y, f0, 1, [1, 2; 3, 4], 20, 1, 4);
%! assert_equal (numel (I.PairValues), 2);

%!test
%! ## The budget is a budget here too: a phase that stops improving says so.
%! x = randn (120, 2);
%! y = double (x(:,1) > 0);
%! M = gamboosttrain (x, y, 1, 100, 1, 1);
%! f0 = gamboostpredict (M.BinEdges, M.ShapeValues, x, M.Intercept);
%! I = gamboostinter (x, y, f0, 1, [1, 2], 100000, 1, 4);
%! assert_equal (I.ReasonForTermination, 'Unable to improve the model fit.');
%! assert_equal (I.NumTrees < 100000, true);

%!test
%! ## A regression pair term works the same way and lowers the residual sum
%! ## of squares it was handed.
%! x = randn (200, 2);
%! y = x(:,1) .* x(:,2);
%! M = gamboosttrain (x, y, 2, 20, 1, 1);
%! f0 = gamboostpredict (M.BinEdges, M.ShapeValues, x, M.Intercept);
%! I = gamboostinter (x, y, f0, 2, [1, 2], 40, 1, 4);
%! assert_equal (I.Deviance < M.Deviance, true);

%!error<Invalid call> gamboostinter (1, 2, 3)
%!error<gamboostinter: X must be a numeric matrix.> ...
%! gamboostinter ('a', [1;0], [0;0], 1, [1,2], 10, 1, 4)
%!error<gamboostinter: X, Y and F0 must have the same number of rows.> ...
%! gamboostinter ([1,2;3,4], [1;0], [0;0;0], 1, [1,2], 10, 1, 4)
%!error<gamboostinter: Method must be 1 or 2.> ...
%! gamboostinter ([1,2;3,4], [1;0], [0;0], 3, [1,2], 10, 1, 4)
%!error<gamboostinter: Pairs must be a numeric matrix with two columns.> ...
%! gamboostinter ([1,2;3,4], [1;0], [0;0], 1, [1,2,3], 10, 1, 4)
%!error<gamboostinter: Pairs must hold predictor indices between 1 and 2.> ...
%! gamboostinter ([1,2;3,4], [1;0], [0;0], 1, [1,5], 10, 1, 4)
%!error<gamboostinter: a pair must name two different predictors.> ...
%! gamboostinter ([1,2;3,4], [1;0], [0;0], 1, [2,2], 10, 1, 4)
%!error<gamboostinter: NumTrees must be a positive integer.> ...
%! gamboostinter ([1,2;3,4], [1;0], [0;0], 1, [1,2], 0, 1, 4)
%!error<gamboostinter: LearnRate must be greater than 0 and at most 1.> ...
%! gamboostinter ([1,2;3,4], [1;0], [0;0], 1, [1,2], 10, 2, 4)
%!error<gamboostinter: MaxNumSplits must be a positive integer.> ...
%! gamboostinter ([1,2;3,4], [1;0], [0;0], 1, [1,2], 10, 1, 0)
%!error<gamboostinter: Y must hold zeros and ones for Method 1.> ...
%! gamboostinter ([1,2;3,4], [1;2], [0;0], 1, [1,2], 10, 1, 4)
*/
