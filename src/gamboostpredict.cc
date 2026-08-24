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

DEFUN_DLD(gamboostpredict, args, ,
          "-*- texinfo -*-\n\
@deftypefn  {statistics} {@var{Y} =} gamboostpredict (@var{BinEdges}, @\n\
@var{ShapeValues}, @var{X}, @var{Intercept})\n\
@deftypefnx {statistics} {@var{Y} =} gamboostpredict (@dots{}, @var{Link})\n\
@deftypefnx {statistics} {@var{Y} =} gamboostpredict (@dots{}, @var{Link}, @\n\
@var{PairBinEdges}, @var{PairValues}, @var{Pairs})\n\
\n\
Predict from a generalized additive model of boosted trees.\n\
\n\
@itemize\n\
@item @var{BinEdges} is a @math{1xP} cell of row vectors and\n\
@var{ShapeValues} a @math{1xP} cell of column vectors, as\n\
@code{gamboosttrain} returns them in the fields of the same names.  Each\n\
shape function is a step function over its predictor's bins, so a term is\n\
evaluated by a lookup.\n\
\n\
@item @var{X} is an @math{NxP} numeric matrix with one column per additive\n\
term, and a count that does not match is an error.  A missing value\n\
contributes nothing from that term rather than making the whole prediction\n\
@qcode{NaN}, which is what a tree does with a value it cannot place.  A\n\
value outside the range the term was fitted over falls in the nearest bin,\n\
so a shape function is constant beyond its data rather than extrapolated.\n\
\n\
@item @var{Intercept} is the model's constant term.\n\
\n\
@item @var{Link}, if given, selects what the additive prediction is mapped\n\
through: @qcode{0} returns it as it stands and @qcode{1} takes it as a\n\
log-odds, returning the @math{Nx2} matrix of class probabilities whose\n\
second column is the logistic function of it.  The default is @qcode{0}.\n\
\n\
@item @var{PairBinEdges}, @var{PairValues} and @var{Pairs} carry the\n\
interaction terms, as @code{gamboostinter} returns the first two and as it\n\
was given the third.  Each pair contributes the value of the cell its two\n\
predictors fall in, and a pair with either predictor missing contributes\n\
nothing.  All three must be given together or none of them.\n\
@end itemize\n\
\n\
@seealso{gamboosttrain, gampredict, ClassificationGAM, RegressionGAM}\n\
@end deftypefn")
{
  octave_idx_type nargin = args.length ();
  if (nargin != 4 && nargin != 5 && nargin != 8)
  {
    print_usage ();
  }

  if (! args(0).iscell ())
  {
    error ("gamboostpredict: BinEdges must be a cell array.");
  }
  if (! args(1).iscell ())
  {
    error ("gamboostpredict: ShapeValues must be a cell array.");
  }

  Cell edges = args(0).cell_value ();
  Cell values = args(1).cell_value ();

  if (edges.numel () != values.numel ())
  {
    error ("gamboostpredict: BinEdges and ShapeValues must have the same "
           "number of elements.");
  }

  if (! args(2).isnumeric () || args(2).iscomplex () || args(2).isempty ())
  {
    error ("gamboostpredict: X must be a numeric matrix.");
  }
  Matrix X = args(2).matrix_value ();

  if (X.columns () != edges.numel ())
  {
    error ("gamboostpredict: X must have one column per additive term.");
  }

  if (! args(3).is_scalar_type () || ! args(3).isnumeric ()
      || args(3).iscomplex ())
  {
    error ("gamboostpredict: Intercept must be a numeric scalar.");
  }
  double intercept = args(3).scalar_value ();

  int link = 0;
  if (nargin > 4)
  {
    if (! args(4).is_scalar_type () || ! args(4).isnumeric ()
        || args(4).iscomplex ())
    {
      error ("gamboostpredict: Link must be a numeric scalar.");
    }
    link = args(4).int_value ();
    if (link != 0 && link != 1)
    {
      error ("gamboostpredict: Link must be 0 or 1.");
    }
  }

  GamBoostFit F;
  F.intercept = intercept;
  octave_idx_type d = edges.numel ();
  F.edges.resize ((std::size_t) d);
  F.value.resize ((std::size_t) d);

  for (octave_idx_type j = 0; j < d; j++)
  {
    if (! edges(j).isnumeric () || edges(j).iscomplex ())
    {
      error ("gamboostpredict: every BinEdges element must be numeric.");
    }
    if (! values(j).isnumeric () || values(j).iscomplex ())
    {
      error ("gamboostpredict: every ShapeValues element must be numeric.");
    }
    F.edges[(std::size_t) j] = edges(j).row_vector_value ();
    F.value[(std::size_t) j] = values(j).column_vector_value ();

    if (F.value[(std::size_t) j].numel ()
        != F.edges[(std::size_t) j].numel () + 1)
    {
      error ("gamboostpredict: term %d has %d values for %d cut points.",
             (int) j + 1, (int) F.value[(std::size_t) j].numel (),
             (int) F.edges[(std::size_t) j].numel ());
    }
  }

  octave_idx_type n = X.rows ();
  ColumnVector y (n);
  for (octave_idx_type i = 0; i < n; i++)
  {
    y(i) = gamb_predict_row (F, X, i);
  }

  if (nargin == 8)
  {
    if (! args(5).iscell () || ! args(6).iscell ())
    {
      error ("gamboostpredict: PairBinEdges and PairValues must be cell "
             "arrays.");
    }
    Cell pedges = args(5).cell_value ();
    Cell pvals = args(6).cell_value ();

    if (pedges.numel () != d)
    {
      error ("gamboostpredict: PairBinEdges must have one element per "
             "additive term.");
    }

    if (! args(7).isnumeric () || args(7).iscomplex ()
        || args(7).columns () != 2)
    {
      error ("gamboostpredict: Pairs must be a numeric matrix with two "
             "columns.");
    }
    Matrix pairs = args(7).matrix_value ();

    if (pairs.rows () != pvals.numel ())
    {
      error ("gamboostpredict: Pairs and PairValues must have the same "
             "number of rows.");
    }

    std::vector<RowVector> pe ((std::size_t) d);
    for (octave_idx_type j = 0; j < d; j++)
    {
      pe[(std::size_t) j] = pedges(j).row_vector_value ();
    }

    for (octave_idx_type q = 0; q < pairs.rows (); q++)
    {
      octave_idx_type j = (octave_idx_type) pairs(q, 0) - 1;
      octave_idx_type k = (octave_idx_type) pairs(q, 1) - 1;
      if (j < 0 || k < 0 || j >= d || k >= d)
      {
        error ("gamboostpredict: Pairs must hold predictor indices between "
               "1 and %d.", (int) d);
      }
      Matrix V = pvals(q).matrix_value ();

      for (octave_idx_type i = 0; i < n; i++)
      {
        octave_idx_type a = gamb_bin_of (pe[(std::size_t) j], X(i, j));
        octave_idx_type b = gamb_bin_of (pe[(std::size_t) k], X(i, k));
        if (a < 0 || b < 0)
        {
          continue;
        }
        if (a >= V.rows () || b >= V.columns ())
        {
          error ("gamboostpredict: pair %d has a %dx%d surface, too small "
                 "for its grid.", (int) q + 1, (int) V.rows (),
                 (int) V.columns ());
        }
        y(i) += V(a, b);
      }
    }
  }

  if (link == 0)
  {
    return ovl (y);
  }

  Matrix score (n, 2);
  for (octave_idx_type i = 0; i < n; i++)
  {
    double pos = 1.0 / (1.0 + std::exp (-y(i)));
    score(i, 0) = 1.0 - pos;
    score(i, 1) = pos;
  }

  return ovl (score);
}

/*
%!test
%! ## The additive prediction is the intercept plus each term's step function.
%! E = {[1.5, 2.5]};
%! V = {[-1; 0; 2]};
%! assert_equal (gamboostpredict (E, V, [1; 2; 3], 0), [-1; 0; 2], 1e-14);
%! assert_equal (gamboostpredict (E, V, [1; 2; 3], 0.5), [-0.5; 0.5; 2.5], ...
%!               1e-14);

%!test
%! ## Terms add.
%! E = {[1.5], [10]};
%! V = {[1; 2], [4; 8]};
%! assert_equal (gamboostpredict (E, V, [1, 5; 2, 20], 0), [5; 10], 1e-14);

%!test
%! ## A value beyond the fitted range falls in the nearest bin, so a shape
%! ## function is constant outside its data rather than extrapolated.
%! E = {[1.5, 2.5]};
%! V = {[-1; 0; 2]};
%! assert_equal (gamboostpredict (E, V, [-100; 100], 0), [-1; 2], 1e-14);

%!test
%! ## A missing value costs that term only.
%! E = {[1.5], [10]};
%! V = {[1; 2], [4; 8]};
%! assert_equal (gamboostpredict (E, V, [NaN, 5], 0), 4, 1e-14);

%!test
%! ## A term with no cut points is a constant.
%! E = {zeros(1,0)};
%! V = {3};
%! assert_equal (gamboostpredict (E, V, [1; 2; 3], 1), [4; 4; 4], 1e-14);

%!test
%! ## The logistic link returns both class probabilities, and they sum to one.
%! E = {[1.5]};
%! V = {[-2; 2]};
%! s = gamboostpredict (E, V, [1; 2], 0, 1);
%! assert_equal (sum (s, 2), [1; 1], 1e-14);
%! assert_equal (s(:,2), 1 ./ (1 + exp ([2; -2])), 1e-14);

%!test
%! ## A model round-trips through its own trainer.
%! x = [1; 2; 3; 4; 5; 6; 7; 8];
%! y = [0; 0; 0; 0; 1; 1; 1; 1];
%! M = gamboosttrain (x, y, 1, 30, 1, 1);
%! f = gamboostpredict (M.BinEdges, M.ShapeValues, x, M.Intercept);
%! p = 1 ./ (1 + exp (-f));
%! assert_equal (M.Residuals, y - p, 1e-12);

%!test
%! ## A pair term adds the value of the cell its two predictors fall in.
%! E = {zeros(1,0), zeros(1,0)};
%! V = {0, 0};
%! PE = {[1.5], [10]};
%! PV = {[1, 2; 3, 4]};
%! X = [1, 5; 1, 20; 2, 5; 2, 20];
%! assert_equal (gamboostpredict (E, V, X, 0, 0, PE, PV, [1, 2]), ...
%!               [1; 2; 3; 4], 1e-14);

%!test
%! ## Interaction terms add to the additive part rather than replacing it.
%! E = {[1.5], zeros(1,0)};
%! V = {[10; 20], 0};
%! PE = {[1.5], [10]};
%! PV = {[1, 2; 3, 4]};
%! X = [1, 5; 2, 20];
%! assert_equal (gamboostpredict (E, V, X, 100, 0, PE, PV, [1, 2]), ...
%!               [111; 124], 1e-14);

%!test
%! ## A pair with either predictor missing contributes nothing, as a tree does
%! ## with a value it cannot place.
%! E = {zeros(1,0), zeros(1,0)};
%! V = {0, 0};
%! PE = {[1.5], [10]};
%! PV = {[1, 2; 3, 4]};
%! assert_equal (gamboostpredict (E, V, [NaN, 5], 0, 0, PE, PV, [1, 2]), 0);

%!test
%! ## The assembled model reproduces what the interaction phase itself
%! ## computed: same residuals, so the stored surfaces and the recentring
%! ## bookkeeping are consistent with the fit they came from.
%! randn ('seed', 42);
%! rand ('seed', 42);
%! X = randn (200, 3);
%! Y = double (rand (200, 1) < 1 ./ (1 + exp (-2 * X(:,1) .* X(:,2))));
%! M = gamboosttrain (X, Y, 1, 50, 1, 1);
%! f0 = gamboostpredict (M.BinEdges, M.ShapeValues, X, M.Intercept);
%! I = gamboostinter (X, Y, f0, 1, [1, 2], 40, 1, 4);
%! f1 = gamboostpredict (M.BinEdges, M.ShapeValues, X, ...
%!                       M.Intercept + I.Intercept, 0, I.PairBinEdges, ...
%!                       I.PairValues, [1, 2]);
%! assert_equal (Y - 1 ./ (1 + exp (-f1)), I.Residuals, 1e-12);

%!error<Invalid call> gamboostpredict ({1}, {1})
%!error<gamboostpredict: BinEdges must be a cell array.> ...
%! gamboostpredict (1, {1}, [1;2], 0)
%!error<gamboostpredict: ShapeValues must be a cell array.> ...
%! gamboostpredict ({1}, 1, [1;2], 0)
%!error<gamboostpredict: BinEdges and ShapeValues must have the same number of elements.> ...
%! gamboostpredict ({1, 2}, {1}, [1;2], 0)
%!error<gamboostpredict: X must be a numeric matrix.> ...
%! gamboostpredict ({[1.5]}, {[1;2]}, 'a', 0)
%!error<gamboostpredict: X must have one column per additive term.> ...
%! gamboostpredict ({[1.5]}, {[1;2]}, [1, 2], 0)
%!error<gamboostpredict: Intercept must be a numeric scalar.> ...
%! gamboostpredict ({[1.5]}, {[1;2]}, [1;2], [1, 2])
%!error<gamboostpredict: Link must be 0 or 1.> ...
%! gamboostpredict ({[1.5]}, {[1;2]}, [1;2], 0, 2)
%!error<gamboostpredict: term 1 has 3 values for 1 cut points.> ...
%! gamboostpredict ({[1.5]}, {[1;2;3]}, [1;2], 0)
%!error<gamboostpredict: PairBinEdges and PairValues must be cell arrays.> ...
%! gamboostpredict ({[1.5], [1.5]}, {[1;2], [1;2]}, [1, 1], 0, 0, 1, {1}, [1, 2])
%!error<gamboostpredict: PairBinEdges must have one element per additive term.> ...
%! gamboostpredict ({[1.5], [1.5]}, {[1;2], [1;2]}, [1, 1], 0, 0, {[1.5]}, ...
%!                  {[1, 2; 3, 4]}, [1, 2])
%!error<gamboostpredict: Pairs must be a numeric matrix with two columns.> ...
%! gamboostpredict ({[1.5], [1.5]}, {[1;2], [1;2]}, [1, 1], 0, 0, ...
%!                  {[1.5], [1.5]}, {[1, 2; 3, 4]}, [1, 2, 3])
%!error<gamboostpredict: Pairs and PairValues must have the same number of rows.> ...
%! gamboostpredict ({[1.5], [1.5]}, {[1;2], [1;2]}, [1, 1], 0, 0, ...
%!                  {[1.5], [1.5]}, {[1, 2; 3, 4]}, [1, 2; 1, 2])
*/
