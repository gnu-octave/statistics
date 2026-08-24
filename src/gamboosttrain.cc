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

DEFUN_DLD(gamboosttrain, args, ,
          "-*- texinfo -*-\n\
@deftypefn  {statistics} {@var{Mdl} =} gamboosttrain (@var{X}, @var{Y}, @\n\
@var{Method}, @var{NumTrees}, @var{LearnRate}, @var{MaxNumSplits})\n\
@deftypefnx {statistics} {@var{Mdl} =} gamboosttrain (@dots{}, @var{Verbose}, @\n\
@var{NumPrint})\n\
\n\
Fit a generalized additive model of boosted trees.\n\
\n\
@code{@var{Mdl} = gamboosttrain (@var{X}, @var{Y}, @var{Method}, @\n\
@var{NumTrees}, @var{LearnRate}, @var{MaxNumSplits})} boosts one tree per\n\
column of @var{X} in each round and returns the additive model as a\n\
structure.  It is the fitting engine shared by @code{ClassificationGAM} and\n\
@code{RegressionGAM}, and it is not meant to be called directly.\n\
\n\
@itemize\n\
@item @var{X} is an @math{NxP} numeric matrix of predictors.  A missing\n\
value is not an error: the observation takes no part in the affected\n\
predictor's trees and that term contributes nothing to its prediction.\n\
\n\
@item @var{Y} is an @math{Nx1} numeric vector of responses.  For\n\
@var{Method} 1 it must hold zeros and ones.\n\
\n\
@item @var{Method} selects what is boosted: @qcode{1} the logistic deviance,\n\
as a classifier is fitted, and @qcode{2} the squared error, as a regression\n\
is fitted.\n\
\n\
@item @var{NumTrees} is the number of rounds, each fitting one tree per\n\
predictor.  It is a budget rather than a count: a fit that stops improving\n\
ends earlier and says so.\n\
\n\
@item @var{LearnRate} is the step a round starts at.  A round that fails to\n\
earn its place is retried at half the step, so this is an initial value and\n\
not a fixed one.\n\
\n\
@item @var{MaxNumSplits} is the largest number of splits any one tree may\n\
make.  @qcode{1} is a stump.\n\
\n\
@item @var{Verbose}, if greater than zero, prints a trace of the fit, and\n\
@var{NumPrint} how often: the first round and then every @var{NumPrint}\n\
rounds.  The @qcode{RelTol} column is the relative improvement the round\n\
bought, which is what the stopping rule reads.  MATLAB prints a column under\n\
the same heading holding a quantity of its own that cannot be derived from\n\
anything else it reports, so the two are not comparable.\n\
@end itemize\n\
\n\
@var{Mdl} is a structure with the following fields.\n\
\n\
@itemize\n\
@item @qcode{Intercept}, the constant term the additive terms are added to.\n\
For a classifier it is fitted rather than fixed: it is seeded with the\n\
log-odds of the response mean and then collects the constant each shape\n\
function gives up when it is recentred.  For a regression it is the response\n\
mean and stays there.\n\
@item @qcode{BinEdges}, a @math{1xP} cell of row vectors, the cut points each\n\
predictor was binned at.\n\
@item @qcode{ShapeValues}, a @math{1xP} cell of column vectors, one value per\n\
bin.  A shape function is a step function, so this is the whole of it however\n\
many trees produced it.\n\
@item @qcode{NumTrees}, the number of rounds actually performed.\n\
@item @qcode{ReasonForTermination}, why fitting stopped.\n\
@item @qcode{Deviance}, the deviance at the last round.\n\
@item @qcode{Residuals}, the @math{Nx1} residual vector at the last round.\n\
@end itemize\n\
\n\
@seealso{gamboostpredict, gamtrain, ClassificationGAM, RegressionGAM}\n\
@end deftypefn")
{
  octave_idx_type nargin = args.length ();
  if (nargin != 6 && nargin != 8)
  {
    print_usage ();
  }

  if (! args(0).isnumeric () || args(0).iscomplex () || args(0).isempty ())
  {
    error ("gamboosttrain: X must be a numeric matrix.");
  }
  if (! args(1).isnumeric () || args(1).iscomplex () || args(1).isempty ()
      || args(1).columns () != 1)
  {
    error ("gamboosttrain: Y must be a numeric column vector.");
  }
  if (args(0).rows () != args(1).rows ())
  {
    error ("gamboosttrain: X and Y must have the same number of rows.");
  }

  Matrix X = args(0).matrix_value ();
  ColumnVector Y = args(1).column_vector_value ();
  octave_idx_type d = X.columns ();

  if (! args(2).is_scalar_type () || ! args(2).isnumeric ()
      || args(2).iscomplex ())
  {
    error ("gamboosttrain: Method must be a numeric scalar.");
  }
  int method = args(2).int_value ();
  if (method != 1 && method != 2)
  {
    error ("gamboosttrain: Method must be 1 or 2.");
  }

  if (! args(3).is_scalar_type () || ! args(3).isnumeric ()
      || args(3).iscomplex ())
  {
    error ("gamboosttrain: NumTrees must be a numeric scalar.");
  }
  octave_idx_type maxtrees = (octave_idx_type) args(3).int_value ();
  if (maxtrees < 1)
  {
    error ("gamboosttrain: NumTrees must be a positive integer.");
  }

  if (! args(4).is_scalar_type () || ! args(4).isnumeric ()
      || args(4).iscomplex ())
  {
    error ("gamboosttrain: LearnRate must be a numeric scalar.");
  }
  double lrate = args(4).scalar_value ();
  if (! (lrate > 0.0) || lrate > 1.0)
  {
    error ("gamboosttrain: LearnRate must be greater than 0 and at most 1.");
  }

  if (! args(5).is_scalar_type () || ! args(5).isnumeric ()
      || args(5).iscomplex ())
  {
    error ("gamboosttrain: MaxNumSplits must be a numeric scalar.");
  }
  octave_idx_type maxsplits = (octave_idx_type) args(5).int_value ();
  if (maxsplits < 1)
  {
    error ("gamboosttrain: MaxNumSplits must be a positive integer.");
  }

  if (method == 1)
  {
    for (octave_idx_type i = 0; i < Y.numel (); i++)
    {
      if (Y(i) != 0.0 && Y(i) != 1.0)
      {
        error ("gamboosttrain: Y must hold zeros and ones for Method 1.");
      }
    }
  }

  int verbose = 0;
  octave_idx_type numprint = 10;
  if (nargin == 8)
  {
    if (! args(6).is_scalar_type () || ! args(6).isnumeric ()
        || args(6).iscomplex ())
    {
      error ("gamboosttrain: Verbose must be a numeric scalar.");
    }
    verbose = args(6).int_value ();
    if (verbose < 0)
    {
      error ("gamboosttrain: Verbose must be 0 or greater.");
    }
    if (! args(7).is_scalar_type () || ! args(7).isnumeric ()
        || args(7).iscomplex ())
    {
      error ("gamboosttrain: NumPrint must be a numeric scalar.");
    }
    numprint = (octave_idx_type) args(7).int_value ();
    if (numprint < 1)
    {
      error ("gamboosttrain: NumPrint must be a positive integer.");
    }
  }

  GamBoostFit F = gamb_boost (X, Y, method, maxtrees, lrate, maxsplits,
                              verbose, numprint);

  Cell edges (1, d);
  Cell values (1, d);
  for (octave_idx_type j = 0; j < d; j++)
  {
    edges(j) = octave_value (F.edges[(std::size_t) j]);
    values(j) = octave_value (F.value[(std::size_t) j]);
  }

  octave_scalar_map Mdl;
  Mdl.assign ("Intercept", octave_value (F.intercept));
  Mdl.assign ("BinEdges", octave_value (edges));
  Mdl.assign ("ShapeValues", octave_value (values));
  Mdl.assign ("NumTrees", octave_value ((double) F.ntrees));
  Mdl.assign ("ReasonForTermination", octave_value (F.reason));
  Mdl.assign ("Deviance", octave_value (F.deviance));
  Mdl.assign ("Residuals", octave_value (F.residuals));

  return ovl (Mdl);
}

/*
%!test
%! ## Every field is present and the shapes follow the predictor count.
%! x = [1 5; 2 4; 3 3; 4 2; 5 1; 6 7; 7 8; 8 9; 9 10; 10 11];
%! y = [0; 0; 0; 0; 0; 1; 1; 1; 1; 1];
%! M = gamboosttrain (x, y, 1, 50, 1, 1);
%! assert_equal (isstruct (M), true);
%! assert_equal (numel (M.BinEdges), 2);
%! assert_equal (numel (M.ShapeValues), 2);
%! assert_equal (numel (M.ShapeValues{1}), numel (M.BinEdges{1}) + 1);
%! assert_equal (numel (M.Residuals), 10);

%!test
%! ## A predictor is cut at every midpoint between its distinct values.
%! x = [1; 2; 3; 4; 2; 3];
%! y = [0; 0; 1; 1; 0; 1];
%! M = gamboosttrain (x, y, 1, 5, 1, 1);
%! assert_equal (M.BinEdges{1}, [1.5, 2.5, 3.5], 1e-12);

%!test
%! ## A constant predictor admits no split: one bin and no cut points.
%! x = [2; 2; 2; 2];
%! y = [0; 1; 0; 1];
%! M = gamboosttrain (x, y, 1, 5, 1, 1);
%! assert_equal (isempty (M.BinEdges{1}), true);
%! assert_equal (numel (M.ShapeValues{1}), 1);

%!test
%! ## A regression intercept is the response mean and boosting leaves it there.
%! x = [1; 2; 3; 4; 5; 6; 7; 8];
%! y = [2; 4; 5; 4; 6; 8; 9; 10];
%! M = gamboosttrain (x, y, 2, 40, 1, 1);
%! assert_equal (M.Intercept, mean (y), 1e-12);

%!test
%! ## A classifier intercept is seeded with the log-odds of the response mean
%! ## and moves away from it, so a balanced response does not stay at zero.
%! x = [1; 2; 3; 4; 5; 6; 7; 8];
%! y = [0; 0; 0; 0; 1; 1; 1; 1];
%! M = gamboosttrain (x, y, 1, 60, 1, 1);
%! assert_equal (M.Intercept != 0, true);

%!test
%! ## The tree budget is a budget: a fit that converges reports so and stops
%! ## short of it.
%! x = [1; 2; 3; 4; 5; 6; 7; 8];
%! y = [0; 0; 0; 0; 1; 1; 1; 1];
%! M = gamboosttrain (x, y, 1, 5000, 1, 1);
%! assert_equal (M.ReasonForTermination, 'Unable to improve the model fit.');
%! assert_equal (M.NumTrees < 5000, true);

%!test
%! ## A budget too small to converge in reports the other reason.
%! x = [1; 2; 3; 4; 5; 6; 7; 8];
%! y = [0; 0; 0; 0; 1; 1; 1; 1];
%! M = gamboosttrain (x, y, 1, 2, 1, 1);
%! assert_equal (M.ReasonForTermination, ...
%!               'Terminated after training the requested number of trees.');
%! assert_equal (M.NumTrees, 2);

%!test
%! ## A single-class response has an infinite log-odds and nothing to fit.
%! x = [1; 2; 3; 4];
%! y = [1; 1; 1; 1];
%! M = gamboosttrain (x, y, 1, 10, 1, 1);
%! assert_equal (M.Intercept, Inf);
%! assert_equal (M.NumTrees, 0);

%!test
%! ## A missing predictor value costs the observation that term, not the fit.
%! x = [1; 2; NaN; 4; 5; 6];
%! y = [0; 0; 0; 1; 1; 1];
%! M = gamboosttrain (x, y, 1, 20, 1, 1);
%! assert_equal (isfinite (M.Intercept), true);
%! assert_equal (any (isnan (M.ShapeValues{1})), false);

%!test
%! ## More splits per tree reach a lower deviance on a shape one split cannot
%! ## follow.
%! x = [1; 2; 3; 4; 5; 6; 7; 8; 9; 10];
%! y = [0; 0; 1; 1; 1; 1; 1; 1; 0; 0];
%! M1 = gamboosttrain (x, y, 1, 30, 1, 1);
%! M2 = gamboosttrain (x, y, 1, 30, 1, 3);
%! assert_equal (M2.Deviance < M1.Deviance, true);

%!test
%! ## Above the cap the cuts are equally spaced through the OBSERVATIONS, not
%! ## through the distinct values, so they crowd where the data is.  Values 1
%! ## to 50 carry 1000 of these 1450 rows and values 51 to 500 carry one row
%! ## each: the equal-frequency grid puts 58 cuts in that dense fifth of the
%! ## range where a grid spread through the distinct values would put 25.
%! x = [repmat((1:50)', 20, 1); (51:500)'];
%! y = double (mod (1:1450, 2))';
%! M = gamboosttrain (x, y, 1, 1, 1, 1);
%! assert_equal (sum (M.BinEdges{1} < 50), 58);

%!test
%! ## Ties collapse cuts: two quantile positions inside one repeated value
%! ## give the same cut, and the repeats are dropped, so a tied predictor ends
%! ## with fewer than the cap allows rather than with duplicate edges.
%! x = [repmat((1:50)', 20, 1); (51:500)'];
%! y = double (mod (1:1450, 2))';
%! M = gamboosttrain (x, y, 1, 1, 1, 1);
%! assert_equal (numel (M.BinEdges{1}), 138);
%! assert_equal (numel (unique (M.BinEdges{1})), 138);

%!test
%! ## The cap binds at 255 cut points however many distinct values there are.
%! x = (1:2000)';
%! y = double (mod (1:2000, 2))';
%! M = gamboosttrain (x, y, 1, 1, 1, 1);
%! assert_equal (numel (M.BinEdges{1}), 255);

%!test
%! ## Below the cap nothing is thinned: one cut per gap between distinct
%! ## values, which is what MATLAB, scikit-learn and the EBM all report.
%! x = (1:200)';
%! y = double (mod (1:200, 2))';
%! M = gamboosttrain (x, y, 1, 1, 1, 1);
%! assert_equal (numel (M.BinEdges{1}), 199);
%! assert_equal (M.BinEdges{1}(1:3), [1.5, 2.5, 3.5], 1e-12);

%!test
%! ## Patience looks across a window rather than at one round: a fit that has
%! ## stopped earning its keep ends and says so, well short of its budget.
%! x = [1; 2; 3; 4; 5; 6; 7; 8];
%! y = [0; 0; 0; 0; 1; 1; 1; 1];
%! M = gamboosttrain (x, y, 1, 100000, 1, 1);
%! assert_equal (M.ReasonForTermination, 'Unable to improve the model fit.');
%! assert_equal (M.NumTrees < 1000, true);

%!test
%! ## A verbose fit prints a trace and returns the same model as a quiet one.
%! x = [1; 2; 3; 4; 5; 6; 7; 8];
%! y = [0; 0; 0; 0; 1; 1; 1; 1];
%! Q = gamboosttrain (x, y, 1, 20, 1, 1);
%! V = evalc ('W = gamboosttrain (x, y, 1, 20, 1, 1, 1, 5);');
%! assert_equal (W.Intercept, Q.Intercept, 1e-12);
%! assert_equal (W.NumTrees, Q.NumTrees);
%! assert_equal (! isempty (strfind (V, 'NumTrees')), true);
%! assert_equal (! isempty (strfind (V, 'LearnRate')), true);

%!test
%! ## NumPrint controls how often a round is reported: the first, then every
%! ## NumPrint after it.
%! x = [1; 2; 3; 4; 5; 6; 7; 8];
%! y = [0; 0; 0; 0; 1; 1; 1; 1];
%! V5 = evalc ('gamboosttrain (x, y, 1, 20, 1, 1, 1, 5);');
%! V1 = evalc ('gamboosttrain (x, y, 1, 20, 1, 1, 1, 1);');
%! assert_equal (numel (strfind (V1, '|    1D|')) > ...
%!               numel (strfind (V5, '|    1D|')), true);

%!test
%! ## Verbose 0 prints nothing at all.
%! x = [1; 2; 3; 4; 5; 6; 7; 8];
%! y = [0; 0; 0; 0; 1; 1; 1; 1];
%! V = evalc ('gamboosttrain (x, y, 1, 20, 1, 1, 0, 5);');
%! assert_equal (isempty (strfind (V, '1D')), true);

%!error<Invalid call> gamboosttrain (1, 2, 3)
%!error<gamboosttrain: X must be a numeric matrix.> ...
%! gamboosttrain ('a', [1;0], 1, 10, 1, 1)
%!error<gamboosttrain: Y must be a numeric column vector.> ...
%! gamboosttrain ([1;2], [1, 0], 1, 10, 1, 1)
%!error<gamboosttrain: X and Y must have the same number of rows.> ...
%! gamboosttrain ([1;2;3], [1;0], 1, 10, 1, 1)
%!error<gamboosttrain: Method must be 1 or 2.> ...
%! gamboosttrain ([1;2], [1;0], 3, 10, 1, 1)
%!error<gamboosttrain: NumTrees must be a positive integer.> ...
%! gamboosttrain ([1;2], [1;0], 1, 0, 1, 1)
%!error<gamboosttrain: LearnRate must be greater than 0 and at most 1.> ...
%! gamboosttrain ([1;2], [1;0], 1, 10, 0, 1)
%!error<gamboosttrain: LearnRate must be greater than 0 and at most 1.> ...
%! gamboosttrain ([1;2], [1;0], 1, 10, 2, 1)
%!error<gamboosttrain: MaxNumSplits must be a positive integer.> ...
%! gamboosttrain ([1;2], [1;0], 1, 10, 1, 0)
%!error<gamboosttrain: Y must hold zeros and ones for Method 1.> ...
%! gamboosttrain ([1;2], [1;2], 1, 10, 1, 1)
%!error<gamboosttrain: Verbose must be 0 or greater.> ...
%! gamboosttrain ([1;2], [1;0], 1, 10, 1, 1, -1, 5)
%!error<gamboosttrain: NumPrint must be a positive integer.> ...
%! gamboosttrain ([1;2], [1;0], 1, 10, 1, 1, 1, 0)
*/
