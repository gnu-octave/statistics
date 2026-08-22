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

#include <cmath>
#include <vector>
#include <string>
#include <octave/oct.h>
#include <octave/ov-struct.h>
#include "gam.cpp"

using namespace std;

DEFUN_DLD(gampredict, args, ,
          "-*- texinfo -*-\n\
 @deftypefn  {statistics} {@var{yFit} =} gampredict (@var{Parameters}, @\n\
 @var{X}, @var{Intercept})\n\
 @deftypefnx {statistics} {@var{score} =} gampredict (@var{Parameters}, @\n\
 @var{X}, @var{Intercept}, @var{Link})\n\
\n\
\n\
Evaluate a generalized additive model on new data.\n\
\n\
@code{@var{yFit} = gampredict (@var{Parameters}, @var{X}, @var{Intercept})}\n\
adds the intercept to the sum of the additive terms evaluated at each row of\n\
@var{X} and returns the @math{Nx1} result.  It is the prediction engine\n\
shared by @code{ClassificationGAM} and @code{RegressionGAM}, and it is not\n\
meant to be called directly.\n\
\n\
@itemize\n\
@item @var{Parameters} is a @math{1xP} structure array of piecewise\n\
polynomials, as returned by @code{gamtrain} in the field of the same name.\n\
\n\
@item @var{X} is an @math{NxP} numeric matrix with one column per additive\n\
term, and a count that does not match is an error.  A model carrying\n\
interaction terms must therefore be given the augmented matrix, not the\n\
predictors alone.\n\
A missing value predicts @qcode{NaN}, since no term of the model is defined\n\
at it.  A value outside the range the term was fitted over is extrapolated\n\
from the nearest piece, as @code{ppval} extrapolates.\n\
\n\
@item @var{Intercept} is the model's constant term.\n\
\n\
@item @var{Link}, if given, selects what the additive prediction is mapped\n\
through: @qcode{0} returns it as it stands and @qcode{1} takes it as a\n\
log-odds, returning the @math{Nx2} matrix of class probabilities whose\n\
second column is the logistic function of it.  The default is @qcode{0}.\n\
@end itemize\n\
\n\
@seealso{gamtrain, ClassificationGAM, RegressionGAM, fitcgam, fitrgam}\n\
@end deftypefn")
{
  octave_idx_type nargin = args.length ();
  if (nargin < 3 || nargin > 4)
  {
    print_usage ();
  }

  if (! args(0).isstruct () || args(0).isempty ())
  {
    error ("gampredict: Parameters must be a structure array.");
  }
  octave_map params = args(0).map_value ();
  if (! params.isfield ("breaks") || ! params.isfield ("coefs"))
  {
    error ("gampredict: Parameters must have 'breaks' and 'coefs' fields.");
  }

  if (! args(1).isnumeric () || args(1).iscomplex () || args(1).isempty ())
  {
    error ("gampredict: X must be a numeric matrix.");
  }
  Matrix X = args(1).matrix_value ();

  octave_idx_type n = X.rows ();
  octave_idx_type d = X.columns ();
  // One column per additive term, exactly.  The m-code this replaces looped
  // over the columns of X and took the term count from them, which tolerated
  // a caller passing a matrix that had not been augmented with its
  // interaction terms: too many terms went unnoticed and evaluated a
  // truncated model, while too few raised an out-of-bound index.  That
  // contract was not self-consistent, and the direction it tolerated is the
  // one that returns a wrong number instead of an error.
  if (params.numel () != d)
  {
    error ("gampredict: X must have one column per additive term.");
  }

  if (! args(2).is_scalar_type () || ! args(2).isnumeric ()
      || args(2).iscomplex ())
  {
    error ("gampredict: Intercept must be a numeric scalar.");
  }
  double intercept = args(2).scalar_value ();

  int link = 0;
  if (nargin > 3)
  {
    if (! args(3).is_scalar_type () || ! args(3).isnumeric ()
        || args(3).iscomplex ())
    {
      error ("gampredict: Link must be a numeric scalar.");
    }
    link = args(3).int_value ();
    if (link != 0 && link != 1)
    {
      error ("gampredict: Link must be either 0 or 1.");
    }
  }

  Cell c_breaks = params.contents ("breaks");
  Cell c_coefs = params.contents ("coefs");

  ColumnVector y (n, intercept);
  for (octave_idx_type j = 0; j < d; j++)
  {
    if (! c_breaks(j).isnumeric () || ! c_coefs(j).isnumeric ())
    {
      error ("gampredict: 'breaks' and 'coefs' must be numeric.");
    }
    RowVector breaks = c_breaks(j).row_vector_value ();
    Matrix coefs = c_coefs(j).matrix_value ();
    if (breaks.numel () != coefs.rows () + 1)
    {
      error ("gampredict: term %d has %d breaks for %d pieces.",
             (int) j + 1, (int) breaks.numel (), (int) coefs.rows ());
    }
    for (octave_idx_type i = 0; i < n; i++)
    {
      y(i) += gam_ppval (breaks, coefs, X(i, j));
    }
  }

  if (link == 0)
  {
    return ovl (y);
  }

  Matrix score (n, 2);
  for (octave_idx_type i = 0; i < n; i++)
  {
    double pos = 1.0 / (1.0 + exp (-y(i)));
    score(i, 0) = 1.0 - pos;
    score(i, 1) = pos;
  }

  return ovl (score);
}

/*
%!test
%! ## The additive prediction is the intercept plus each term's spline
%! x = linspace (0, 1, 40)';
%! pp = splinefit (x, cos (3*x), 5, 'order', 3);
%! P = struct ('form', 'pp', 'breaks', pp.breaks, 'coefs', pp.coefs, ...
%!             'pieces', pp.pieces, 'order', pp.order, 'dim', pp.dim);
%! assert_equal (gampredict (P, x, 0), ppval (pp, x), 1e-14);
%! assert_equal (gampredict (P, x, 2.5), ppval (pp, x) + 2.5, 1e-14);

%!test
%! ## Two terms add
%! x = linspace (0, 1, 30)';
%! p1 = splinefit (x, cos (3*x), 5, 'order', 3);
%! p2 = splinefit (x, x.^2, 4, 'order', 2);
%! P(1) = struct ('form', 'pp', 'breaks', p1.breaks, 'coefs', p1.coefs, ...
%!                'pieces', p1.pieces, 'order', p1.order, 'dim', p1.dim);
%! P(2) = struct ('form', 'pp', 'breaks', p2.breaks, 'coefs', p2.coefs, ...
%!                'pieces', p2.pieces, 'order', p2.order, 'dim', p2.dim);
%! y = gampredict (P, [x, x], 1);
%! assert_equal (y, 1 + ppval (p1, x) + ppval (p2, x), 1e-13);

%!test
%! ## Fewer columns than terms is an error, not a truncated model.  A caller
%! ## that passes the predictors alone to a model carrying interaction terms
%! ## used to be answered with the prediction of its leading terms.
%! x = linspace (0, 1, 30)';
%! p1 = splinefit (x, cos (3*x), 5, 'order', 3);
%! p2 = splinefit (x, x.^2, 4, 'order', 2);
%! P2(1) = struct ('form', 'pp', 'breaks', p1.breaks, 'coefs', p1.coefs, ...
%!                 'pieces', p1.pieces, 'order', p1.order, 'dim', p1.dim);
%! P2(2) = struct ('form', 'pp', 'breaks', p2.breaks, 'coefs', p2.coefs, ...
%!                 'pieces', p2.pieces, 'order', p2.order, 'dim', p2.dim);
%! fail ('gampredict (P2, x, 0)', ...
%!       'gampredict: X must have one column per additive term.');
%! assert_equal (gampredict (P2, [x, x], 0), ...
%!               ppval (p1, x) + ppval (p2, x), 1e-13);

%!test
%! ## The logistic link returns both class probabilities, and they sum to one
%! x = linspace (0, 1, 20)';
%! pp = splinefit (x, 4*x - 2, 5, 'order', 3);
%! P = struct ('form', 'pp', 'breaks', pp.breaks, 'coefs', pp.coefs, ...
%!             'pieces', pp.pieces, 'order', pp.order, 'dim', pp.dim);
%! s = gampredict (P, x, 0, 1);
%! assert_equal (size (s), [20, 2]);
%! assert_equal (sum (s, 2), ones (20, 1), 1e-14);
%! assert_equal (s(:,2), 1 ./ (1 + exp (- gampredict (P, x, 0))), 1e-14);

%!test
%! ## A point outside the fitted range is extrapolated from the nearest piece,
%! ## exactly as ppval extrapolates
%! x = linspace (0, 1, 30)';
%! pp = splinefit (x, cos (3*x), 5, 'order', 3);
%! P = struct ('form', 'pp', 'breaks', pp.breaks, 'coefs', pp.coefs, ...
%!             'pieces', pp.pieces, 'order', pp.order, 'dim', pp.dim);
%! xq = [-0.5; 1.5];
%! assert_equal (gampredict (P, xq, 0), ppval (pp, xq), 1e-13);

%!test
%! ## A missing predictor predicts NaN, and leaves the other rows alone
%! x = linspace (0, 1, 10)';
%! pp = splinefit (x, cos (3*x), 4, 'order', 3);
%! P = struct ('form', 'pp', 'breaks', pp.breaks, 'coefs', pp.coefs, ...
%!             'pieces', pp.pieces, 'order', pp.order, 'dim', pp.dim);
%! xq = [0.2; NaN; 0.8];
%! y = gampredict (P, xq, 0);
%! assert (isnan (y(2)));
%! assert_equal (y([1, 3]), ppval (pp, [0.2; 0.8]), 1e-14);

%!shared P
%! pp = splinefit (linspace (0, 1, 20)', linspace (0, 1, 20)', 4, 'order', 3);
%! P = struct ('form', 'pp', 'breaks', pp.breaks, 'coefs', pp.coefs, ...
%!             'pieces', pp.pieces, 'order', pp.order, 'dim', pp.dim);
%!error<Invalid call> gampredict ()
%!error<Invalid call> gampredict (P, ones (5, 1))
%!error<Invalid call> gampredict (P, ones (5, 1), 0, 1, 2)
%!error<gampredict: Parameters must be a structure array.> ...
%! gampredict (5, ones (5, 1), 0)
%!error<gampredict: X must be a numeric matrix.> gampredict (P, 'a', 0)
%!error<gampredict: X must have one column per additive term.> ...
%! gampredict (P, ones (5, 3), 0)
%!error<gampredict: Intercept must be a numeric scalar.> ...
%! gampredict (P, ones (5, 1), [1, 2])
%!error<gampredict: Link must be a numeric scalar.> ...
%! gampredict (P, ones (5, 1), 0, 'a')
%!error<gampredict: Link must be either 0 or 1.> ...
%! gampredict (P, ones (5, 1), 0, 2)
*/
