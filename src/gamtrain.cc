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

DEFUN_DLD(gamtrain, args, ,
          "-*- texinfo -*-\n\
 @deftypefn  {statistics} {@var{Mdl} =} gamtrain (@var{X}, @var{Y}, @\n\
 @var{Knots}, @var{Order}, @var{Method}, @var{Inter}, @var{P1}, @var{P2})\n\
\n\
\n\
Fit a generalized additive model of smoothing splines.\n\
\n\
@code{@var{Mdl} = gamtrain (@var{X}, @var{Y}, @var{Knots}, @var{Order}, @\n\
@var{Method}, @var{Inter}, @var{P1}, @var{P2})} fits one univariate spline\n\
per column of @var{X} and returns the additive model as a structure.  It is\n\
the fitting engine shared by @code{ClassificationGAM} and\n\
@code{RegressionGAM}, and it is not meant to be called directly.\n\
\n\
@itemize\n\
@item @var{X} is an @math{NxP} numeric matrix of predictors.  A missing\n\
value is not an error: the observation is left out of the affected\n\
predictor's spline and its prediction from that term is @qcode{NaN}.\n\
\n\
@item @var{Y} is an @math{Nx1} numeric vector of responses.  For\n\
@var{Method} 1 it must hold zeros and ones.\n\
\n\
@item @var{Knots} is a @math{1xP} vector giving the number of spline pieces\n\
for each predictor, and @var{Order} a @math{1xP} vector giving the degree of\n\
the polynomial on each piece.  A spline of @math{K} pieces and degree\n\
@math{D} spans a space of @math{K + D} dimensions.\n\
\n\
@item @var{Method} selects the fitting scheme: @qcode{1} boosts the log-odds\n\
by gradient descent, as a classifier is fitted, and @qcode{2} backfits the\n\
partial residuals, as a regression is fitted.\n\
\n\
@item @var{Inter} is the intercept the fit starts from: a proportion for\n\
@var{Method} 1, which is stored as its log-odds, and the response mean for\n\
@var{Method} 2.  A proportion of zero or one is not an error: its log-odds is\n\
infinite, the gradient is zero throughout and every additive term stays at\n\
zero, which is the fit a single-class response has.\n\
\n\
@item @var{P1} and @var{P2} are the scheme's two parameters.  For\n\
@var{Method} 1 they are the learning rate and the number of boosting\n\
iterations; for @var{Method} 2 the convergence tolerance and the maximum\n\
number of backfitting cycles.\n\
@end itemize\n\
\n\
@var{Mdl} is a structure with the following fields.\n\
\n\
@itemize\n\
@item @qcode{Intercept}, the constant term the additive terms are added to.\n\
@item @qcode{Parameters}, a @math{1xP} structure array of piecewise\n\
polynomials in the form @code{ppval} consumes, one per predictor.\n\
@item @qcode{Iterations}, the number of iterations performed.\n\
@item @qcode{Residuals}, the @math{Nx1} residual vector at the last\n\
iteration.\n\
@item @qcode{RSS}, the scalar residual sum of squares for @var{Method} 1 and\n\
the @math{1xP} per-term criterion the backfitting stops on for @var{Method}\n\
2.\n\
@end itemize\n\
\n\
@seealso{gampredict, ClassificationGAM, RegressionGAM, fitcgam, fitrgam}\n\
@end deftypefn")
{
  if (args.length () != 8)
  {
    print_usage ();
  }

  // X and Y
  if (! args(0).isnumeric () || args(0).iscomplex () || args(0).isempty ())
  {
    error ("gamtrain: X must be a numeric matrix.");
  }
  if (! args(1).isnumeric () || args(1).iscomplex () || args(1).isempty ()
      || args(1).columns () != 1)
  {
    error ("gamtrain: Y must be a numeric column vector.");
  }
  if (args(0).rows () != args(1).rows ())
  {
    error ("gamtrain: X and Y must have the same number of rows.");
  }

  Matrix X = args(0).matrix_value ();
  ColumnVector Y = args(1).column_vector_value ();
  octave_idx_type d = X.columns ();

  // Knots and Order, one per predictor
  if (! args(2).isnumeric () || args(2).iscomplex ()
      || args(2).numel () != d)
  {
    error ("gamtrain: Knots must have one element per column of X.");
  }
  if (! args(3).isnumeric () || args(3).iscomplex ()
      || args(3).numel () != d)
  {
    error ("gamtrain: Order must have one element per column of X.");
  }

  Matrix knots = args(2).matrix_value ();
  Matrix order = args(3).matrix_value ();
  for (octave_idx_type j = 0; j < d; j++)
  {
    if (knots(j) < 1 || floor (knots(j)) != knots(j))
    {
      error ("gamtrain: Knots must be positive integers.");
    }
    // Order zero is the m-code's own lower bound and gives piecewise
    // constants: splinefit takes the degree, so a degree of zero is a spline
    // of order one.  Rejecting it here would narrow what the learners accept.
    if (order(j) < 0 || floor (order(j)) != order(j))
    {
      error ("gamtrain: Order must be non-negative integers.");
    }
  }

  // Method, intercept and the scheme's two parameters
  if (! args(4).is_scalar_type () || ! args(4).isnumeric ()
      || args(4).iscomplex ())
  {
    error ("gamtrain: Method must be a numeric scalar.");
  }
  int method = args(4).int_value ();
  if (method != 1 && method != 2)
  {
    error ("gamtrain: Method must be either 1 or 2.");
  }

  if (! args(5).is_scalar_type () || ! args(5).isnumeric ()
      || args(5).iscomplex ())
  {
    error ("gamtrain: Inter must be a numeric scalar.");
  }
  double inter = args(5).scalar_value ();

  if (! args(6).is_scalar_type () || ! args(6).isnumeric ()
      || args(6).iscomplex () || args(6).scalar_value () <= 0)
  {
    error ("gamtrain: P1 must be a positive scalar.");
  }
  if (! args(7).is_scalar_type () || ! args(7).isnumeric ()
      || args(7).iscomplex () || args(7).scalar_value () < 1)
  {
    error ("gamtrain: P2 must be a scalar not less than 1.");
  }
  double p1 = args(6).scalar_value ();
  octave_idx_type p2 = (octave_idx_type) args(7).scalar_value ();

  GamFit fit;
  if (method == 1)
  {
    // Inter is not restricted to the open unit interval.  A single-class fit
    // is allowed and gives a proportion of exactly zero or one, whose log-odds
    // is infinite; the gradient is then zero everywhere and every additive
    // term stays at zero, which is the answer such a fit has.
    fit = gam_boost (X, Y, inter, knots, order, p1, p2);
  }
  else
  {
    fit = gam_backfit (X, Y, inter, knots, order, p1, p2);
  }

  octave_scalar_map Mdl;
  Mdl.assign ("Intercept", fit.intercept);
  Mdl.assign ("Parameters", fit.params);
  Mdl.assign ("Iterations", fit.iterations);
  Mdl.assign ("Residuals", fit.res);
  Mdl.assign ("RSS", fit.RSS);

  return ovl (Mdl);
}

/*
%!test
%! ## The boosted model returns the fields the learners store
%! X = [linspace(0, 1, 40)', linspace(1, 2, 40)'];
%! Y = double ([1:40]' > 20);
%! Mdl = gamtrain (X, Y, [5, 5], [3, 3], 1, 0.5, 0.1, 100);
%! assert_equal (fieldnames (Mdl), ...
%!               {'Intercept'; 'Parameters'; 'Iterations'; 'Residuals'; 'RSS'});
%! assert_equal (size (Mdl.Parameters), [1, 2]);
%! assert_equal (Mdl.Iterations, 100);
%! assert_equal (size (Mdl.Residuals), [40, 1]);
%! assert_equal (size (Mdl.RSS), [1, 1]);
%! assert_equal (Mdl.Intercept, log (0.5 / 0.5), 1e-14);

%!test
%! ## A spline of K pieces and degree D is a piecewise polynomial of K pieces
%! ## and K + D coefficients per term
%! X = linspace (0, 1, 30)';
%! Y = double ([1:30]' > 15);
%! Mdl = gamtrain (X, Y, 6, 3, 1, 0.5, 0.1, 20);
%! P = Mdl.Parameters(1);
%! assert_equal (P.form, 'pp');
%! assert_equal (P.pieces, 6);
%! assert_equal (P.order, 4);
%! assert_equal (P.dim, 1);
%! assert_equal (size (P.coefs), [6, 4]);
%! assert_equal (size (P.breaks), [1, 7]);

%!test
%! ## The engine reproduces splinefit on a single fitted term.  One boosting
%! ## round at a learning rate of one is the spline of the gradient, which for
%! ## an intercept of one half is the response less one half.
%! x = linspace (0, 1, 50)';
%! Y = double (x > 0.4);
%! Mdl = gamtrain (x, Y, 5, 3, 1, 0.5, 1, 1);
%! pp = splinefit (x, Y - 0.5, 5, 'order', 3);
%! assert_equal (Mdl.Parameters(1).coefs, pp.coefs, 1e-9);
%! assert_equal (Mdl.Parameters(1).breaks, pp.breaks, 1e-14);

%!test
%! ## Backfitting a single term converges to the spline of the centred
%! ## response, which is what one cycle already fits
%! x = linspace (0, 2*pi, 60)';
%! Y = cos (x);
%! Mdl = gamtrain (x, Y, 5, 3, 2, mean (Y), 1e-3, 1000);
%! pp = splinefit (x, Y - mean (Y), 5, 'order', 3);
%! assert_equal (Mdl.Parameters(1).coefs, pp.coefs, 1e-9);
%! assert_equal (Mdl.Intercept, mean (Y), 1e-14);
%! assert_equal (size (Mdl.RSS), [1, 1]);

%!test
%! ## Backfitting stops on the tolerance, and a looser one stops sooner
%! X = [linspace(0, 1, 50)', linspace(0, 2, 50)'];
%! Y = cos (3 * X(:,1)) + X(:,2);
%! M1 = gamtrain (X, Y, [5, 5], [3, 3], 2, mean (Y), 1e-12, 1000);
%! M2 = gamtrain (X, Y, [5, 5], [3, 3], 2, mean (Y), 1, 1000);
%! assert (M2.Iterations <= M1.Iterations);
%! assert_equal (size (M1.RSS), [1, 2]);

%!test
%! ## Fewer observations than the spline space has dimensions: the fit is the
%! ## minimum norm one and does not raise
%! x = [1; 2; 3; 4; 5];
%! Y = [2; 1; 4; 3; 6];
%! Mdl = gamtrain (x, Y, 5, 3, 2, mean (Y), 1e-3, 1000);
%! assert_equal (size (Mdl.Parameters(1).coefs), [5, 4]);
%! assert (all (isfinite (Mdl.Parameters(1).coefs(:))));

%!test
%! ## A missing predictor drops the observation from that term's fit and
%! ## predicts NaN for it, leaving the other observations finite
%! x = linspace (0, 1, 40)';
%! x(7) = NaN;
%! Y = double ([1:40]' > 20);
%! Mdl = gamtrain (x, Y, 5, 3, 1, 0.5, 0.1, 10);
%! assert (isnan (Mdl.Residuals(7)));
%! assert (all (isfinite (Mdl.Residuals([1:6, 8:40]))));

%!test
%! ## Tied and unsorted predictor values are handled by the break placement
%! x = [0.5; 0.5; 0.1; 0.9; 0.1; 0.7; 0.3; 0.5; 0.2; 0.6; 0.4; 0.8];
%! Y = [0; 0; 0; 1; 0; 1; 0; 1; 0; 1; 0; 1];
%! Mdl = gamtrain (x, Y, 5, 3, 1, 0.5, 0.1, 20);
%! assert (all (isfinite (Mdl.Parameters(1).coefs(:))));
%! assert (issorted (Mdl.Parameters(1).breaks));

%!error<Invalid call> gamtrain ()
%!error<Invalid call> gamtrain (ones (5, 2), ones (5, 1))
%!error<gamtrain: X must be a numeric matrix.> ...
%! gamtrain ('a', ones (5, 1), 5, 3, 1, 0.5, 0.1, 10)
%!error<gamtrain: Y must be a numeric column vector.> ...
%! gamtrain (ones (5, 1), ones (5, 2), 5, 3, 1, 0.5, 0.1, 10)
%!error<gamtrain: X and Y must have the same number of rows.> ...
%! gamtrain (ones (5, 1), ones (4, 1), 5, 3, 1, 0.5, 0.1, 10)
%!error<gamtrain: Knots must have one element per column of X.> ...
%! gamtrain (ones (5, 2), ones (5, 1), 5, [3, 3], 1, 0.5, 0.1, 10)
%!error<gamtrain: Order must have one element per column of X.> ...
%! gamtrain (ones (5, 2), ones (5, 1), [5, 5], 3, 1, 0.5, 0.1, 10)
%!error<gamtrain: Knots must be positive integers.> ...
%! gamtrain (ones (5, 1), ones (5, 1), 2.5, 3, 1, 0.5, 0.1, 10)
%!error<gamtrain: Order must be non-negative integers.> ...
%! gamtrain (ones (5, 1), ones (5, 1), 5, -1, 1, 0.5, 0.1, 10)
%!error<gamtrain: Order must be non-negative integers.> ...
%! gamtrain (ones (5, 1), ones (5, 1), 5, 1.5, 1, 0.5, 0.1, 10)
%!error<gamtrain: Method must be a numeric scalar.> ...
%! gamtrain (ones (5, 1), ones (5, 1), 5, 3, 'a', 0.5, 0.1, 10)
%!error<gamtrain: Method must be either 1 or 2.> ...
%! gamtrain (ones (5, 1), ones (5, 1), 5, 3, 3, 0.5, 0.1, 10)
%!error<gamtrain: Inter must be a numeric scalar.> ...
%! gamtrain (ones (5, 1), ones (5, 1), 5, 3, 1, [1, 2], 0.1, 10)
%!error<gamtrain: P1 must be a positive scalar.> ...
%! gamtrain (ones (5, 1), ones (5, 1), 5, 3, 1, 0.5, 0, 10)
%!error<gamtrain: P2 must be a scalar not less than 1.> ...
%! gamtrain (ones (5, 1), ones (5, 1), 5, 3, 1, 0.5, 0.1, 0)

%!test
%! ## Order zero is piecewise constants, which the learners accept and the
%! ## m-code fitted as splines of order one
%! x = linspace (0, 1, 30)';
%! Y = double (x > 0.5);
%! Mdl = gamtrain (x, Y, 5, 0, 1, 0.5, 0.1, 20);
%! P = Mdl.Parameters(1);
%! assert_equal (P.pieces, 5);
%! assert_equal (P.order, 1);
%! assert_equal (size (P.coefs), [5, 1]);
%! pp = splinefit (x, Y - 0.5, 5, 'order', 0);
%! M1 = gamtrain (x, Y, 5, 0, 1, 0.5, 1, 1);
%! assert_equal (M1.Parameters(1).coefs, pp.coefs, 1e-12);

%!test
%! ## A single-class response gives an infinite intercept and terms that stay
%! ## at zero, which is what the m-code it replaces produced
%! x = linspace (0, 1, 20)';
%! Mdl = gamtrain (x, ones (20, 1), 5, 3, 1, 1, 0.1, 10);
%! assert_equal (Mdl.Intercept, Inf);
%! assert_equal (Mdl.Parameters(1).coefs, zeros (5, 4));
%! assert_equal (Mdl.Residuals, zeros (20, 1));
*/
