## Copyright (C) 2026 Andreas Bertsatos <abertsatos@biol.uoa.gr>
##
## This file is part of the statistics package for GNU Octave.
##
## Octave is free software; you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or
## (at your option) any later version.
##
## Octave is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with Octave; see the file COPYING.  If not,
## see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn  {statistics} {@var{x0} =} invpred (@var{x}, @var{y}, @var{y0})
## @deftypefnx {statistics} {[@var{x0}, @var{dxlo}, @var{dxup}] =} invpred (@var{x}, @var{y}, @var{y0})
## @deftypefnx {statistics} {[@dots{}] =} invpred (@dots{}, @var{name}, @var{value})
##
## Inverse prediction from a simple linear regression.
##
## @code{@var{x0} = invpred (@var{x}, @var{y}, @var{y0})} fits the simple
## linear regression of @var{y} on @var{x} and returns, for each element of
## @var{y0}, the value of the predictor at which the fitted line takes that
## response.  @var{x} and @var{y} must be vectors of real values of the same
## length; @var{y0} may be of any size and @var{x0} is returned with the same
## size.  Observations where either @var{x} or @var{y} is @code{NaN} are
## dropped in pairs before the fit.
##
## @code{[@var{x0}, @var{dxlo}, @var{dxup}] = invpred (@dots{})} also returns
## the width of a confidence interval on either side of @var{x0}, so that the
## interval is @code{[@var{x0} - @var{dxlo}, @var{x0} + @var{dxup}]}.  The
## bounds follow Fieller's theorem and are therefore not symmetric about
## @var{x0}.  They are not simultaneous over the elements of @var{y0}, and
## they need not be finite: when the slope is not significantly different
## from zero at the requested level, the interval either covers the entire
## real line, with @var{dxlo} and @var{dxup} both @code{Inf}, or it excludes
## a finite range around the fitted line and extends to @code{Inf} on one
## side only, with the other side of @var{dxlo}/@var{dxup} finite.
##
## @code{[@dots{}] = invpred (@dots{}, @var{name}, @var{value})} accepts the
## following name-value pairs:
##
## @itemize
## @item
## @qcode{"alpha"} is the significance level of the interval, a scalar
## strictly between 0 and 1, so that the interval has confidence
## @math{100 * (1 - @var{alpha})%}.  The default is @qcode{0.05}.
##
## @item
## @qcode{"predopt"} selects what the interval covers.  With
## @qcode{"observation"}, the default, it covers a new observation whose
## response is @var{y0}.  With @qcode{"curve"}, it covers the point at which
## the true regression line takes the value @var{y0}, and is narrower because
## it carries no new-observation variance.
## @end itemize
##
## @seealso{regress, fitlm, polyfit}
## @end deftypefn

function [x0, dxlo, dxup] = invpred (x, y, y0, varargin)

  ## Check for valid number of input arguments
  if (nargin < 3)
    print_usage ();
  endif

  if (! (isvector (x) && isnumeric (x) && isreal (x)))
    error ("invpred: X must be a vector of real values.");
  endif
  if (! (isvector (y) && isnumeric (y) && isreal (y)))
    error ("invpred: Y must be a vector of real values.");
  endif
  if (numel (x) != numel (y))
    error ("invpred: X and Y must have the same length.");
  endif
  if (! (isnumeric (y0) && isreal (y0)))
    error ("invpred: Y0 must be numeric and real.");
  endif

  ## Parse optional name-value pairs
  alpha = 0.05;
  predopt = "observation";
  if (mod (numel (varargin), 2) != 0)
    error ("invpred: optional arguments must be name-value pairs.");
  endif
  for i = 1:2:numel (varargin)
    name = varargin{i};
    value = varargin{i+1};
    if (! (ischar (name) && isrow (name)))
      error ("invpred: parameter name must be a character vector.");
    endif
    switch (lower (name))
      case 'alpha'
        alpha = value;
      case 'predopt'
        predopt = value;
      otherwise
        error ("invpred: invalid parameter name: %s.", name);
    endswitch
  endfor
  if (! (isnumeric (alpha) && isreal (alpha) && isscalar (alpha) ...
         && alpha > 0 && alpha < 1))
    error ("invpred: ALPHA must be a scalar strictly between 0 and 1.");
  endif
  if (! (ischar (predopt) && isrow (predopt)) ...
      || ! any (strcmpi (predopt, {"curve", "observation"})))
    error ("invpred: PREDOPT value must be 'curve' or 'observation'.");
  endif

  ## Drop observations that are missing in either variable
  x = x(:);
  y = y(:);
  keep = ! (isnan (x) | isnan (y));
  x = x(keep);
  y = y(keep);
  n = numel (x);

  xbar = mean (x);
  ybar = mean (y);
  sxx = sum ((x - xbar) .^ 2);
  if (sxx == 0)
    error ("invpred: cannot compute inverse predictions if X is constant.");
  endif

  ## Fit the line and take the residual variance about it
  slope = sum ((x - xbar) .* (y - ybar)) / sxx;
  s2 = sum ((y - (ybar + slope * (x - xbar))) .^ 2) / (n - 2);

  ## The point estimate simply inverts the fitted line
  offset = (y0 - ybar) / slope;
  x0 = xbar + offset;

  if (nargout < 2)
    return;
  endif

  ## Fieller's interval.  The factor g compares the sampling variability of
  ## the slope with the slope itself; once it reaches 1 the line cannot be
  ## distinguished from a horizontal one and the interval is unbounded.
  t_crit = tinv (1 - alpha / 2, n - 2);
  g = t_crit ^ 2 * s2 / (slope ^ 2 * sxx);
  if (g == 1)
    dxlo = Inf (size (y0));
    dxup = Inf (size (y0));
    return;
  endif

  ## A new observation carries its own variance; a point on the curve does not
  extra = double (strcmpi (predopt, "observation"));
  centre = offset / (1 - g);
  term = offset .^ 2 / sxx + (1 - g) * (1 / n + extra);
  halfwidth = (t_crit * sqrt (s2) / (abs (slope) * abs (1 - g))) ...
              .* sqrt (max (term, 0));

  if (g < 1)
    dxlo = offset - centre + halfwidth;
    dxup = centre + halfwidth - offset;
  else
    ## Outside the significance level for the slope: the confidence region
    ## is the exterior of the two roots, so one side of x0 is unbounded.
    dxlo = Inf (size (y0));
    dxup = Inf (size (y0));
    real_roots = term >= 0;
    lower_ray = real_roots & (offset < centre);
    upper_ray = real_roots & (offset >= centre);
    dxup(lower_ray) = centre - halfwidth - offset(lower_ray);
    dxlo(upper_ray) = offset(upper_ray) - centre - halfwidth;
  endif

endfunction

%!demo
%! ## Estimate the predictor value at which a fitted line reaches a response
%! ## of 20, with a 95% confidence interval either side of it.
%!
%! x = (1:10)';
%! y = 2 + 3 * x + [0.5; -0.3; 0.2; 0.8; -0.6; 0.1; -0.4; 0.7; -0.2; 0.3];
%! [x0, dxlo, dxup] = invpred (x, y, 20)
%! interval = [x0 - dxlo, x0 + dxup]

## The point estimate inverts the fitted line.
%!test
%! x = (1:10)';
%! y = 2 + 3 * x + [0.5; -0.3; 0.2; 0.8; -0.6; 0.1; -0.4; 0.7; -0.2; 0.3];
%! x0 = invpred (x, y, 20);
%! assert_equal (x0, 5.964741641337386, 1e-12);

## A confidence interval is returned on either side of the estimate.
%!test
%! x = (1:10)';
%! y = 2 + 3 * x + [0.5; -0.3; 0.2; 0.8; -0.6; 0.1; -0.4; 0.7; -0.2; 0.3];
%! [x0, dxlo, dxup] = invpred (x, y, 20);
%! assert_equal (x0, 5.964741641337386, 1e-12);
%! assert_equal (dxlo, 0.408566704275405, 1e-12);
%! assert_equal (dxup, 0.410279499114103, 1e-12);

## The interval is not symmetric about the estimate.
%!test
%! x = (1:10)';
%! y = 2 + 3 * x + [0.5; -0.3; 0.2; 0.8; -0.6; 0.1; -0.4; 0.7; -0.2; 0.3];
%! [~, dxlo, dxup] = invpred (x, y, 20);
%! assert_equal (dxlo != dxup, true);

## Y0 may hold several values, and the output follows its shape.
%!test
%! x = (1:10)';
%! y = 2 + 3 * x + [0.5; -0.3; 0.2; 0.8; -0.6; 0.1; -0.4; 0.7; -0.2; 0.3];
%! [x0, dxlo, dxup] = invpred (x, y, [10; 20; 30]);
%! assert_equal (x0, [2.621276595744681; 5.964741641337386; ...
%!                    9.308206686930092], 1e-12);
%! assert_equal (dxlo, [0.432537159571881; 0.408566704275405; ...
%!                      0.433439039304474], 1e-12);
%! assert_equal (dxup, [0.421927689383959; 0.410279499114103; ...
%!                      0.447474099169791], 1e-12);
%! [x0r, dxlor, dxupr] = invpred (x, y, [10, 20, 30]);
%! assert_equal (x0r, x0');
%! assert_equal (dxlor, dxlo');
%! assert_equal (dxupr, dxup');

## A matrix of responses is returned as a matrix.
%!test
%! x = (1:10)';
%! y = 2 + 3 * x + [0.5; -0.3; 0.2; 0.8; -0.6; 0.1; -0.4; 0.7; -0.2; 0.3];
%! x0 = invpred (x, y, [10, 20; 30, 40]);
%! assert_equal (size (x0), [2, 2]);
%! assert_equal (x0, [2.621276595744681, 5.964741641337386; ...
%!                    9.308206686930092, 12.651671732522798], 1e-12);

## An empty Y0 gives an empty result.
%!test
%! x = (1:10)';
%! y = 2 + 3 * x + [0.5; -0.3; 0.2; 0.8; -0.6; 0.1; -0.4; 0.7; -0.2; 0.3];
%! x0 = invpred (x, y, []);
%! assert_equal (size (x0), [0, 0]);

## ALPHA widens the interval without moving the estimate.
%!test
%! x = (1:10)';
%! y = 2 + 3 * x + [0.5; -0.3; 0.2; 0.8; -0.6; 0.1; -0.4; 0.7; -0.2; 0.3];
%! [x0, dxlo, dxup] = invpred (x, y, 20, 'alpha', 0.01);
%! assert_equal (x0, 5.964741641337386, 1e-12);
%! assert_equal (dxlo, 0.594536204654897, 1e-12);
%! assert_equal (dxup, 0.598170042325828, 1e-12);

## An interval on the curve is narrower than one on a new observation.
%!test
%! x = (1:10)';
%! y = 2 + 3 * x + [0.5; -0.3; 0.2; 0.8; -0.6; 0.1; -0.4; 0.7; -0.2; 0.3];
%! [x0, dxlo, dxup] = invpred (x, y, 20, 'predopt', 'curve');
%! assert_equal (x0, 5.964741641337386, 1e-12);
%! assert_equal (dxlo, 0.124048892442359, 1e-12);
%! assert_equal (dxup, 0.125761687281057, 1e-12);
%! [~, obslo, obsup] = invpred (x, y, 20, 'predopt', 'observation');
%! assert_equal (dxlo < obslo && dxup < obsup, true);

## A slope indistinguishable from zero leaves the interval unbounded.
%!test
%! [x0, dxlo, dxup] = invpred ((1:5)', [1; 5; 2; 8; 3], 4);
%! assert_equal (x0, 3.285714285714286, 1e-12);
%! assert_equal (dxlo, Inf);
%! assert_equal (dxup, Inf);

## Rows and columns are accepted alike.
%!test
%! x = (1:10)';
%! y = 2 + 3 * x + [0.5; -0.3; 0.2; 0.8; -0.6; 0.1; -0.4; 0.7; -0.2; 0.3];
%! [x0, dxlo, dxup] = invpred (x', y', 20);
%! assert_equal (x0, 5.964741641337386, 1e-12);
%! assert_equal (dxlo, 0.408566704275405, 1e-12);
%! assert_equal (dxup, 0.410279499114103, 1e-12);

## Observations missing in either variable are dropped in pairs.
%!test
%! x = (1:10)';
%! y = 2 + 3 * x + [0.5; -0.3; 0.2; 0.8; -0.6; 0.1; -0.4; 0.7; -0.2; 0.3];
%! xn = x; xn(3) = NaN;
%! assert_equal (invpred (xn, y, 20), 5.967084254482929, 1e-12);
%! yn = y; yn(4) = NaN;
%! assert_equal (invpred (x, yn, 20), 5.988352745424295, 1e-12);

## When the slope is not significant, one side of the interval may be
## unbounded while the other stays finite.
%!test
%! [x0, dxlo, dxup] = invpred ((1:5)', [1; 5; 2; 8; 3], 100);
%! assert_equal (x0, 140.42857142857144, 1e-12);
%! assert_equal (dxlo, 111.30773209684287, 1e-12);
%! assert_equal (dxup, Inf);

%!error<Invalid call to invpred> invpred ((1:10)', (1:10)')
%!error<invpred: X must be a vector of real values.> ...
%! invpred (ones (10, 2), (1:10)', 5)
%!error<invpred: Y must be a vector of real values.> ...
%! invpred ((1:10)', ones (10, 2), 5)
%!error<invpred: X and Y must have the same length.> ...
%! invpred ((1:10)', (1:5)', 5)
%!error<invpred: ALPHA must be a scalar strictly between 0 and 1.> ...
%! invpred ((1:10)', (1:10)', 5, 'alpha', 5)
%!error<invpred: ALPHA must be a scalar strictly between 0 and 1.> ...
%! invpred ((1:10)', (1:10)', 5, 'alpha', [0.05, 0.01])
%!error<invpred: PREDOPT value must be 'curve' or 'observation'.> ...
%! invpred ((1:10)', (1:10)', 5, 'predopt', 'bogus')
%!error<invpred: invalid parameter name: badname.> ...
%! invpred ((1:10)', (1:10)', 5, 'badname', 1)
%!error<invpred: optional arguments must be name-value pairs.> ...
%! invpred ((1:10)', (1:10)', 5, 'alpha')
%!error<invpred: cannot compute inverse predictions if X is constant.> ...
%! invpred (ones (10, 1), (1:10)', 5)
