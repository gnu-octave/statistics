## Copyright (C) 2026 Andreas Bertsatos <abertsatos@biol.uoa.gr>
##
## This file is part of the statistics package for GNU Octave.
##
## This program is free software; you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation; either version 3 of the License, or (at your option) any later
## version.
##
## This program is distributed in the hope that it will be useful, but WITHOUT
## ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
## FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
## more details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn  {statistics} {@var{h} =} rjtest (@var{x})
## @deftypefnx {statistics} {@var{h} =} rjtest (@var{x}, @var{name}, @var{value})
## @deftypefnx {statistics} {[@var{h}, @var{p}] =} rjtest (@dots{})
## @deftypefnx {statistics} {[@var{h}, @var{p}, @var{rjstat}, @var{critval}] =} rjtest (@dots{})
##
## Ryan-Joiner hypothesis test of composite normality.
##
## @code{@var{h} = rjtest (@var{x})} performs the Ryan-Joiner test of the null
## hypothesis that the sample in the vector @var{x} comes from a normal
## distribution with unknown mean and variance, against the alternative that it
## does not come from a normal distribution.  The result @var{h} is 1 if the
## test rejects the null hypothesis at the 5% significance level, and 0
## otherwise.  @var{x} must be a vector of real values with at least four
## non-missing elements; @qcode{NaN} values are treated as missing and removed.
##
## The test statistic is the correlation between the sorted sample and the
## approximate expected normal order statistics of Blom,
## @tex
## $$ R = \frac{\sum_{i=1}^n x_{(i)} b_i}
##             {\sqrt{\sum_{i=1}^n (x_i - \bar{x})^2 \sum_{i=1}^n b_i^2}},
## \qquad b_i = \Phi^{-1}\left( \frac{i - 3/8}{n + 1/4} \right). $$
##
## @end tex
## @ifnottex
## @code{R = sum (sort (x) .* b) / sqrt (sum ((x - mean (x)).^2) * sum (b.^2))},
## where @code{b = norminv (((1:n) - 3/8) / (n + 1/4))}.
## @end ifnottex
## @math{R} lies in @math{(0,1]}, and small values of @math{R} are evidence
## against normality.  @math{R^2} is the statistic @math{W'} of the
## Shapiro-Francia test in @code{swtest}.  To test for a lognormal
## distribution, test @code{log (@var{x})}.
##
## The following @qcode{Name-Value} pairs are supported:
##
## @multitable @columnfractions 0.32 0.68
## @headitem @var{Name} @tab @var{Value}
##
## @item @qcode{'Alpha'} @tab The significance level, a double or single scalar
## in the range @math{(0,1)}.  The default is @math{0.05}.
##
## @item @qcode{'MonteCarloTolerance'} @tab A positive finite double or single
## scalar.  When given, @var{p} and @var{critval} are estimated by Monte Carlo
## simulation of normal samples of the same size, and enough samples are drawn
## for the standard error of @var{p} to fall below this value.  The name may be
## shortened to @qcode{'MonteCarloTol'}.
## @end multitable
##
## Without @qcode{'MonteCarloTolerance'} the critical value at the
## significance levels @math{0.01}, @math{0.05} and @math{0.1} is the
## approximation of Ryan and Joiner,
## @tex
## $$ c_\alpha (n) = \min \left( b_0 + \frac{b_1}{\sqrt{n}} + \frac{b_2}{n}
##                    + \frac{b_3}{n^2}, \; 1 \right), $$
##
## @end tex
## @ifnottex
## @code{min (b0 + b1 / sqrt (n) + b2 / n + b3 / n^2, 1)},
## @end ifnottex
## with the coefficients
##
## @multitable @columnfractions 0.16 0.16 0.16 0.16 0.16
## @headitem Alpha @tab b0 @tab b1 @tab b2 @tab b3
## @item 0.01 @tab 0.9963 @tab -0.0211 @tab -1.4106 @tab 3.1791
## @item 0.05 @tab 1.0063 @tab -0.1288 @tab -0.6118 @tab 1.3505
## @item 0.10 @tab 1.0071 @tab -0.1371 @tab -0.3682 @tab 0.7780
## @end multitable
##
## At other significance levels the critical value is interpolated between
## these three with a monotone piecewise cubic (@code{pchip}), extended below
## @math{0.01} along the line through the values at @math{0.01} and
## @math{0.05}, and above @math{0.1} along the line to 1 at significance level
## 1.  The p-value is interpolated linearly in @math{R} through the points
## (0, 0), (c01, 0.01), (c05, 0.05), (c10, 0.1) and (1, 1), where c01, c05 and
## c10 are the three critical values, so away from those levels it is a rough
## guide only.  The approximation is calibrated for samples of up to 100
## values, and a larger sample gives a warning; use
## @qcode{'MonteCarloTolerance'} for such samples, or for an accurate p-value.
##
## @code{[@var{h}, @var{p}] = rjtest (@dots{})} also returns the p-value
## @var{p} of the test, the probability of observing a statistic as small as
## @var{rjstat} under the null hypothesis.
##
## @code{[@var{h}, @var{p}, @var{rjstat}, @var{critval}] = rjtest (@dots{})}
## also returns the test statistic @var{rjstat} and the critical value
## @var{critval} at significance level @var{alpha}.  The null hypothesis is
## rejected when @code{@var{rjstat} < @var{critval}}.
##
## A sample whose values are all equal has no defined statistic and is refused,
## as is a sample holding an infinite value, where MATLAB returns @qcode{NaN}
## for @var{p}, @var{rjstat} and @var{critval} and 1 for @var{h}.  The option
## names are not case sensitive but must be given in full, bar
## @qcode{'MonteCarloTol'}, where MATLAB accepts any shortened name.
##
## References:
## @enumerate
## @item
## T. A. Ryan and B. L. Joiner.  Normal probability plots and tests for
## normality.  Technical report, Statistics Department, The Pennsylvania State
## University, 1976.
## @item
## G. Blom.  @emph{Statistical Estimates and Transformed Beta-Variables}.
## Wiley, New York, 1958.
## @end enumerate
##
## @seealso{swtest, adtest, jbtest, kstest, lillietest}
## @end deftypefn

function [h, p, rjstat, critval] = rjtest (x, varargin)

  ## Input validation
  if (nargin < 1)
    print_usage ();
  endif
  if (! ((isnumeric (x) || islogical (x)) && isreal (x) && isvector (x)))
    error ("rjtest: X must be a vector of real values.");
  endif

  ## Parse optional Name-Value paired arguments, 'MonteCarloTol' being the
  ## one short form taken
  k = 2 * find (strcmpi (varargin(1:2:end), 'MonteCarloTol')) - 1;
  varargin(k) = {'MonteCarloTolerance'};
  optNames = {'Alpha', 'MonteCarloTolerance'};
  dfValues = {0.05, []};
  [alpha, mctol, args] = parsePairedArguments (optNames, dfValues, ...
                                               varargin(:));
  if (! isempty (args))
    if (mod (numel (args), 2) == 0 && ischar (args{1}))
      error ("rjtest: unknown parameter name '%s'.", args{1});
    endif
    error ("rjtest: optional arguments must be Name-Value pairs.");
  endif
  if (! (isfloat (alpha) && isreal (alpha) && isscalar (alpha)
         && alpha > 0 && alpha < 1))
    error ("rjtest: 'Alpha' must be a scalar between 0 and 1.");
  endif
  domc = any (strcmpi (varargin(1:2:end), 'MonteCarloTolerance'));
  if (domc && ! (isfloat (mctol) && isreal (mctol) && isscalar (mctol)
                 && mctol > 0 && isfinite (mctol)))
    error ("rjtest: 'MonteCarloTolerance' must be a positive finite scalar.");
  endif
  mctol = double (mctol);

  ## Remove missing values and check the sample
  x = sort (double (x(! isnan (x))(:)));
  n = numel (x);
  if (any (isinf (x)))
    error ("rjtest: X must not contain infinite values.");
  endif
  if (n < 4)
    error ("rjtest: X must contain at least four non-missing values.");
  endif
  if (x(n) == x(1))
    error ("rjtest: X must not be constant.");
  endif

  ## Correlation with Blom's approximate expected normal order statistics
  b = norminv (((1:n)' - 3/8) / (n + 1/4));
  rjstat = sum (x .* b) / sqrt (sum ((x - mean (x)) .^ 2) * sum (b .^ 2));

  if (domc)
    ## Draw normal samples until the standard error of P is below MCTOL
    sim = zeros (0, 1);
    reps = 1000;
    while (true)
      sim = [sim; rjtest_simulate_(b, reps - numel (sim))];
      p = mean (sim <= rjstat);
      pq = p * (1 - p);
      if (pq / reps < mctol ^ 2)
        break;
      endif
      reps = max (ceil (pq / mctol ^ 2), reps + 1);
    endwhile
    critval = quantile (sim, double (alpha));
  else
    if (n > 100)
      warning (strcat ("rjtest: the approximation is calibrated for at", ...
                       " most 100 values; use 'MonteCarloTolerance' for a", ...
                       " more reliable p-value."));
    endif
    ## Critical values of Ryan and Joiner at ALPHA = 0.01, 0.05 and 0.1
    C = [0.9963, -0.0211, -1.4106, 3.1791;
         1.0063, -0.1288, -0.6118, 1.3505;
         1.0071, -0.1371, -0.3682, 0.7780];
    cv = min (C * [1; 1 / sqrt(n); 1 / n; 1 / n ^ 2], 1);
    ## A single ALPHA is compared with the three levels in its own precision,
    ## so single (0.05) takes the tabulated value as MATLAB does
    istab = (alpha == [0.01, 0.05, 0.1]);
    alpha = double (alpha);
    if (any (istab))
      critval = cv(istab);
    elseif (alpha < 0.01)
      critval = cv(1) + (alpha - 0.01) * (cv(2) - cv(1)) / 0.04;
    elseif (alpha <= 0.1)
      critval = interp1 ([0.01, 0.05, 0.1], cv, alpha, "pchip");
    else
      critval = cv(3) + (alpha - 0.1) * (1 - cv(3)) / 0.9;
    endif
    ## P-value linear in RJSTAT through the three critical values.  A
    ## statistic of 1 has P 1 even where the cap at 1 has met C05 or C10.
    r = min (rjstat, 1);
    if (r == 1)
      p = 1;
    elseif (r <= cv(1))
      p = 0.01 * r / cv(1);
    elseif (r <= cv(2))
      p = 0.01 + 0.04 * (r - cv(1)) / (cv(2) - cv(1));
    elseif (r <= cv(3))
      p = 0.05 + 0.05 * (r - cv(2)) / (cv(3) - cv(2));
    else
      p = 0.1 + 0.9 * (r - cv(3)) / (1 - cv(3));
    endif
  endif

  h = double (rjstat < critval);

endfunction

## Ryan-Joiner statistics of REPS standard normal samples of size numel (B),
## drawn in blocks of at most 1e6 values to bound the memory used.
function r = rjtest_simulate_ (b, reps)
  n = numel (b);
  r = zeros (reps, 1);
  block = max (floor (1e6 / n), 1);
  sb = sqrt (sum (b .^ 2));
  for k = 1:block:reps
    m = min (block, reps - k + 1);
    X = sort (randn (n, m));
    r(k:k+m-1) = (b' * X) ./ (sb * sqrt (sum ((X - mean (X)) .^ 2)));
  endfor
endfunction

%!demo
%! ## Test whether the sodium content of breakfast cereals is normal
%! load cereal
%! [h, p, rjstat, critval] = rjtest (Sodium)

%!demo
%! ## A Monte Carlo p-value for a sample of more than 100 values
%! load carbig
%! [h, p, rjstat, critval] = rjtest (MPG, 'MonteCarloTolerance', 0.001)

## Test output
## Values checked against MATLAB R2026b Update 1
%!test
%! b = norminv (((1:20) - 3/8) / (20 + 1/4));
%! [h, p, rjstat, critval] = rjtest (exp (0.6 * b));
%! assert_equal (h, 1);
%! assert_equal (p, 0.046583958270526, 1e-14);
%! assert_equal (rjstat, 0.948467843754262, 1e-14);
%! assert_equal (critval, 0.950285694449803, 1e-14);
%!test
%! b = norminv (((1:20) - 3/8) / (20 + 1/4));
%! [h, p, rjstat, critval] = rjtest (exp (0.6 * b), 'Alpha', 0.01);
%! assert_equal (h, 0);
%! assert_equal (p, 0.046583958270526, 1e-14);
%! assert_equal (critval, 0.928999646567475, 1e-14);
%!test
%! b = norminv (((1:20) - 3/8) / (20 + 1/4));
%! [h, ~, ~, critval] = rjtest (exp (0.6 * b), 'Alpha', 0.2);
%! assert_equal (h, 1);
%! assert_equal (critval, 0.964425340469758, 1e-14);
%!test
%! b = norminv (((1:20) - 3/8) / (20 + 1/4));
%! [h, p, rjstat] = rjtest (b);
%! assert_equal (h, 0);
%! assert_equal (p, 1);
%! assert_equal (rjstat, 1);
%!test
%! b = norminv (((1:20) - 3/8) / (20 + 1/4));
%! [h, p, rjstat] = rjtest (exp (0.1 * b));
%! assert_equal (h, 0);
%! assert_equal (p, 0.964707970026192, 1e-13);
%! assert_equal (rjstat, 0.998430622561938, 1e-14);
%!test
%! b = norminv (((1:20) - 3/8) / (20 + 1/4));
%! [h, p, rjstat] = rjtest (exp (0.55 * b));
%! assert_equal (h, 0);
%! assert_equal (p, 0.079811798723863, 1e-14);
%! assert_equal (rjstat, 0.956064898599311, 1e-14);
%!test
%! b = norminv (((1:20) - 3/8) / (20 + 1/4));
%! [h, p, rjstat] = rjtest (exp (0.7 * b));
%! assert_equal (h, 1);
%! assert_equal (p, 0.015775669183456, 1e-14);
%! assert_equal (rjstat, 0.932073175837263, 1e-14);
%!test
%! b = norminv (((1:20) - 3/8) / (20 + 1/4));
%! [h, p, rjstat] = rjtest (exp (3 * b));
%! assert_equal (h, 1);
%! assert_equal (p, 0.006366180860427, 1e-14);
%! assert_equal (rjstat, 0.591417976932142, 1e-14);
%!test
%! b = norminv (((1:100) - 3/8) / (100 + 1/4));
%! [h, p, rjstat, critval] = rjtest (exp (0.2 * b));
%! assert_equal (h, 0);
%! assert_equal (p, 0.251186660905392, 1e-12);
%! assert_equal (rjstat, 0.991501634213133, 1e-14);
%! assert_equal (critval, 0.98743705, 1e-14);
%!test
%! b = norminv (((1:5) - 3/8) / (5 + 1/4));
%! [h, p, rjstat, critval] = rjtest (exp (2 * b));
%! assert_equal (h, 1);
%! assert_equal (p, 0.026740104695873, 1e-14);
%! assert_equal (rjstat, 0.852184703528202, 1e-14);
%! assert_equal (critval, 0.880358888899605, 1e-14);
%!test
%! [h, p, rjstat, critval] = rjtest ([1 2 3 5]);
%! assert_equal (h, 0);
%! assert_equal (p, 0.871931434395105, 1e-14);
%! assert_equal (rjstat, 0.985076454646874, 1e-14);
%! assert_equal (critval, 0.87335625, 1e-14);
%!test
%! [h, p, rjstat, critval] = rjtest (int8 ([1 2 3 5 8 13]));
%! assert_equal (h, 0);
%! assert_equal (p, 0.520710875770748, 1e-14);
%! assert_equal (rjstat, 0.952802591324644, 1e-14);
%! assert_equal (critval, 0.889264842410477, 1e-14);
%!test
%! [h, p, rjstat] = rjtest (logical ([1 0 1 1 0]));
%! assert_equal (h, 1);
%! assert_equal (p, 0.021235058557769, 1e-14);
%! assert_equal (rjstat, 0.845516565587550, 1e-14);
%!test
%! b = norminv (((1:20) - 3/8) / (20 + 1/4));
%! [~, p, rjstat] = rjtest (single (exp (0.6 * b)));
%! assert_equal (class (rjstat), 'double');
%! assert_equal (p, 0.046583962851523, 1e-14);
%! assert_equal (rjstat, 0.948467846192045, 1e-14);
%!test
%! b = norminv (((1:20) - 3/8) / (20 + 1/4));
%! [h, p, rjstat] = rjtest ([NaN, exp(0.6 * b), NaN]');
%! assert_equal (h, 1);
%! assert_equal (p, 0.046583958270526, 1e-14);
%! assert_equal (rjstat, 0.948467843754262, 1e-14);
%!test
%! load cereal
%! [h, p, rjstat, critval] = rjtest (Sodium);
%! assert_equal (h, 1);
%! assert_equal (p, 0.009908718822304, 1e-14);
%! assert_equal (rjstat, 0.967202079663432, 1e-14);
%! assert_equal (critval, 0.983904201921336, 1e-14);
%!test
%! load cereal
%! [~, ~, ~, critval] = rjtest (Sodium, 'Alpha', 0.01);
%! assert_equal (critval, 0.976112146291078, 1e-14);
%!warning <rjtest: the approximation is calibrated for at most 100 values> ...
%! load carbig
%! [h, p, rjstat, critval] = rjtest (MPG);
%! assert_equal (h, 1);
%! assert_equal (p, 0.009928312778520, 1e-14);
%! assert_equal (rjstat, 0.984608848799090, 1e-14);
%! assert_equal (critval, 0.998315179111227, 1e-14);

%!warning <rjtest: the approximation is calibrated for at most 100 values> ...
%! b = norminv (((1:500) - 3/8) / (500 + 1/4));
%! [h, p, rjstat, critval] = rjtest (b);
%! assert_equal (h, 0);
%! assert_equal (p, 1);
%! assert_equal (rjstat, 1);
%! assert_equal (critval, 0.999321690889961, 1e-14);
%!warning <rjtest: the approximation is calibrated for at most 100 values> ...
%! b = norminv (((1:500) - 3/8) / (500 + 1/4));
%! [h, p, rjstat] = rjtest (exp (0.02 * b));
%! assert_equal (h, 0);
%! assert_equal (p, 0.092980498877550, 1e-11);
%! assert_equal (rjstat, 0.999904772168814, 1e-13);
%!warning <rjtest: the approximation is calibrated for at most 100 values> ...
%! b = norminv (((1:500) - 3/8) / (500 + 1/4));
%! [h, p, rjstat] = rjtest (exp (0.08 * b));
%! assert_equal (h, 1);
%! assert_equal (p, 0.045015541119476, 1e-12);
%! assert_equal (rjstat, 0.998477598301906, 1e-13);
%!warning <rjtest: the approximation is calibrated for at most 100 values> ...
%! b = norminv (((1:1000) - 3/8) / (1000 + 1/4));
%! [h, p, rjstat, critval] = rjtest (exp (0.02 * b));
%! assert_equal (h, 1);
%! assert_equal (p, 0.049327273854406, 1e-12);
%! assert_equal (rjstat, 0.999902880855905, 1e-13);
%! assert_equal (critval, 1);
%!warning <rjtest: the approximation is calibrated for at most 100 values> ...
%! b = norminv (((1:1000) - 3/8) / (1000 + 1/4));
%! [h, p] = rjtest (b);
%! assert_equal (h, 0);
%! assert_equal (p, 1);

## Test the critical values
## Values checked against MATLAB R2026b Update 1
%!test
%! [~, ~, ~, critval] = rjtest ([1 2 3 5], 'Alpha', 0.01);
%! assert_equal (critval, 0.83179375, 1e-15);
%!test
%! [~, ~, ~, critval] = rjtest ([1 2 3 5], 'Alpha', 0.1);
%! assert_equal (critval, 0.895125, 1e-15);
%!test
%! b = norminv (((1:7) - 3/8) / (7 + 1/4));
%! [~, ~, ~, critval] = rjtest (b, 'Alpha', 0.05);
%! assert_equal (critval, 0.897779400366207, 1e-14);
%!test
%! b = norminv (((1:50) - 3/8) / (50 + 1/4));
%! [~, ~, ~, critval] = rjtest (b, 'Alpha', 0.1);
%! assert_equal (critval, 0.980658332059865, 1e-14);
%!test
%! b = norminv (((1:20) - 3/8) / (20 + 1/4));
%! [~, ~, ~, critval] = rjtest (b, 'Alpha', 0.001);
%! assert_equal (critval, 0.924210285793952, 1e-14);
%!test
%! b = norminv (((1:20) - 3/8) / (20 + 1/4));
%! [~, ~, ~, critval] = rjtest (b, 'Alpha', 0.03);
%! assert_equal (critval, 0.941609306183096, 1e-14);
%!test
%! b = norminv (((1:20) - 3/8) / (20 + 1/4));
%! [~, ~, ~, critval] = rjtest (b, 'Alpha', 0.075);
%! assert_equal (critval, 0.956902493214627, 1e-14);
%!test
%! b = norminv (((1:20) - 3/8) / (20 + 1/4));
%! [~, ~, ~, critval] = rjtest (b, 'Alpha', 0.9);
%! assert_equal (critval, 0.995553167558720, 1e-14);
%!test
%! b = norminv (((1:100) - 3/8) / (100 + 1/4));
%! [~, ~, ~, critval] = rjtest (b, 'Alpha', 0.2);
%! assert_equal (critval, 0.990920711111111, 1e-14);
%!test
%! b = norminv (((1:20) - 3/8) / (20 + 1/4));
%! [~, ~, ~, critval] = rjtest (b, 'Alpha', single (0.05));
%! assert_equal (critval, 0.950285694449803, 1e-15);

%!warning <rjtest: the approximation is calibrated for at most 100 values> ...
%! b = norminv (((1:500) - 3/8) / (500 + 1/4));
%! [~, ~, ~, critval] = rjtest (b, 'Alpha', 0.075);
%! assert_equal (critval, 0.999822959048033, 1e-14);
%!warning <rjtest: the approximation is calibrated for at most 100 values> ...
%! b = norminv (((1:1000) - 3/8) / (1000 + 1/4));
%! [~, ~, ~, critval] = rjtest (b, 'Alpha', 0.03);
%! assert_equal (critval, 0.998155316469655, 1e-14);
%!warning <rjtest: the approximation is calibrated for at most 100 values> ...
%! b = norminv (((1:1000) - 3/8) / (1000 + 1/4));
%! [~, ~, ~, critval] = rjtest (b, 'Alpha', 0.005);
%! assert_equal (critval, 0.993503505827917, 1e-14);
%!warning <rjtest: the approximation is calibrated for at most 100 values> ...
%! b = norminv (((1:1000) - 3/8) / (1000 + 1/4));
%! [~, ~, ~, critval] = rjtest (b, 'Alpha', 0.2);
%! assert_equal (critval, 1);

## Test the Monte Carlo p-value
%!test
%! randn ('state', 1);
%! b = norminv (((1:20) - 3/8) / (20 + 1/4));
%! [~, p] = rjtest (exp (3 * b), 'MonteCarloTol', 1e-4);
%! assert_equal (p, 0);
%!test
%! randn ('state', 1);
%! b = norminv (((1:20) - 3/8) / (20 + 1/4));
%! [h, p, rjstat, critval] = rjtest (exp (0.6 * b), 'MonteCarloTolerance', ...
%!                                   0.001);
%! assert_equal (h, 1);
%! assert_equal (p, 0.042123485285632, 0.005);
%! assert_equal (rjstat, 0.948467843754262, 1e-14);
%! assert_equal (critval, 0.950586583658806, 0.002);
%!test
%! randn ('state', 1);
%! b = norminv (((1:20) - 3/8) / (20 + 1/4));
%! [h, p] = rjtest (exp (0.3 * b), 'MonteCarloTolerance', 0.001);
%! assert_equal (h, 0);
%! assert_equal (p, 0.718153708746535, 0.005);
%!test
%! randn ('state', 1);
%! b = norminv (((1:20) - 3/8) / (20 + 1/4));
%! [h, p] = rjtest (exp (3 * b), 'MonteCarloTolerance', 1e-4);
%! assert_equal (h, 1);
%! assert_equal (p, 0);
%!test
%! randn ('state', 1);
%! b = norminv (((1:20) - 3/8) / (20 + 1/4));
%! [~, p] = rjtest (b, 'MonteCarloTolerance', 1);
%! assert_equal (p, 1);
%!test
%! randn ('state', 1);
%! b = norminv (((1:20) - 3/8) / (20 + 1/4));
%! [~, ~, ~, critval] = rjtest (exp (0.6 * b), 'MonteCarloTolerance', 0.01, ...
%!                              'alpha', 0.3);
%! assert_equal (critval, 0.974724689853778, 0.003);
%!test
%! randn ('state', 1);
%! b = norminv (((1:1000) - 3/8) / (1000 + 1/4));
%! [h, p, ~, critval] = rjtest (exp (0.3 * b), 'MonteCarloTolerance', 0.01);
%! assert_equal (h, 1);
%! assert_equal (p, 0);
%! assert_equal (critval, 0.998460671929262, 2e-4);

## Test input validation
%!error <Invalid call to rjtest> rjtest ()
%!error <rjtest: X must be a vector of real values.> rjtest ([])
%!error <rjtest: X must be a vector of real values.> rjtest (magic (4))
%!error <rjtest: X must be a vector of real values.> rjtest ([1 2 3 4 5] + 1i)
%!error <rjtest: X must be a vector of real values.> rjtest ('abcde')
%!error <rjtest: X must be a vector of real values.> rjtest ({1, 2, 3, 4})
%!error <rjtest: unknown parameter name 'Foo'.> rjtest (1:10, 'Foo', 0.05)
%!error <rjtest: unknown parameter name 'Alp'.> rjtest (1:10, 'Alp', 0.05)
%!error <rjtest: unknown parameter name 'MonteCarlo'.> ...
%! rjtest (1:10, 'MonteCarlo', 0.01)
%!error <rjtest: 'MonteCarloTolerance' must be a positive finite scalar.> ...
%! rjtest (1:10, 'montecarlotol', 0)
%!error <rjtest: optional arguments must be Name-Value pairs.> ...
%! rjtest (1:10, 'Alpha')
%!error <rjtest: 'Alpha' must be a scalar between 0 and 1.> ...
%! rjtest (1:10, 'Alpha', 0)
%!error <rjtest: 'Alpha' must be a scalar between 0 and 1.> ...
%! rjtest (1:10, 'Alpha', 1)
%!error <rjtest: 'Alpha' must be a scalar between 0 and 1.> ...
%! rjtest (1:10, 'Alpha', NaN)
%!error <rjtest: 'Alpha' must be a scalar between 0 and 1.> ...
%! rjtest (1:10, 'Alpha', [0.05 0.1])
%!error <rjtest: 'Alpha' must be a scalar between 0 and 1.> ...
%! rjtest (1:10, 'Alpha', int8 (1))
%!error <rjtest: 'MonteCarloTolerance' must be a positive finite scalar.> ...
%! rjtest (1:10, 'MonteCarloTolerance', 0)
%!error <rjtest: 'MonteCarloTolerance' must be a positive finite scalar.> ...
%! rjtest (1:10, 'MonteCarloTolerance', -1)
%!error <rjtest: 'MonteCarloTolerance' must be a positive finite scalar.> ...
%! rjtest (1:10, 'MonteCarloTolerance', NaN)
%!error <rjtest: 'MonteCarloTolerance' must be a positive finite scalar.> ...
%! rjtest (1:10, 'MonteCarloTolerance', Inf)
%!error <rjtest: 'MonteCarloTolerance' must be a positive finite scalar.> ...
%! rjtest (1:10, 'MonteCarloTolerance', [])
%!error <rjtest: X must not contain infinite values.> rjtest ([1 2 3 Inf 5 6])
%!error <rjtest: X must contain at least four non-missing values.> ...
%! rjtest ([1 2 3 NaN])
%!error <rjtest: X must not be constant.> rjtest (ones (1, 10))
