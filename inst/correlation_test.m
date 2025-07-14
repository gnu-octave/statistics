## Copyright (C) 1995-2017 Kurt Hornik
## Copyright (C) 2023 Andreas Bertsatos <abertsatos@biol.uoa.gr>
##
## This file is part of the statistics package for GNU Octave.
##
## This program is free software: you can redistribute it and/or
## modify it under the terms of the GNU General Public License as
## published by the Free Software Foundation, either version 3 of the
## License, or (at your option) any later version.
##
## This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; see the file COPYING.  If not, see
## <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn  {statistics} {@var{h} =} correlation_test (@var{x}, @var{y})
## @deftypefnx {statistics} {[@var{h}, @var{pval}] =} correlation_test (@var{y}, @var{x})
## @deftypefnx {statistics} {[@var{h}, @var{pval}, @var{stats}] =} correlation_test (@var{y}, @var{x})
## @deftypefnx {statistics} {[@dots{}] =} correlation_test (@var{y}, @var{x}, @var{Name}, @var{Value})
##
## Perform a correlation coefficient test to determine whether two samples
## @var{x} and @var{y} come from uncorrelated populations.
##
## @code{@var{h} = correlation_test (@var{y}, @var{x})} tests the null
## hypothesis that the two samples @var{x} and @var{y} come from uncorrelated
## populations.  The result is @var{h} = 0 if the null hypothesis cannot be
### rejected at the 5% significance level, or @var{h} = 1 if the null hypothesis
## can be rejected at the 5% level.  @var{y} and @var{x} must be vectors of
## equal length with finite real numbers.
##
## The p-value of the test is returned in @var{pval}.  @var{stats} is a
## structure with the following fields:
## @multitable @columnfractions 0.05 0.2 0.05 0.70
## @headitem @tab Field @tab @tab Value
## @item @tab @qcode{method} @tab @tab the type of correlation coefficient used
## for the test
## @item @tab @qcode{df} @tab @tab the degrees of freedom (where applicable)
## @item @tab @qcode{corrcoef} @tab @tab the correlation coefficient
## @item @tab @qcode{stat} @tab @tab the test's statistic
## @item @tab @qcode{dist} @tab @tab the respective distribution for the test
## @item @tab @qcode{alt} @tab @tab the alternative hypothesis for the test
## @end multitable
##
##
## @code{[@dots{}] = correlation_test (@dots{}, @var{name}, @var{value})}
## specifies one or more of the following name/value pairs:
##
## @multitable @columnfractions 0.05 0.2 0.75
## @headitem @tab Name @tab Value
## @item @tab @qcode{"alpha"} @tab the significance level. Default is 0.05.
##
## @item @tab @qcode{"tail"} @tab a string specifying the alternative hypothesis
## @end multitable
## @multitable @columnfractions 0.1 0.25 0.65
## @item @tab @qcode{"both"} @tab @math{corrcoef} is not 0 (two-tailed, default)
## @item @tab @qcode{"left"} @tab @math{corrcoef} is less than 0 (left-tailed)
## @item @tab @qcode{"right"} @tab @math{corrcoef} is greater than 0
## (right-tailed)
## @end multitable
##
## @multitable @columnfractions 0.05 0.2 0.75
## @item @tab @qcode{"method"} @tab a string specifying the correlation
## coefficient used for the test
## @end multitable
## @multitable @columnfractions 0.1 0.25 0.65
## @item @tab @qcode{"pearson"} @tab Pearson's product moment correlation
## (Default)
## @item @tab @qcode{"kendall"} @tab Kendall's rank correlation tau
## @item @tab @qcode{"spearman"} @tab Spearman's rank correlation rho
## @end multitable
##
## @seealso{regression_ftest, regression_ttest}
## @end deftypefn

function [h, pval, stats] = correlation_test (x, y, varargin)

  if (nargin < 2)
    print_usage ();
  endif

  if (! isvector (x) || ! isvector (y) || length (x) != length (y))
    error ("correlation_test: X and Y must be vectors of equal length.");
  endif

  ## Force to column vectors
  x = x(:);
  y = y(:);

  ## Check for finite real numbers in X and Y
  if (! all (isfinite (x)) || ! isreal (x))
    error ("correlation_test: X must contain finite real numbers.");
  endif
  if (! all (isfinite (y(:))) || ! isreal (y))
    error ("correlation_test: Y must contain finite real numbers.");
  endif

  ## Set default arguments
  alpha = 0.05;
  tail = "both";
  method = "pearson";

  ## Check additional options
  i = 1;
  while (i <= length (varargin))
    switch lower (varargin{i})
      case "alpha"
        i = i + 1;
        alpha = varargin{i};
        ## Check for valid alpha
        if (! isscalar (alpha) || ! isnumeric (alpha) || ...
                    alpha <= 0 || alpha >= 1)
          error ("correlation_test: invalid value for alpha.");
        endif
      case "tail"
        i = i + 1;
        tail = varargin{i};
        if (! any (strcmpi (tail, {"both", "left", "right"})))
          error ("correlation_test: invalid value for tail.");
        endif
      case "method"
        i = i + 1;
        method = varargin{i};
        if (! any (strcmpi (method, {"pearson", "kendall", "spearman"})))
          error ("correlation_test: invalid value for method.");
        endif
      otherwise
        error ("correlation_test: invalid Name argument.");
    endswitch
    i = i + 1;
  endwhile

  n = length (x);

  if (strcmpi (method, "pearson"))
    r = corr (x, y);
    stats.method = "Pearson's product moment correlation";
    stats.df = n - 2;
    stats.corrcoef = r;
    stats.stat = sqrt (stats.df) .* r / sqrt (1 - r.^2);
    stats.dist = "Student's t";
    cdf = tcdf (stats.stat, stats.df);
  elseif (strcmpi (method, "kendall"))
    tau = kendall (x, y);
    stats.method = "Kendall's rank correlation tau";
    stats.df = [];
    stats.corrcoef = tau;
    stats.stat = tau / sqrt ((2 * (2*n+5)) / (9*n*(n-1)));
    stats.dist = "standard normal";
    cdf = stdnormal_cdf (stats.stat);
  else  # spearman
    rho = spearman (x, y);
    stats.method = "Spearman's rank correlation rho";
    stats.df = [];
    stats.corrcoef = rho;
    stats.stat = sqrt (n-1) * (rho - 6/(n^3-n));
    stats.dist = "standard normal";
    cdf = stdnormal_cdf (stats.stat);
  endif

  ## Based on the "tail" argument determine the P-value
  switch lower (tail)
    case "both"
      pval = 2 * min (cdf, 1 - cdf);
    case "right"
      pval = 1 - cdf;
    case "left"
      pval = cdf;
  endswitch

  stats.alt = tail;

  ## Determine the test outcome
  h = double (pval < alpha);

endfunction

## Test input validation
%!error<Invalid call to correlation_test.  Correct usage> correlation_test ();
%!error<Invalid call to correlation_test.  Correct usage> correlation_test (1);
%!error<correlation_test: X must contain finite real numbers.> ...
%! correlation_test ([1 2 NaN]', [2 3 4]');
%!error<correlation_test: X must contain finite real numbers.> ...
%! correlation_test ([1 2 Inf]', [2 3 4]');
%!error<correlation_test: X must contain finite real numbers.> ...
%! correlation_test ([1 2 3+i]', [2 3 4]');
%!error<correlation_test: Y must contain finite real numbers.> ...
%! correlation_test ([1 2 3]', [2 3 NaN]');
%!error<correlation_test: Y must contain finite real numbers.> ...
%! correlation_test ([1 2 3]', [2 3 Inf]');
%!error<correlation_test: Y must contain finite real numbers.> ...
%! correlation_test ([1 2 3]', [3 4 3+i]');
%!error<correlation_test: X and Y must be vectors of equal length.> ...
%! correlation_test ([1 2 3]', [3 4 4 5]');
%!error<correlation_test: invalid value for alpha.> ...
%! correlation_test ([1 2 3]', [2 3 4]', "alpha", 0);
%!error<correlation_test: invalid value for alpha.> ...
%! correlation_test ([1 2 3]', [2 3 4]', "alpha", 1.2);
%!error<correlation_test: invalid value for alpha.> ...
%! correlation_test ([1 2 3]', [2 3 4]', "alpha", [.02 .1]);
%!error<correlation_test: invalid value for alpha.> ...
%! correlation_test ([1 2 3]', [2 3 4]', "alpha", "a");
%!error<correlation_test: invalid Name argument.> ...
%! correlation_test ([1 2 3]', [2 3 4]', "some", 0.05);
%!error<correlation_test: invalid value for tail.>  ...
%! correlation_test ([1 2 3]', [2 3 4]', "tail", "val");
%!error<correlation_test: invalid value for tail.>  ...
%! correlation_test ([1 2 3]', [2 3 4]', "alpha", 0.01, "tail", "val");
%!error<correlation_test: invalid value for method.>  ...
%! correlation_test ([1 2 3]', [2 3 4]', "method", 0.01);
%!error<correlation_test: invalid value for method.>  ...
%! correlation_test ([1 2 3]', [2 3 4]', "method", "some");

%!test
%! x = [6 7 7 9 10 12 13 14 15 17];
%! y = [19 22 27 25 30 28 30 29 25 32];
%! [h, pval, stats] = correlation_test (x, y);
%! assert (stats.corrcoef, corr (x', y'), 1e-14);
%! assert (pval, 0.0223, 1e-4);
%!test
%! x = [6 7 7 9 10 12 13 14 15 17]';
%! y = [19 22 27 25 30 28 30 29 25 32]';
%! [h, pval, stats] = correlation_test (x, y);
%! assert (stats.corrcoef, corr (x, y), 1e-14);
%! assert (pval, 0.0223, 1e-4);
