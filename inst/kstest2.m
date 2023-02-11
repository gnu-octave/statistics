## Copyright (C) 2022 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
## details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn  {statistics} {@var{h} =} kstest2 (@var{x1}, @var{x2})
## @deftypefnx {statistics} {@var{h} =} kstest2 (@var{x1}, @var{x2}, @var{name}, @var{value})
## @deftypefnx {statistics} {[@var{h}, @var{p}] =} kstest2 (@dots{})
## @deftypefnx {statistics} {[@var{h}, @var{p}, @var{ks2stat}] =} kstest2 (@dots{})
##
## Two-sample Kolmogorov-Smirnov goodness-of-fit hypothesis test.
##
## @code{@var{h} = kstest2 (@var{x1}, @var{x2})} returns a test decision for the
## null hypothesis that the data in vectors @var{x1} and @var{x2} are from the
## same continuous distribution, using the two-sample Kolmogorov-Smirnov test.
## The alternative hypothesis is that @var{x1} and @var{x2} are from different
## continuous distributions.  The result @var{h} is 1 if the test rejects the
## null hypothesis at the 5% significance level, and 0 otherwise.
##
## @code{@var{h} = kstest2 (@var{x1}, @var{x2}, @var{name}, @var{value})}
## returns a test decision for a two-sample Kolmogorov-Smirnov test with
## additional options specified by one or more name-value pair arguments as
## shown below.
##
## @multitable @columnfractions 0.20 0.8
## @item "alpha" @tab A value @var{alpha} between 0 and 1 specifying the
## significance level.  Default is 0.05 for 5% significance.
##
## @item "tail" @tab A string indicating the type of test:
## @end multitable
##
## @multitable @columnfractions 0.03 0.2 0.77
## @item @tab "unequal" @tab "F(X1) not equal to F(X2)" (two-sided) [Default]
##
## @item @tab "larger" @tab "F(X1) > F(X2)" (one-sided)
##
## @item @tab "smaller" @tab "F(X1) < F(X2)" (one-sided)
## @end multitable
##
## The two-sided test uses the maximum absolute difference between the cdfs of
## the distributions of the two data vectors.  The test statistic is
## @code{D* = max(|F1(x) - F2(x)|)}, where F1(x) is the proportion of @var{x1}
## values less or equal to x and F2(x) is the proportion of @var{x2} values less
## than or equal to x.  The one-sided test uses the actual value of the
## difference between the cdfs of the distributions of the two data vectors
## rather than the absolute value. The test statistic is
## @code{D* = max(F1(x) - F2(x))} or @code{D* = max(F2(x) - F1(x))} for
## @code{tail} = "larger" or "smaller", respectively.
##
## @code{[@var{h}, @var{p}] = kstest2 (@dots{})} also returns the
## asymptotic p-value @var{p}.
##
## @code{[@var{h}, @var{p}, @var{ks2stat}] = kstest2 (@dots{})} also returns
## the Kolmogorov-Smirnov test statistic @var{ks2stat} defined above for the
## test type indicated by @code{tail}.
##
## @seealso{kstest, cdfplot}
## @end deftypefn

function [H, pValue, ks2stat] = kstest2 (x1, x2, varargin)
  ## Check input parameters
  if nargin < 2
    error ("kstest2: Too few inputs.");
  endif
  if ! isvector (x1) || ! isreal (x1) || ! isvector (x2) || ! isreal (x2)
    error ("kstest2: X1 and X2 must be vectors of real numbers.");
  endif
  ## Add defaults
  alpha = 0.05;
  tail = "unequal";
  ## Parse extra parameters
  if nargin > 2 && mod (numel (varargin), 2) == 0
    [~, prop] = parseparams (varargin);
    while (!isempty (prop))
      switch (lower (prop{1}))
        case "alpha"
          alpha = prop{2};
        case "tail"
          tail = prop{2};
        otherwise
          error ("kstest2: Unknown option %s", prop{1});
      endswitch
      prop = prop(3:end);
    endwhile
  elseif nargin > 2
    error ("kstest2: optional parameters must be in name/value pairs.");
  endif
  ## Check for valid alpha and tail parameters
  if (! isnumeric (alpha) || isnan (alpha) || ! isscalar (alpha) ...
                          || alpha <= 0 || alpha >= 1)
    error ("kstest2: alpha must be a numeric scalar in the range (0,1).");
  endif
  if ! isa (tail, 'char')
    error ("kstest2: tail argument must be a string");
  elseif sum (strcmpi (tail, {"unequal", "larger", "smaller"})) < 1
    error ("kstest2: tail value must be either 'both', right' or 'left'.");
  endif
  ## Make x1 and x2 column vectors
  x1 = x1(:);
  x2 = x2(:);
  ## Remove missing values (NaN)
  x1(isnan (x1)) = [];
  x2(isnan (x2)) = [];
  ## Check for remaining data in both vectors
  if isempty (x1)
    error ("kstest2: Not enough data in X1");
  elseif isempty (x2)
    error ("kstest2: Not enough data in X2");
  endif
  ## Calculate F1(x) and F2(x)
  binEdges    =  [-inf; sort([x1;x2]); inf];
  binCounts1  =  histc (x1 , binEdges, 1);
  binCounts2  =  histc (x2 , binEdges, 1);
  sumCounts1  =  cumsum (binCounts1) ./ sum (binCounts1);
  sumCounts2  =  cumsum (binCounts2) ./ sum (binCounts2);
  sampleCDF1  =  sumCounts1(1:end - 1);
  sampleCDF2  =  sumCounts2(1:end - 1);
  ## Calculate the suitable KS statistic according to tail
  switch tail
    case "unequal"    # 2-sided test: T = max|F1(x) - F2(x)|.
      deltaCDF  =  abs (sampleCDF1 - sampleCDF2);
    case "smaller"    # 1-sided test: T = max[F2(x) - F1(x)].
      deltaCDF  =  sampleCDF2 - sampleCDF1;
    case "larger"     # 1-sided test: T = max[F1(x) - F2(x)].
      deltaCDF  =  sampleCDF1 - sampleCDF2;
  endswitch
  ks2stat = max (deltaCDF);
  ## Compute the asymptotic P-value approximation
  n_x1 = length(x1);
  n_x2 = length(x2);
  n =  n_x1 * n_x2 /(n_x1 + n_x2);
  lambda = max ((sqrt (n) + 0.12 + 0.11 / sqrt (n)) * ks2stat, 0);
  if strcmpi (tail, "unequal")    # 2-sided test
    v = [1:101];
    pValue = 2 * sum ((-1) .^ (v-1) .* exp (-2 * lambda * lambda * v .^ 2));
    pValue = min (max (pValue, 0), 1);
  else                            # 1-sided test
    pValue  =  exp(-2 * lambda * lambda);
  endif
  ## Return hypothesis test
  H = (alpha >= pValue);
endfunction

## Test input
%!error kstest2 ([1,2,3,4,5,5])
%!error kstest2 (ones(2,4), [1,2,3,4,5,5])
%!error kstest2 ([2,3,5,7,3+3i], [1,2,3,4,5,5])
%!error kstest2 ([2,3,4,5,6],[3;5;7;8;7;6;5],"tail")
%!error kstest2 ([2,3,4,5,6],[3;5;7;8;7;6;5],"tail", "whatever")
%!error kstest2 ([2,3,4,5,6],[3;5;7;8;7;6;5],"badoption", 0.51)
%!error kstest2 ([2,3,4,5,6],[3;5;7;8;7;6;5],"tail", 0)
%!error kstest2 ([2,3,4,5,6],[3;5;7;8;7;6;5],"alpha", 0)
%!error kstest2 ([2,3,4,5,6],[3;5;7;8;7;6;5],"alpha", NaN)
%!error kstest2 ([NaN,NaN,NaN,NaN,NaN],[3;5;7;8;7;6;5],"tail", "unequal")

## Test results
%!test
%! load examgrades
%! [h, p] = kstest2 (grades(:,1), grades(:,2));
%! assert (h, false);
%! assert (p, 0.1222791870137312, 1e-14);
%!test
%! load examgrades
%! [h, p] = kstest2 (grades(:,1), grades(:,2), "tail", "larger");
%! assert (h, false);
%! assert (p, 0.1844421391011258, 1e-14);
%!test
%! load examgrades
%! [h, p] = kstest2 (grades(:,1), grades(:,2), "tail", "smaller");
%! assert (h, false);
%! assert (p, 0.06115357930171663, 1e-14);
%!test
%! load examgrades
%! [h, p] = kstest2 (grades(:,1), grades(:,2), "tail", "smaller", "alpha", 0.1);
%! assert (h, true);
%! assert (p, 0.06115357930171663, 1e-14);
