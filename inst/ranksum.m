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
## @deftypefn  {statistics} @var{p} = ranksum (@var{x}, @var{y})
## @deftypefnx {statistics} @var{p} = ranksum (@var{x}, @var{y}, @var{alpha})
## @deftypefnx {statistics} @var{p} = ranksum (@var{x}, @var{y}, @var{alpha}, @var{Name}, @var{Value})
## @deftypefnx {statistics} @var{p} = ranksum (@var{x}, @var{y}, @var{Name}, @var{Value})
## @deftypefnx {statistics} [@var{p}, @var{h}] = ranksum (@var{x}, @var{y}, @dots{})
## @deftypefnx {statistics} [@var{p}, @var{h}, @var{stats}] = ranksum (@var{x}, @var{y}, @dots{})
##
## Wilcoxon rank sum test for equal medians.  This test is equivalent to a
## Mann-Whitney U-test.
##
## @code{@var{p} = ranksum (@var{x}, @var{y})} returns the p-value of a
## two-sided Wilcoxon rank sum test.  It tests the null hypothesis that two
## independent samples, in the vectors X and Y, come from continuous
## distributions with equal medians, against the alternative hypothesis that
## they are not.  @var{x} and @var{y} can have different lengths and the test
## assumes that they are independent.
##
## @code{ranksum} treats NaN in @var{x}, @var{y} as missing values.
## The two-sided p-value is computed by doubling the most significant one-sided
## value.
##
## @code{[@var{p}, @var{h}] = ranksum (@var{x}, @var{y})} also returns the
## result of the hypothesis test with @code{@var{h} = 1} indicating a rejection
## of the null hypothesis at the default alpha = 0.05 significance level, and
## @code{@var{h} = 0} indicating a failure to reject the null hypothesis at the
## same significance level.
##
## @code{[@var{p}, @var{h}, @var{stats}] = ranksum (@var{x}, @var{y})} also
## returns the structure @var{stats} with information about the test statistic.
## It contains the field @code{ranksum} with the value of the rank sum test
## statistic and if computed with the "approximate" method it also contains the
## value of the z-statistic in the field @code{zval}.
##
## @code{[@dots{}] = ranksum (@var{x}, @var{y}, @var{alpha})} or alternatively
## @code{[@dots{}] = ranksum (@var{x}, @var{y}, "alpha", @var{alpha})} returns
## the result of the hypothesis test performed at the significance level ALPHA.
##
## @code{[@dots{}] = ranksum (@var{x}, @var{y}, "method", @var{M})} defines the
## computation method of the p-value specified in @var{M}, which can be "exact",
## "approximate", or "oldexact". @var{M} must be a single string.  When "method"
## is unspecified, the default is: "exact" when
## @code{min (length (@var{x}), length (@var{y})) < 10} and
## @code{length (@var{x}) + length (@var{y}) < 10}, otherwise the "approximate"
## method is used.
##
## @itemize
## @item
## "exact" method uses full enumeration for small total sample size (< 10),
## otherwise the network algorithm is used for larger samples.
## @item
## "approximate" uses normal approximation method for computing the p-value.
## @item
## "oldexact" uses full enumeration for any sample size.  Note, that this option
## can lead to out of memory error for large samples.  Use with caution!
## @end itemize
##
## @code{[@dots{}] = ranksum (@var{x}, @var{y}, "tail", @var{tail})} defines the
## type of test, which can be "both", "right", or "left".  @var{tail} must be a
## single string.
##
## @itemize
## @item
## "both" -- "medians are not equal" (two-tailed test, default)
## @item
## "right" -- "median of X is greater than median of Y" (right-tailed test)
## @item
## "left" -- "median of X is less than median of Y" (left-tailed test)
## @end itemize
##
## Note: the rank sum statistic is based on the smaller sample of vectors
## @var{x} and @var{y}.
##
## @end deftypefn

function [p, h, stats] = ranksum(x, y, varargin)

  ## Check that x and y are vectors
  if ! isvector (x) || ! isvector (y)
     error ("X and Y must be vectors");
  endif
  ## Remove missing data and make column vectors
  x = x(! isnan (x))(:);
  y = y(! isnan (y))(:);
  if isempty (x)
    error ("Not enough data in X");
  endif
  if isempty (y)
    error ("Not enough data in Y");
  endif

  ## Check for extra input arguments
  alpha = 0.05;
  method = [];
  tail = "both";
  ## Old syntax: ranksum (x, y, alpha)
  if nargin > 2 && isnumeric (varargin{1}) && isscalar (varargin{1})
    alpha = varargin{1};
    varargin(1) = [];
    if isnan (alpha) || alpha <= 0 || alpha >= 1
      error ("Alpha does not have a valid value");
    endif
  end
  ## Check for Name:Value pairs
  arg_pairs = length (varargin);
  if ! (int16 (arg_pairs / 2) == arg_pairs / 2)
    error ("Extra arguments are not in Name:Value pairs");
  endif
  num_pair = 1;
  while (arg_pairs)
    name = varargin{num_pair};
    value = varargin{num_pair + 1};
    switch (lower (name))
      case "alpha"
        alpha = value;
        if (isnan (alpha) || alpha <= 0 || alpha >= 1 || ! isnumeric (alpha) ...
            || ! isscalar (alpha))
          error ("Alpha does not have a valid value");
        endif
      case "method"
        method = value;
        if ! any (strcmpi (method, {"exact", "approximate", "oldexact"}))
          error ("Wrong value for method option");
        endif
      case "tail"
        tail = value;
        if ! any (strcmpi (tail, {"both", "right", "left"}))
          error ("Wrong value for tail option");
        endif
    endswitch
    arg_pairs -= 2;
    num_pair += 2;
  endwhile

  ## Determine method
  nx = length (x);
  ny = length (y);
  ns = min (nx, ny);
  if isempty (method)
    if (ns < 10)  &&  ((nx + ny) < 20)
      method = "exact";
    else
      method = "approximate";
    endif
  endif

  % Determine computational technique
  switch method
    case "approximate"
      technique = "approximation";
    case "oldexact"
      technique = "exact";
    case "exact"
      if (nx + ny) < 10
        technique = "exact";
      else
        technique = "network_algorithm";
      endif
  endswitch

  % Compute the rank sum statistic based on the smaller sample
  if nx <= ny
    [ranks, tieadj] = tiedrank ([x; y]);
    x_y = true;
  else
    [ranks, tieadj] = tiedrank ([y; x]);
    x_y = false;
  endif
  srank = ranks(1:ns);
  ranksumstat = sum (srank);

  ## Calculate p-value according to selected technique
  switch technique
    case "exact"
      allpos = nchoosek (ranks, ns);
      sumranks = sum (allpos, 2);
      np = size (sumranks, 1);
      switch tail
        case "both"
          p_low = sum (sumranks <= ranksumstat) / np;
          p_high = sum (sumranks >= ranksumstat) / np;
          p = 2 * min (p_low, p_high);
          if p > 1
            p = 1;
          endif
        case "right"
          if x_y
            p = sum (sumranks >= ranksumstat) / np;
          else
            p = sum (sumranks <= ranksumstat) / np;
          endif
        case "left"
          if x_y
            p = sum (sumranks <= ranksumstat) / np;
          else
            p = sum (sumranks >= ranksumstat) / np;
          endif
      endswitch
    case "network_algorithm"
      ## Calculate contingency table
      u = unique ([x; y]);
      ct = zeros (2, length (u));
      if x_y
        ct(1,:) = histc (x,u)';
        ct(2,:) = histc (y,u)';
      else
        ct(1,:) = histc (y,u)';
        ct(2,:) = histc (x,u)';
      endif
      ## Calculate weights for wmw test
      colsum = sum (ct,1);
      tmp = cumsum (colsum);
      weights = [0 tmp(1:end - 1)] + .5 * (1 + diff ([0 tmp]));
      ## Compute p-value using network algorithm for contingency tables
      [p_net, p_val] = exact2xkCT (ct, weights, ranksumstat);
      ## Check if p = NaN
      if any (isnan (p_net)) || any (isnan (p_val))
        p = NaN;
      else
        switch tail
          case "both"
            p = 2 * p_net;
            if p > 1
              p = 1;
            endif
          case "right"
            if x_y
              p = p_val(2) + p_val(3);
            else
              p = p_val(2) + p_val(1);
            endif
          case "left"
            if x_y
              p = p_val(2) + p_val(1);
            else
              p = p_val(2) + p_val(3);
            endif
        endswitch
      endif
    case "approximation"
      wmean = ns * (nx + ny + 1) / 2;
      tiescores = 2 * tieadj / ((nx + ny) * (nx + ny - 1));
      wvar  = nx * ny * ((nx + ny + 1) - tiescores) / 12;
      wc = ranksumstat - wmean;
      ## compute z-value, including continuity correction
      switch tail
        case "both"
          z = (wc - 0.5 * sign (wc)) / sqrt (wvar);
          if ! x_y
            z = -z;
          endif
          p = 2 * normcdf (-abs(z));
        case "right"
          if x_y
            z = (wc - 0.5) / sqrt (wvar);
          else
            z = -(wc + 0.5) / sqrt (wvar);
          endif
          p = normcdf (-z);
        case "left"
          if x_y
            z = (wc + 0.5) / sqrt (wvar);
          else
            z = -(wc - 0.5) / sqrt (wvar);
          endif
          p = normcdf (z);
      endswitch
      ## For additional output argument
      if (nargout > 2)
        stats.zval = z;
      endif
  endswitch

  ## For additional output arguments
  if nargout > 1,
     h = (p <= alpha);
     if (nargout > 2)
       if x_y
         stats.ranksum = ranksumstat;
       else
         stats.ranksum = sum (ranks(ns+1:end));
       endif
     endif
  endif
endfunction

## testing against mileage data and results from Matlab
%!test
%! mileage = [33.3, 34.5, 37.4; 33.4, 34.8, 36.8; ...
%!            32.9, 33.8, 37.6; 32.6, 33.4, 36.6; ...
%!            32.5, 33.7, 37.0; 33.0, 33.9, 36.7];
%! [p,h,stats] = ranksum(mileage(:,1),mileage(:,2));
%! assert (p, 0.004329004329004329, 1e-14);
%! assert (h, true);
%! assert (stats.ranksum, 21.5);
%!test
%! year1 = [51 52 62 62 52 52 51 53 59 63 59 56 63 74 68 86 82 70 69 75 73 ...
%!          49 47 50 60 59 60 62 61 71]';
%! year2 = [54 53 64 66 57 53 54 54 62 66 59 59 67 76 75 86 82 67 74 80 75 ...
%!          54 50 53 62 62 62 72 60 67]';
%! [p,h,stats] = ranksum(year1, year2, "alpha", 0.01, "tail", "left");
%! assert (p, 0.1270832752950605, 1e-14);
%! assert (h, false);
%! assert (stats.ranksum, 837.5);
%! assert (stats.zval, -1.140287483634606, 1e-14);
%! [p,h,stats] = ranksum(year1, year2, "alpha", 0.01, "tail", "left", ...
%!                       "method", "exact");
%! assert (p, 0.127343916432862, 1e-14);
%! assert (h, false);
%! assert (stats.ranksum, 837.5);
