## Copyright (C) 2014 Tony Richardson
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
## @deftypefn  {Function File} {[@var{pval}, @var{h}, @var{stats}] =} signtest (@var{x})
## @deftypefnx {Function File} {[@var{pval}, @var{h}, @var{stats}] =} signtest (@var{x}, @var{m})
## @deftypefnx {Function File} {[@var{pval}, @var{h}, @var{stats}] =} signtest (@var{x}, @var{y})
## @deftypefnx {Function File} {[@var{pval}, @var{h}, @var{stats}] =} signtest (@var{x}, @var{y}, @var{Name}, @var{Value})
##
## Test for median.
##
## Perform a signtest of the null hypothesis that @var{x} is from a distribution
## that has a zero median.  @var{x} must be a vector.
##
## If the second argument @var{m} is a scalar, the null hypothesis is that
## X has median m.
##
## If the second argument @var{y} is a vector, the null hypothesis is that
## the distribution of @code{@var{x} - @var{y}} has zero median.
##
## The argument @qcode{"alpha"} can be used to specify the significance level
## of the test (the default value is 0.05).  The string
## argument @qcode{"tail"}, can be used to select the desired alternative
## hypotheses.  If @qcode{"tail"} is @qcode{"both"} (default) the null is
## tested against the two-sided alternative @code{median (@var{x}) != @var{m}}.
## If @qcode{"tail"} is @qcode{"right"} the one-sided
## alternative @code{median (@var{x}) > @var{m}} is considered.
## Similarly for @qcode{"left"}, the one-sided alternative @code{median
## (@var{x}) < @var{m}} is considered.
##
## When @qcode{"method"} is @qcode{"exact"} the p-value is computed using an
## exact method.  When @qcode{"method"} is @qcode{"approximate"} a normal
## approximation is used for the test statistic. When @qcode{"method"} is not
## defined as an optional input argument, then for @code{length (@var{x}) < 100}
## the @qcode{"exact"} method is computed, otherwise the @qcode{"approximate"}
## method is used.
##
## The p-value of the test is returned in @var{pval}. If @var{h} is 0 the
## null hypothesis is accepted, if it is 1 the null hypothesis is rejected.
## @var{stats} is a structure containing the value of the test statistic
## (@var{sign}) and the value of the z statistic (@var{zval}) (only computed
## when the 'method' is 'approximate'.
##
## @end deftypefn

function [p, h, stats] = signtest (x, my, varargin)

  ## Check X being a vector
  if ! isvector (x)
    error ("signtest: X must be a vector.");
  endif
  ## Add defaults
  my_default = 0;
  alpha = 0.05;
  tail  = "both";
  ## Matlab compliant default method selection
  method_present = false;
  if length (x) < 100
    method = "exact";
  else
    method = "approximate";
  endif
  ## When called with a single input argument of second argument is empty
  if nargin == 1 || isempty (my)
    my = zeros (size (x));
  endif
  ## If second argument is a scalar convert to vector or check for Y being a
  ## vector and that X and Y have equal lengths
  if isscalar (my)
    my = repmat (my, size (x));
  elseif ! isvector (my)
    error ("signtest: Y must be either a scalar of a vector.");
  elseif numel (x) != numel (my)
    error ("signtest: X and Y vectors have different lengths.");
  endif
  ## Get optional input arguments
  i = 1;
  while (i <= length (varargin))
    switch lower (varargin{i})
      case "alpha"
        i = i + 1;
        alpha = varargin{i};
      case "tail"
        i = i + 1;
        tail = varargin{i};
      case "method"
        i = i + 1;
        method = varargin{i};
        method_present = true;
      otherwise
        error ("signtest: Invalid Name argument.");
    endswitch
    i = i + 1;
  endwhile
  ## Check values for optional input arguments
  if (! isnumeric (alpha) || isnan (alpha) || ! isscalar (alpha) ...
                          || alpha <= 0 || alpha >= 1)
    error ("signtest: alpha must be a numeric scalar in the range (0,1).");
  endif
  if ! isa (tail, 'char')
    error ("signtest: tail argument must be a string");
  elseif sum (strcmpi (tail, {"both", "right", "left"})) < 1
    error ("signtest: tail value must be either 'both', right' or 'left'.");
  endif
  if ! isa (method, 'char')
    error("signtest: method argument must be a string");
  elseif sum (strcmpi (method, {"exact", "approximate"})) < 1
    error ("signtest: method value must be either 'exact' or 'approximate'.");
  end
  ## Calculate differences between X and Y vectors: remove equal values of NaNs
  XY_diff = x(:) - my(:);
  NO_diff = (XY_diff == 0);
  XY_diff(NO_diff | isnan (NO_diff)) = [];
  ## Recalculate remaining length of X vector (after equal or NaNs removal)
  n = length (XY_diff);
  ## Check for identical X and Y input arguments
  if n == 0
    p = 1;
    h = 0;
    stats.sign = 0;
    stats.zval = NaN;
    return;
  endif
  ## Re-evaluate method selection
  if ! method_present && n < 100
    method = "exact";
  else
    method = "approximate";
  endif
  ## Get the number of positive and negative elements from X-Y differences
  pos_n = length (find (XY_diff > 0));
  neg_n = n - pos_n;
  ## Calculate stats according to selected method and tail
  switch lower(method)
    case "exact"
      stats.zval = nan;
      switch lower(tail)
        case "both"
          p = 2 * binocdf (min (neg_n, pos_n), n, 0.5);
          p = min (1, p);
        case "left"
          p = binocdf (pos_n, n, 0.5);
        case "right"
          p = binocdf (neg_n, n, 0.5);
      endswitch
      stats.zval = NaN;
    case 'approximate'
      switch lower(tail)
        case 'both'
          z_value = (pos_n - neg_n - sign (pos_n - neg_n)) / sqrt (n);
          p = 2 * normcdf (- abs (z_value));
        case 'left'
          z_value = (pos_n - neg_n + 1) / sqrt (n);
          p = normcdf (z_value);
        case 'right'
          z_value = (pos_n - neg_n - 1) / sqrt (n);
          p = normcdf (- z_value);
      endswitch
      stats.zval = z_value;
  endswitch
  stats.sign = pos_n;
  h = double (p < alpha);
endfunction

## Test suite
%!error signtest ();
%!error signtest ([]);
%!error signtest (ones(1,10), ones(1,8));
%!error signtest (ones(1,10), ones(2,10));
%!error signtest (ones(2,10), 0);
%!error signtest (ones(1,10), zeros(1,10), "alpha", 1.4)
%!error signtest (ones(1,10), zeros(1,10), "tail", "<")
%!error signtest (ones(1,10), zeros(1,10), "method", "some")
%!test
%! [pval, h, stats] = signtest ([-ones(1, 1000) 1], 0, "tail", "left");
%! assert (pval, 1.091701889420221e-218, 1e-14);
%! assert (h, 1);
%! assert (stats.zval, -31.5437631079266, 1e-14);
%!test
%! [pval, h, stats] = signtest ([-2 -1 0 2 1 3 1], 0);
%! assert (pval, 0.6875000000000006, 1e-14);
%! assert (h, 0);
%! assert (stats.zval, NaN);
%! assert (stats.sign, 4);
%!test
%! [pval, h, stats] = signtest ([-2 -1 0 2 1 3 1], 0, "method", "approximate");
%! assert (pval, 0.6830913983096086, 1e-14);
%! assert (h, 0);
%! assert (stats.zval, 0.4082482904638631, 1e-14);
%! assert (stats.sign, 4);
