## Copyright (C) 2014 Tony Richardson
## Copyright (C) 2022-2024 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {statistics} {@var{pval} =} signtest (@var{x})
## @deftypefnx {statistics} {@var{pval} =} signtest (@var{x}, @var{my})
## @deftypefnx {statistics} {@var{pval} =} signtest (@var{x}, @var{my}, @var{Name}, @var{Value})
## @deftypefnx {statistics} {[@var{pval}, @var{h}] =} signtest (@dots{})
## @deftypefnx {statistics} {[@var{pval}, @var{h}, @var{stats}] =} signtest (@dots{})
##
## Signed test for median.
##
## @code{@var{pval} = signtest (@var{x})} returns the @math{p}-value of a
## two-sided sign test. It tests the null hypothesis that data in @var{x} come
## from a distribution with zero median at the 5% significance level.  @var{x}
## must be a vector.
##
## If the second argument @var{my} is a scalar, the null hypothesis is that
## @var{x} has median @var{my}, whereas if @var{my} is a vector, the null
## hypothesis is that the distribution of @code{@var{x} - @var{my}} has zero
## median.
##
## @code{@var{pval} = signtest (@dots{}, @var{Name}, @var{Value})} performs the
## Wilcoxon signed rank test with additional options specified by one or more of
## the following @var{Name}, @var{Value} pair arguments:
##
## @multitable @columnfractions 0.18 0.02 0.8
## @headitem @var{Name} @tab @tab @var{Value}
##
## @item @qcode{"alpha"} @tab @tab A scalar value for the significance level of
## the test.  Default is 0.05.
##
## @item @qcode{"tail"} @tab @tab A character vector specifying the alternative
## hypothesis.  It can take one of the following values:
## @end multitable
##
## @multitable @columnfractions 0.05 0.2 0.75
## @headitem @tab @var{Value} @tab @var{Description}
##
## @item @tab @qcode{"both"} @tab For one-sample test (@var{my} is empty or a
## scalar), the data in @var{x} come from a continuous distribution with median
## different than zero or @var{my}.  For two-sample test (@var{my} is a vector),
## the data in @qcode{@var{x} - @var{my}} come from a continuous distribution
## with median different than zero.
##
## @item @tab @qcode{"left"} @tab For one-sample test (@var{my} is empty or a
## scalar), the data in @var{x} come from a continuous distribution with median
## less than zero or @var{my}.  For two-sample test (@var{my} is a vector), the
## data in @qcode{@var{x} - @var{my}} come from a continuous distribution with
## median less than zero.
##
## @item @tab @qcode{"right"} @tab For one-sample test (@var{my} is empty or a
## scalar), the data in @var{x} come from a continuous distribution with median
## greater than zero or @var{my}.  For two-sample test (@var{my} is a vector),
## the data in @qcode{@var{x} - @var{my}} come from a continuous distribution
## with median greater than zero.
## @end multitable
##
## @multitable @columnfractions 0.18 0.02 0.8
## @headitem @var{Name} @tab @tab @var{Value}
##
## @item @qcode{"method"} @tab @tab A character vector specifying the method for
## computing the @math{p}-value.  It can take one of the following values:
## @end multitable
##
## @multitable @columnfractions 0.05 0.2 0.75
## @headitem @tab @var{Value} @tab @var{Description}
##
## @item @tab @qcode{"exact"} @tab Exact computation of the @math{p}-value.  It
## is the default value for fewer than 100 observations when @qcode{"method"} is
## not specified.
##
## @item @tab @qcode{"approximate"} @tab Using normal approximation for
## computing the @math{p}-value.  It is the default value for 100 or more
## observations when @qcode{"method"} is not specified.
## @end multitable
##
## @code{[@var{pval}, @var{h}] = signtest (@dots{})} also returns a logical
## value indicating the test decision.  If @var{h} is 0, the null hypothesis is
## accepted, whereas if @var{h} is 1, the null hypothesis is rejected.
##
## @code{[@var{pval}, @var{h}, @var{stats}] = signtest (@dots{})} also returns
## the structure @var{stats} containing the following fields:
##
## @multitable @columnfractions 0.18 0.02 0.8
## @headitem @var{Field} @tab @tab @var{Value}
## @item @qcode{sign} @tab @tab Value of the sign test statistic.
##
## @item @qcode{zval} @tab @tab Value of the @math{z}-statistic (only computed
## when the @qcode{"method"} is @qcode{"approximate"}).
## @end multitable
##
## @seealso{signrank, tiedrank, runstest}
## @end deftypefn

function [p, h, stats] = signtest (x, my, varargin)

  ## Check X being a vector
  if (! isvector (x))
    error ("signtest: X must be a vector.");
  endif

  ## Add defaults
  alpha = 0.05;
  tail  = "both";
  if (numel (x) < 100)
    method = "exact";
  else
    method = "approximate";
  endif
  method_present = false;

  ## When called with a single input argument of second argument is empty
  if (nargin == 1 || isempty (my))
    my = zeros (size (x));
  endif

  ## If second argument is a scalar convert to vector or check for Y being a
  ## vector and that X and Y have equal lengths
  if (isscalar (my))
    my = repmat (my, size (x));
  elseif (! isvector (my))
    error ("signtest: Y must be either a scalar of a vector.");
  elseif (numel (x) != numel (my))
    error ("signtest: X and Y vectors have different lengths.");
  endif

  ## Get optional input arguments
  if (mod (numel (varargin), 2) != 0)
    error ("signtest: optional arguments must be in pairs.");
  endif
  while (numel (varargin) > 0)
    switch (lower (varargin{1}))
      case "alpha"
        alpha = varargin{2};
      case "tail"
        tail = varargin{2};
      case "method"
        method = varargin{2};
        method_present = true;
      otherwise
        error ("signtest: invalid Name argument.");
    endswitch
    varargin([1:2]) = [];
  endwhile

  ## Check values for optional input arguments
  if (! isnumeric (alpha) || isnan (alpha) || ! isscalar (alpha) ...
                          || alpha <= 0 || alpha >= 1)
    error ("signtest: 'alpha' must be a numeric scalar in the range 0 to 1.");
  endif
  if (! ischar (tail))
    error ("signtest: 'tail' argument must be a character vector.");
  elseif (sum (strcmpi (tail, {"both", "right", "left"})) != 1)
    error ("signtest: 'tail' value must be either 'both', right' or 'left'.");
  endif
  if (! ischar (method))
    error("signtest: 'method' argument must be a character vector.");
  elseif (sum (strcmpi (method, {"exact", "approximate"})) != 1)
    error ("signtest: 'method' value must be either 'exact' or 'approximate'.");
  endif

  ## Calculate differences between X and Y vectors: remove equal values of NaNs
  XY_diff = x(:) - my(:);
  NO_diff = (XY_diff == 0);
  XY_diff(NO_diff | isnan (NO_diff)) = [];

  ## Recalculate remaining length of X vector (after equal or NaNs removal)
  n = length (XY_diff);

  ## Check for identical X and Y input arguments
  if (n == 0)
    p = 1;
    h = 0;
    stats.sign = 0;
    stats.zval = NaN;
    return;
  endif

  ## Re-evaluate method selection
  if (! method_present)
    if (n < 100)
      method = "exact";
    else
      method = "approximate";
    endif
  endif

  ## Get the number of positive and negative elements from X-Y differences
  pos_n = length (find (XY_diff > 0));
  neg_n = n - pos_n;

  ## Calculate stats according to selected method and tail
  switch (lower (method))

    case "exact"
      switch (lower (tail))
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
      switch (lower (tail))
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

## Test output
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

## Test input validation
%!error <signtest: X must be a vector.> signtest (ones (2))
%!error <signtest: Y must be either a scalar of a vector.> ...
%! signtest ([1, 2, 3, 4], ones (2))
%!error <signtest: X and Y vectors have different lengths.> ...
%! signtest ([1, 2, 3, 4], [1, 2, 3])
%!error <signtest: optional arguments must be in pairs.> ...
%! signtest ([1, 2, 3, 4], [], 'tail')
%!error <signtest: 'alpha' must be a numeric scalar in the range 0 to 1.> ...
%! signtest ([1, 2, 3, 4], [], 'alpha', 1.2)
%!error <signtest: 'alpha' must be a numeric scalar in the range 0 to 1.> ...
%! signtest ([1, 2, 3, 4], [], 'alpha', 0)
%!error <signtest: 'alpha' must be a numeric scalar in the range 0 to 1.> ...
%! signtest ([1, 2, 3, 4], [], 'alpha', -0.05)
%!error <signtest: 'alpha' must be a numeric scalar in the range 0 to 1.> ...
%! signtest ([1, 2, 3, 4], [], 'alpha', "a")
%!error <signtest: 'alpha' must be a numeric scalar in the range 0 to 1.> ...
%! signtest ([1, 2, 3, 4], [], 'alpha', [0.01, 0.05])
%!error <signtest: 'tail' argument must be a character vector.> ...
%! signtest ([1, 2, 3, 4], [], 'tail', 0.01)
%!error <signtest: 'tail' argument must be a character vector.> ...
%! signtest ([1, 2, 3, 4], [], 'tail', {"both"})
%!error <signtest: 'tail' value must be either 'both', right' or 'left'.> ...
%! signtest ([1, 2, 3, 4], [], 'tail', "some")
%!error <signtest: 'tail' value must be either 'both', right' or 'left'.> ...
%! signtest ([1, 2, 3, 4], [], 'method', 'exact', 'tail', "some")
%!error <signtest: 'method' argument must be a character vector.> ...
%! signtest ([1, 2, 3, 4], [], 'method', 0.01)
%!error <signtest: 'method' argument must be a character vector.> ...
%! signtest ([1, 2, 3, 4], [], 'method', {"exact"})
%!error <signtest: 'method' value must be either 'exact' or 'approximate'.> ...
%! signtest ([1, 2, 3, 4], [], 'method', "some")
%!error <signtest: 'method' value must be either 'exact' or 'approximate'.> ...
%! signtest ([1, 2, 3, 4], [], 'tail', "both", 'method', "some")
