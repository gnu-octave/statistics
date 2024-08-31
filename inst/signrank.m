## Copyright (C) 2024 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {statistics} {@var{pval} =} signrank (@var{x})
## @deftypefnx {statistics} {@var{pval} =} signrank (@var{x}, @var{my})
## @deftypefnx {statistics} {@var{pval} =} signrank (@var{x}, @var{my}, @var{Name}, @var{Value})
## @deftypefnx {statistics} {[@var{pval}, @var{h}] =} signrank (@dots{})
## @deftypefnx {statistics} {[@var{pval}, @var{h}, @var{stats}] =} signrank (@dots{})
##
## Wilcoxon signed rank test for median.
##
## @code{@var{pval} = signrank (@var{x})} returns the @math{p}-value of a
## two-sided Wilcoxon signed rank test. It tests the null hypothesis that data
## in @var{x} come from a distribution with zero median at the 5% significance
## level under the assumption that the distribution is symmetric about its
## median.  @var{x} must be a vector.
##
## If the second argument @var{my} is a scalar, the null hypothesis is that
## @var{x} has median @var{my}, whereas if @var{my} is a vector, the null
## hypothesis is that the distribution of @code{@var{x} - @var{my}} has zero
## median.
##
## @code{@var{pval} = signrank (@dots{}, @var{Name}, @var{Value})} performs the
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
## is the default value for 15 of fewer observations when @qcode{"method"} is
## not specified.
##
## @item @tab @qcode{"approximate"} @tab Using normal approximation for
## computing the @math{p}-value.  It is the default value for more than 15
## observations when @qcode{"method"} is not specified.
## @end multitable
##
## @code{[@var{pval}, @var{h}] = signrank (@dots{})} also returns a logical
## value indicating the test decision.  If @var{h} is 0, the null hypothesis is
## accepted, whereas if @var{h} is 1, the null hypothesis is rejected.
##
## @code{[@var{pval}, @var{h}, @var{stats}] = signrank (@dots{})} also returns
## the structure @var{stats} containing the following fields:
##
## @multitable @columnfractions 0.18 0.02 0.8
## @headitem @var{Field} @tab @tab @var{Value}
## @item @qcode{signedrank} @tab @tab Value of the sign rank test statistic.
##
## @item @qcode{zval} @tab @tab Value of the @math{z}-statistic (only computed
## when the @qcode{"method"} is @qcode{"approximate"}).
## @end multitable
##
## @seealso{tiedrank, signtest, runstest}
## @end deftypefn

function [p, h, stats] = signrank (x, my, varargin)

  ## Check X being a vector
  if (! isvector (x))
    error ("signrank: X must be a vector.");
  endif

  ## Add defaults
  alpha = 0.05;
  tail  = "both";
  if (numel (x) <= 15)
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
    error ("signrank: Y must be either a scalar of a vector.");
  elseif (numel (x) != numel (my))
    error ("signrank: X and Y vectors have different lengths.");
  endif

  ## Get optional input arguments
  if (mod (numel (varargin), 2) != 0)
    error ("signrank: optional arguments must be in pairs.");
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
        error ("signrank: invalid Name argument.");
    endswitch
    varargin([1:2]) = [];
  endwhile

  ## Check values for optional input arguments
  if (! isnumeric (alpha) || isnan (alpha) || ! isscalar (alpha) ...
                          || alpha <= 0 || alpha >= 1)
    error ("signrank: 'alpha' must be a numeric scalar in the range 0 to 1.");
  endif
  if (! ischar (tail))
    error ("signrank: 'tail' argument must be a character vector.");
  elseif (sum (strcmpi (tail, {"both", "right", "left"})) != 1)
    error ("signrank: 'tail' value must be either 'both', right' or 'left'.");
  endif
  if (! ischar (method))
    error("signrank: 'method' argument must be a character vector.");
  elseif (sum (strcmpi (method, {"exact", "approximate"})) != 1)
    error ("signrank: 'method' value must be either 'exact' or 'approximate'.");
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
    if (n <= 15)
      method = "exact";
    else
      method = "approximate";
    endif
  endif

  ## Compute signed rank statistic
  [tie_rank, tieadj] = tiedrank (abs (XY_diff));
  w = sum (tie_rank(XY_diff > 0));
  stats.signedrank = w;

  ## Calculate stats according to selected method and tail
  switch (lower (method))

    case "exact"
      w_max = n * (n + 1) / 2;
      ## Always compute lower tail
      switch_tail = false;
      if (w > w_max / 2)
        w = w_max - w;
        switch_tail = true;
      endif
      ## Avoid integers in tied ranks
      double_ties = any (tie_rank != fix (tie_rank));
      if (double_ties)
        tie_rank = round (2 * tie_rank);
        w = round (2 * w);
      endif
      ## Loop through all combinations of ranks
      C = zeros (w + 1,1);
      C(1) = 1;
      curr = 1;
      tie_rank = sort (tie_rank);
      w_tr = tie_rank(tie_rank <= w);
      for tr = 1:numel (w_tr)
        next = min (curr + w_tr(tr), w + 1);
        C_hi = min (w_tr(tr), w + 1) + 1:next;
        C_lo = 1:length (C_hi);
        C(C_hi) = C(C_hi) + C(C_lo);
        curr = next;
      endfor
      ## Fix rank statistic
      if (double_ties)
        w = w / 2;
      endif
      ## Compute tail probability
      C = C / (2 ^ n);
      p = sum (C);
      switch (lower (tail))
        case 'both'
        p = min (1, 2 * p);  # two-sided
        case 'right'
          if (! switch_tail) # right tail is larger
            p =  1 - p + C(end);
          endif
        case 'left'
          if (switch_tail)   # left tail is larger
            p =  1 - p + C(end);
          endif
      endswitch
      stats.zval = NaN;

    case "approximate"
      ## Compute z-value
      z_nom = w - n * (n + 1) / 4;
      z_den = sqrt ((n * (n + 1) * (2 * n + 1) - tieadj) / 24);
      switch (lower (tail))
        case 'both'
          z = z_nom / z_den;
          p = 2 * normcdf (-abs (z));
        case 'right'
          z = (z_nom - 0.5) / z_den;
          p = normcdf (-z);
        case 'left'
          z = (z_nom + 0.5) / z_den;
          p = normcdf (z);
      endswitch
      stats.zval = z;

  endswitch
  h = p <= alpha;
endfunction

## Test output
%!test
%! load gradespaired.mat
%! [p, h, stats] = signrank (gradespaired(:,1), ...
%!                           gradespaired(:,2), 'tail', 'left');
%! assert (p, 0.0047, 1e-4);
%! assert (h, true);
%! assert (stats.zval, -2.5982, 1e-4);
%! assert (stats.signedrank, 2017.5);
%!test
%! load ('gradespaired.mat');
%! [p, h, stats] = signrank (gradespaired(:,1), gradespaired(:,2), ...
%!                           'tail', 'left', 'method', 'exact');
%! assert (p, 0.0045, 1e-4);
%! assert (h, true);
%! assert (stats.zval, NaN);
%! assert (stats.signedrank, 2017.5);
%!test
%! load mileage
%! [p, h, stats] = signrank (mileage(:,2), 33);
%! assert (p, 0.0312, 1e-4);
%! assert (h, true);
%! assert (stats.zval, NaN);
%! assert (stats.signedrank, 21);
%!test
%! load mileage
%! [p, h, stats] = signrank (mileage(:,2), 33, 'tail', 'right');
%! assert (p, 0.0156, 1e-4);
%! assert (h, true);
%! assert (stats.zval, NaN);
%! assert (stats.signedrank, 21);
%!test
%! load mileage
%! [p, h, stats] = signrank (mileage(:,2), 33, 'tail', 'right', ...
%!                           'alpha', 0.01, 'method', 'approximate');
%! assert (p, 0.0180, 1e-4);
%! assert (h, false);
%! assert (stats.zval, 2.0966, 1e-4);
%! assert (stats.signedrank, 21);

## Test input validation
%!error <signrank: X must be a vector.> signrank (ones (2))
%!error <signrank: Y must be either a scalar of a vector.> ...
%! signrank ([1, 2, 3, 4], ones (2))
%!error <signrank: X and Y vectors have different lengths.> ...
%! signrank ([1, 2, 3, 4], [1, 2, 3])
%!error <signrank: optional arguments must be in pairs.> ...
%! signrank ([1, 2, 3, 4], [], 'tail')
%!error <signrank: 'alpha' must be a numeric scalar in the range 0 to 1.> ...
%! signrank ([1, 2, 3, 4], [], 'alpha', 1.2)
%!error <signrank: 'alpha' must be a numeric scalar in the range 0 to 1.> ...
%! signrank ([1, 2, 3, 4], [], 'alpha', 0)
%!error <signrank: 'alpha' must be a numeric scalar in the range 0 to 1.> ...
%! signrank ([1, 2, 3, 4], [], 'alpha', -0.05)
%!error <signrank: 'alpha' must be a numeric scalar in the range 0 to 1.> ...
%! signrank ([1, 2, 3, 4], [], 'alpha', "a")
%!error <signrank: 'alpha' must be a numeric scalar in the range 0 to 1.> ...
%! signrank ([1, 2, 3, 4], [], 'alpha', [0.01, 0.05])
%!error <signrank: 'tail' argument must be a character vector.> ...
%! signrank ([1, 2, 3, 4], [], 'tail', 0.01)
%!error <signrank: 'tail' argument must be a character vector.> ...
%! signrank ([1, 2, 3, 4], [], 'tail', {"both"})
%!error <signrank: 'tail' value must be either 'both', right' or 'left'.> ...
%! signrank ([1, 2, 3, 4], [], 'tail', "some")
%!error <signrank: 'tail' value must be either 'both', right' or 'left'.> ...
%! signrank ([1, 2, 3, 4], [], 'method', 'exact', 'tail', "some")
%!error <signrank: 'method' argument must be a character vector.> ...
%! signrank ([1, 2, 3, 4], [], 'method', 0.01)
%!error <signrank: 'method' argument must be a character vector.> ...
%! signrank ([1, 2, 3, 4], [], 'method', {"exact"})
%!error <signrank: 'method' value must be either 'exact' or 'approximate'.> ...
%! signrank ([1, 2, 3, 4], [], 'method', "some")
%!error <signrank: 'method' value must be either 'exact' or 'approximate'.> ...
%! signrank ([1, 2, 3, 4], [], 'tail', "both", 'method', "some")
