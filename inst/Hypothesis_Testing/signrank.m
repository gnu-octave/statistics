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
## @multitable @columnfractions 0.18 0.8
## @headitem @var{Name} @tab @var{Value}
##
## @item @qcode{'alpha'} @tab A scalar value for the significance level of
## the test.  Default is 0.05.
##
## @item @qcode{'tail'} @tab A character vector specifying the alternative
## hypothesis.  It can take one of the following values:
## @end multitable
##
## @multitable @columnfractions 0.2 0.75
## @headitem @var{Value} @tab @var{Description}
##
## @item @qcode{'both'} @tab For one-sample test (@var{my} is empty or a
## scalar), the data in @var{x} come from a continuous distribution with median
## different than zero or @var{my}.  For two-sample test (@var{my} is a vector),
## the data in @qcode{@var{x} - @var{my}} come from a continuous distribution
## with median different than zero.
##
## @item @qcode{'left'} @tab For one-sample test (@var{my} is empty or a
## scalar), the data in @var{x} come from a continuous distribution with median
## less than zero or @var{my}.  For two-sample test (@var{my} is a vector), the
## data in @qcode{@var{x} - @var{my}} come from a continuous distribution with
## median less than zero.
##
## @item @qcode{'right'} @tab For one-sample test (@var{my} is empty or a
## scalar), the data in @var{x} come from a continuous distribution with median
## greater than zero or @var{my}.  For two-sample test (@var{my} is a vector),
## the data in @qcode{@var{x} - @var{my}} come from a continuous distribution
## with median greater than zero.
## @end multitable
##
## @multitable @columnfractions 0.18 0.8
## @headitem @var{Name} @tab @var{Value}
##
## @item @qcode{'method'} @tab A character vector specifying the method for
## computing the @math{p}-value.  It can take one of the following values:
## @end multitable
##
## @multitable @columnfractions 0.2 0.75
## @headitem @var{Value} @tab @var{Description}
##
## @item @qcode{'exact'} @tab Exact computation of the @math{p}-value.  It
## is the default value for 15 of fewer observations when @qcode{'method'} is
## not specified.
##
## @item @qcode{'approximate'} @tab Using normal approximation for
## computing the @math{p}-value.  It is the default value for more than 15
## observations when @qcode{'method'} is not specified.
## @end multitable
##
## @code{[@var{pval}, @var{h}] = signrank (@dots{})} also returns a logical
## value indicating the test decision.  If @var{h} is 0, the null hypothesis is
## accepted, whereas if @var{h} is 1, the null hypothesis is rejected.
##
## @code{[@var{pval}, @var{h}, @var{stats}] = signrank (@dots{})} also returns
## the structure @var{stats} containing the following fields:
##
## @multitable @columnfractions 0.18 0.8
## @headitem @var{Field} @tab @var{Value}
## @item @qcode{signedrank} @tab Value of the sign rank test statistic.
##
## @item @qcode{zval} @tab Value of the @math{z}-statistic (only computed
## when the @qcode{'method'} is @qcode{'approximate'}).
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
  tail  = 'both';
  if (numel (x) <= 15)
    method = 'exact';
  else
    method = 'approximate';
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
      case 'alpha'
        alpha = varargin{2};
      case 'tail'
        tail = varargin{2};
      case 'method'
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
  elseif (sum (strcmpi (tail, {'both', 'right', 'left'})) != 1)
    error ("signrank: 'tail' value must be either 'both', right' or 'left'.");
  endif
  if (! ischar (method))
    error ("signrank: 'method' argument must be a character vector.");
  elseif (sum (strcmpi (method, {'exact', 'approximate'})) != 1)
    error ("signrank: 'method' value must be either 'exact' or 'approximate'.");
  endif

  ## Calculate differences between X and Y vectors: remove equal values of NaNs.
  ## A difference smaller than the combined resolution of the two values it came
  ## from is not a real difference, so it counts as equal, and the same
  ## tolerance decides which differences rank as tied.  MATLAB defines it as
  ## EPS (X) + EPS (Y) per pair.
  XY_diff = x(:) - my(:);
  if (isfloat (x) && isfloat (my))
    epsdiff = eps (x(:)) + eps (my(:));
  else
    epsdiff = zeros (size (XY_diff));
  endif
  drop = abs (XY_diff) < epsdiff | XY_diff == 0 | isnan (XY_diff);
  XY_diff(drop) = [];
  epsdiff(drop) = [];

  ## Recalculate remaining length of X vector (after equal or NaNs removal)
  n = length (XY_diff);

  ## Check for identical X and Y input arguments
  if (n == 0)
    p = 1;
    h = 0;
    stats.signedrank = 0;
    stats.zval = [];
    return;
  endif

  ## Re-evaluate method selection
  if (! method_present)
    if (n <= 15)
      method = 'exact';
    else
      method = 'approximate';
    endif
  endif

  ## Compute signed rank statistic
  [tie_rank, tieadj] = tiedrank (abs (XY_diff), 0, 0, epsdiff);
  w = sum (tie_rank(XY_diff > 0));
  stats.signedrank = w;

  ## Calculate stats according to selected method and tail
  switch (lower (method))

    case 'exact'
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
      ## No Z-statistic exists for the exact test, so the field is created
      ## empty rather than holding a value.  R2024a omits it altogether, but
      ## R2026a and the documentation both give an empty one.
      stats.zval = [];

    case 'approximate'
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
## Field layouts below are R2024a's, measured 2026-08-17.
%!test
%! ## the exact test has no z-statistic, so no ZVAL field is created at all
%! x = [1.83 0.50 1.62 2.48 1.68 1.88 1.55 3.06 1.30];
%! y = [0.878 0.647 0.598 2.05 1.06 1.29 1.06 3.14 1.29];
%! [p, h, stats] = signrank (x, y, 'method', 'exact');
%! assert_equal (fieldnames (stats), {'signedrank'; 'zval'});
%! assert_equal (stats.signedrank, 40);
%! assert_equal (isempty (stats.zval), true);
%! assert_equal (p, 0.039062500000000, 1e-14);
%!test
%! ## the default method for a small sample is the exact one
%! x = [1.83 0.50 1.62 2.48 1.68 1.88 1.55 3.06 1.30];
%! [~, ~, stats] = signrank (x, 1);
%! assert_equal (fieldnames (stats), {'signedrank'; 'zval'});
%! assert_equal (stats.signedrank, 43);
%! assert_equal (isempty (stats.zval), true);
%!test
%! ## identical inputs give a zero statistic and no z-value
%! [p, h, stats] = signrank ([1 2 3], [1 2 3]);
%! assert_equal (p, 1);
%! assert_equal (fieldnames (stats), {'signedrank'; 'zval'});
%! assert_equal (stats.signedrank, 0);
%! assert_equal (isempty (stats.zval), true);

%!test
%! load gradespaired.mat
%! [p, h, stats] = signrank (gradespaired(:,1), ...
%!                           gradespaired(:,2), 'tail', 'left');
%! assert_equal (p, 0.0047, 1e-4);
%! assert_equal (h, true);
%! assert_equal (stats.zval, -2.5982, 1e-4);
%! assert_equal (stats.signedrank, 2017.5);
%!test
%! load ('gradespaired.mat');
%! [p, h, stats] = signrank (gradespaired(:,1), gradespaired(:,2), ...
%!                           'tail', 'left', 'method', 'exact');
%! assert_equal (p, 0.0045, 1e-4);
%! assert_equal (h, true);
%! assert_equal (isempty (stats.zval), true);
%! assert_equal (stats.signedrank, 2017.5);
%!test
%! load mileage
%! [p, h, stats] = signrank (mileage(:,2), 33);
%! assert_equal (p, 0.0312, 1e-4);
%! assert_equal (h, true);
%! assert_equal (isempty (stats.zval), true);
%! assert_equal (stats.signedrank, 21);
%!test
%! load mileage
%! [p, h, stats] = signrank (mileage(:,2), 33, 'tail', 'right');
%! assert_equal (p, 0.0156, 1e-4);
%! assert_equal (h, true);
%! assert_equal (isempty (stats.zval), true);
%! assert_equal (stats.signedrank, 21);
%!test
%! load mileage
%! [p, h, stats] = signrank (mileage(:,2), 33, 'tail', 'right', ...
%!                           'alpha', 0.01, 'method', 'approximate');
%! assert_equal (p, 0.0180, 1e-4);
%! assert_equal (h, false);
%! assert_equal (stats.zval, 2.0966, 1e-4);
%! assert_equal (stats.signedrank, 21);
%!test
%! x = [1, 2, 3, NaN, 4, 5];
%! p_clean = signrank ([1, 2, 3, 4, 5]);
%! p_nan   = signrank (x);
%! assert_equal (p_nan, p_clean);

%!test
%! ## Differences equal to within the precision of the values they came from
%! ## rank as tied, as MATLAB ranks them.  Two of these sit 2 ulps apart.
%! big = [2.1 3.4 1.2 5.6 4.3 2.2 6.7 3.3 4.4 5.5 1.1 2.9 3.8 4.9 5.1 6.2 ...
%!        7.3 2.4];
%! [p, h, stats] = signrank (big, 3, 'method', 'approximate');
%! assert_equal (stats.signedrank, 132);
%! assert_equal (stats.zval, 2.025571814690999, 1e-12);

%!test
%! ## A difference below that precision counts as no difference at all, and is
%! ## dropped exactly as an exact zero is.  Both forms measured against R2024a.
%! x = [5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20];
%! y = x + 1;
%! y(16) = x(16) + eps (x(16));
%! [p, h, stats] = signrank (x, y, 'method', 'approximate');
%! assert_equal (stats.signedrank, 0);
%! assert_equal (stats.zval, -3.872983346207417, 1e-12);
%! y(16) = x(16);
%! [p2, h2, stats2] = signrank (x, y, 'method', 'approximate');
%! assert_equal (p2, p, 1e-15);
%! assert_equal (stats2.zval, stats.zval, 1e-15);

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
%! signrank ([1, 2, 3, 4], [], 'alpha', 'a')
%!error <signrank: 'alpha' must be a numeric scalar in the range 0 to 1.> ...
%! signrank ([1, 2, 3, 4], [], 'alpha', [0.01, 0.05])
%!error <signrank: 'tail' argument must be a character vector.> ...
%! signrank ([1, 2, 3, 4], [], 'tail', 0.01)
%!error <signrank: 'tail' argument must be a character vector.> ...
%! signrank ([1, 2, 3, 4], [], 'tail', {'both'})
%!error <signrank: 'tail' value must be either 'both', right' or 'left'.> ...
%! signrank ([1, 2, 3, 4], [], 'tail', 'some')
%!error <signrank: 'tail' value must be either 'both', right' or 'left'.> ...
%! signrank ([1, 2, 3, 4], [], 'method', 'exact', 'tail', 'some')
%!error <signrank: 'method' argument must be a character vector.> ...
%! signrank ([1, 2, 3, 4], [], 'method', 0.01)
%!error <signrank: 'method' argument must be a character vector.> ...
%! signrank ([1, 2, 3, 4], [], 'method', {'exact'})
%!error <signrank: 'method' value must be either 'exact' or 'approximate'.> ...
%! signrank ([1, 2, 3, 4], [], 'method', 'some')
%!error <signrank: 'method' value must be either 'exact' or 'approximate'.> ...
%! signrank ([1, 2, 3, 4], [], 'tail', 'both', 'method', 'some')
