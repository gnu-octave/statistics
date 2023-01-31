## Copyright (C) 2014 Tony Richardson
## Copyright (C) 2022 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {statistics} @var{h} = vartest2 (@var{x}, @var{y})
## @deftypefnx {statistics} @var{h} = vartest2 (@var{x}, @var{y}, @var{name}, @var{value})
## @deftypefnx {statistics} [@var{h}, @var{pval}] = vartest2 (@dots{})
## @deftypefnx {statistics} [@var{h}, @var{pval}, @var{ci}] = vartest2 (@dots{})
## @deftypefnx {statistics} [@var{h}, @var{pval}, @var{ci}, @var{stats}] = vartest2 (@dots{})
##
## Two-sample F test for equal variances.
##
## @code{@var{h} = vartest2 (@var{x}, @var{y})} performs an F test of the
## hypothesis that the independent data in vectors @var{x} and @var{y} come from
## normal distributions with equal variance, against the alternative that they
## come from normal distributions with different variances.  The result is
## @var{h} = 0 if the null hypothesis ("variance are equal") cannot be rejected
## at the 5% significance level, or @var{h} = 1 if the null hypothesis can be
## rejected at the 5% level.
##
## @var{x} and @var{y} may also be matrices or N-D arrays.  For matrices,
## @code{vartest2} performs separate tests along each column and returns a
## vector of results.  For N-D arrays, @code{vartest2} works along the first
## non-singleton dimension and @var{x} and @var{y} must have the same size along
## all the remaining dimensions.
##
## @code{vartest} treats NaNs as missing values, and ignores them.
##
## @code{[@var{h}, @var{pval}] = vartest (@dots{})} returns the p-value.  That
## is the probability of observing the given result, or one more extreme, by
## chance if the null hypothesisis true.
##
## @code{[@var{h}, @var{pval}, @var{ci}] = vartest (@dots{})} returns a 100 *
## (1 - @var{alpha})% confidence interval for the true ratio var(X)/var(Y).
##
## @code{[@var{h}, @var{pval}, @var{ci}, @var{stats}] = vartest (@dots{})}
## returns a structure with the following fields:
##
## @multitable @columnfractions 0.05 0.2 0.75
## @item @tab "fstat" @tab the value of the test statistic
## @item @tab "df1" @tab the numerator degrees of freedom of the test
## @item @tab "df2" @tab the denominator degrees of freedom of the test
## @end multitable
##
## @code{[@dots{}] = vartest (@dots{}, @var{name}, @var{value}), @dots{}}
## specifies one or more of the following name/value pairs:
##
## @multitable @columnfractions 0.05 0.2 0.75
## @headitem @tab Name @tab Value
## @item @tab "alpha" @tab the significance level. Default is 0.05.
##
## @item @tab "dim" @tab dimension to work along a matrix or an N-D array.
##
## @item @tab "tail" @tab a string specifying the alternative hypothesis:
## @end multitable
## @multitable @columnfractions 0.1 0.15 0.75
## @item @tab "both" @tab "variance is not @var{v}" (two-tailed, default)
## @item @tab "left" @tab "variance is less than @var{v}" (left-tailed)
## @item @tab "right" @tab "variance is greater than @var{v}" (right-tailed)
## @end multitable
##
## @seealso{ttest2, kstest2, bartlett_test, levene_test}
## @end deftypefn

function [h, pval, ci, stats] = vartest2 (x, y, varargin)

  ## Validate input arguments
  if (nargin < 2)
    error ("vartest2: too few input arguments.");
  endif
  if (isscalar (x) || isscalar(y))
    error ("vartest2: X and Y must be vectors or matrices or N-D arrays.");
  endif
  ## If X and Y are vectors make them the same orientation
  if (isvector (x) && isvector (y))
    if (size (x, 1) == 1)
        y = y(:)';
    else
        y = y(:);
    endif
  endif
  ## Add defaults
  alpha = 0.05;
  tail = "both";
  dim = [];
  if (nargin > 2 && mod (numel (varargin(:)), 2) == 0)
    for idx = 3:2:nargin
      name = varargin{idx-2};
      value = varargin{idx-1};
      switch (lower (name))
        case "alpha"
          alpha = value;
          if (! isscalar (alpha) || ! isnumeric (alpha) || ...
                alpha <= 0 || alpha >= 1)
            error ("vartest2: invalid value for alpha.");
          endif
        case "tail"
          tail = value;
          if (! any (strcmpi (tail, {"both", "left", "right"})))
            error ("vartest2: invalid value for tail.");
          endif
        case "dim"
          dim = value;
          if (! isscalar (dim) || ! ismember (dim, 1:ndims (x)))
            error ("vartest2: invalid value for operating dimension.");
          endif
        otherwise
          error ("vartest2: invalid name for optional arguments.");
      endswitch
    endfor
  elseif (nargin > 2 && mod (numel (varargin(:)), 2) != 0)
    error ("vartest2: optional arguments must be in name/value pairs.");
  endif
  ## Figure out which dimension mean will work along
  if (isempty (dim))
    dim = find (size (x) != 1, 1);
  endif
  ## Check that all non-working dimensions of X and Y are of equal size
  x_size = size (x);
  y_size = size (y);
  x_size(dim) = 1;
  y_size(dim) = 1;
  if (! isequal (x_size, y_size))
    error ("vartestt2: input size mismatch.");
  endif
  ## Compute statistics for each sample
  [df1, x_var] = getstats(x,dim);
  [df2, y_var] = getstats(y,dim);
  ## Compute F statistic
  F = NaN (size (x_var));
  t1 = (y_var > 0);
  F(t1) = x_var(t1) ./ y_var(t1);
  t2 = (x_var > 0) & ! t1;
  F(t2) = Inf;
  ## Calculate p-value for the test and confidence intervals (if requested)
  if (strcmpi (tail, "both"))
    pval = 2 * min (fcdf (F, df1, df2), 1 - fcdf (F, df1, df2));
    if (nargout > 2)
      ci = cat (dim, F .* finv (alpha / 2, df2, df1), ...
                     F ./ finv (alpha / 2, df1, df2));
    endif
  elseif (strcmpi (tail, "right"))
    pval = 1 - fcdf (F, df1, df2);
    if (nargout > 2)
      ci = cat (dim, F .* finv (alpha, df2, df1), Inf (size (F)));
    endif
  elseif (strcmpi (tail, "left"))
    pval = fcdf (F, df1, df2);
    if (nargout > 2)
      ci = cat (dim, zeros (size (F)), F ./ finv (alpha, df1, df2));
    endif
  endif
  ## Determine the test outcome
  h = double (pval < alpha);
  h(isnan (pval)) = NaN;
  ## Create stats output structure (if requested)
  if (nargout > 3)
    stats = struct ("fstat", F, "df1", df1, "df2", df2);
  endif

endfunction

## Compute statistics for one sample
function [df, data_var] = getstats (data, dim)
  ## Calculate sample size and df by ignoring NaNs
  is_nan = isnan (data);
  n_data = sum (! is_nan, dim);
  df = max (n_data - 1, 0);
  ## Calculate mean
  data(is_nan) = 0;
  m_data = sum (data, dim) ./ max (1, n_data);
  ## Calculate variance
  if (isscalar (m_data))
     c_data = data - m_data;
  else
     rep = ones (1, ndims (data));
     rep(dim) = size (data, dim);
     c_data = data - repmat (m_data, rep);
  end
  c_data(is_nan) = 0;
  data_var = sum (abs (c_data) .^ 2,dim);
  t = (df > 0);
  data_var(t) = data_var(t) ./ df(t);
  data_var(! t) = NaN;
  ## Make df a scalar if possible
  if (numel (df) > 1 && all (df(:) == df(1)))
     df = df(1);
  end
endfunction

## Test input validation
%!error<vartest2: too few input arguments.> vartest2 ();
%!error<vartest2: too few input arguments.> vartest2 (ones (20,1));
%!error<vartest2: X and Y must be vectors or matrices or N-D arrays.> ...
%! vartest2 (rand (20,1), 5);
%!error<vartest2: invalid value for alpha.> ...
%! vartest2 (rand (20,1), rand (25,1)*2, "alpha", 0);
%!error<vartest2: invalid value for alpha.> ...
%! vartest2 (rand (20,1), rand (25,1)*2, "alpha", 1.2);
%!error<vartest2: invalid value for alpha.> ...
%! vartest2 (rand (20,1), rand (25,1)*2, "alpha", "some");
%!error<vartest2: invalid value for alpha.> ...
%! vartest2 (rand (20,1), rand (25,1)*2, "alpha", [0.05, 0.001]);
%!error<vartest2: invalid value for tail.> ...
%! vartest2 (rand (20,1), rand (25,1)*2, "tail", [0.05, 0.001]);
%!error<vartest2: invalid value for tail.> ...
%! vartest2 (rand (20,1), rand (25,1)*2, "tail", "some");
%!error<vartest2: invalid value for operating dimension.> ...
%! vartest2 (rand (20,1), rand (25,1)*2, "dim", 3);
%!error<vartest2: invalid value for operating dimension.> ...
%! vartest2 (rand (20,1), rand (25,1)*2, "alpha", 0.001, "dim", 3);
%!error<vartest2: invalid name for optional arguments.> ...
%! vartest2 (rand (20,1), rand (25,1)*2, "some", 3);
%!error<vartest2: optional arguments must be in name/value pairs.> ...
%! vartest2 (rand (20,1), rand (25,1)*2, "some");
## Test results
%!test
%! load carsmall
%! [h, pval, ci, stat] = vartest2 (MPG(Model_Year==82), MPG(Model_Year==76));
%! assert (h, 0);
%! assert (pval, 0.6288022362718455, 1e-13);
%! assert (ci, [0.4139; 1.7193], 1e-4);
%! assert (stat.fstat, 0.8384, 1e-4);
%! assert (stat.df1, 30);
%! assert (stat.df2, 33);

