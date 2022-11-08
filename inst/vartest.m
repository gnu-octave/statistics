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
## @deftypefn {Function File} @var{h} = vartest (@var{x}, @var{v})
## @deftypefnx {Function File} @var{h} = vartest (@var{x}, @var{v}, @var{name}, @var{value})
## @deftypefnx {Function File} [@var{h}, @var{pval}] = vartest (@dots{})
## @deftypefnx {Function File} [@var{h}, @var{pval}, @var{ci}] = vartest (@dots{})
## @deftypefnx {Function File} [@var{h}, @var{pval}, @var{ci}, @var{stats}] = vartest (@dots{})
##
## One-sample test of variance.
##
## @code{@var{h} = vartest (@var{x}, @var{v})} performs a chi-square test of the
## hypothesis that the data in the vector @var{x} come from a normal
## distribution with variance @var{v}, against the alternative that @var{x}
## comes from a normal distribution with a different variance.  The result is
## @var{h} = 0 if the null hypothesis ("variance is V") cannot be rejected at
## the 5% significance level, or @var{h} = 1 if the null hypothesis can be
## rejected at the 5% level.
##
## @var{x} may also be a matrix or an N-D array.  For matrices, @code{vartest}
## performs separate tests along each column of @var{x}, and returns a vector of
## results.  For N-D arrays, @code{vartest} works along the first non-singleton
## dimension of @var{x}.  @var{v} must be a scalar.
##
## @code{vartest} treats NaNs as missing values, and ignores them.
##
## @code{[@var{h}, @var{pval}] = vartest (@dots{})} returns the p-value.  That
## is the probability of observing the given result, or one more extreme, by
## chance if the null hypothesisis true.
##
## @code{[@var{h}, @var{pval}, @var{ci}] = vartest (@dots{})} returns a
## 100 * (1 - @var{alpha})% confidence interval for the true variance.
##
## @code{[@var{h}, @var{pval}, @var{ci}, @var{stats}] = vartest (@dots{})}
## returns a structure with the following fields:
##
## @multitable @columnfractions 0.05 0.2 0.75
## @item @tab "chisqstat" @tab the value of the test statistic
## @item @tab "df" @tab the degrees of freedom of the test
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
## @seealso{ttest, ztest, kstest}
## @end deftypefn

function [h, pval, ci, stats] = vartest (x, v, varargin)

  ## Validate input arguments
  if (nargin < 2)
    error ("vartest: too few input arguments.");
  endif
  if (! isscalar (v) || ! isnumeric(v) || ! isreal(v) || v < 0)
    error ("vartest: invalid value for variance.");
  endif
  ## Add defaults
  alpha = 0.05;
  tail = "both";
  dim = [];
  if (nargin > 2)
    for idx = 3:2:nargin
      name = varargin{idx-2};
      value = varargin{idx-1};
      switch (lower (name))
        case "alpha"
          alpha = value;
          if (! isscalar (alpha) || ! isnumeric (alpha) || ...
                alpha <= 0 || alpha >= 1)
            error ("vartest: invalid value for alpha.");
          endif
        case "tail"
          tail = value;
          if (! any (strcmpi (tail, {"both", "left", "right"})))
            error ("vartest: invalid value for tail.");
          endif
        case "dim"
          dim = value;
          if (! isscalar (dim) || ! ismember (dim, 1:ndims (x)))
            error ("vartest: invalid value for operating dimension.");
          endif
        otherwise
          error ("vartest: invalid name for optional arguments.");
      endswitch
    endfor
  endif
  ## Figure out which dimension mean will work along
  if (isempty (dim))
    dim = find (size (x) != 1, 1);
  endif
  ## Replace all NaNs with zeros
  is_nan = isnan (x);
  x_dims = ndims (x);
  x(is_nan) = 0;
  ## Find sample size for each group (if more than one)
  if (any (is_nan(:)))
    sz = sum (! is_nan, dim);
  else
    sz = size (x, dim);
  endif
  ## Find degrees of freedom for each group (if more than one)
  df = max (sz - 1, 0);
  ## Calculate mean for each group (if more than one)
  x_mean = sum (x, dim) ./ max (1, sz);
  ## Center data
  if (isscalar (x_mean))
    x_centered = x - x_mean;
  else
    rep = ones (1, x_dims);
    rep(dim) = size (x, dim);
    x_centered = x - repmat (x_mean, rep);
  endif
  ## Replace all NaNs with zeros
  x_centered(is_nan) = 0;
  ## Calculate chi-square statistic
  sumsq = sum (abs (x_centered) .^ 2, dim);
  if (v > 0)
    chisqstat = sumsq ./ v;
  else
    chisqstat = Inf (size (sumsq));
    chisqstat(sumsq == 0) = NaN;
  endif
  ## Calculate p-value for the test and confidence intervals (if requested)
  if (strcmpi (tail, "both"))
    pval = chi2cdf (chisqstat, df);
    pval = 2 * min (pval, 1 - pval);
    if (nargout > 2)
      ci = cat (dim, sumsq ./ chi2inv (1 - alpha / 2, df), ...
                     sumsq ./ chi2inv (alpha / 2, df));
    endif
  elseif (strcmpi (tail, "right"))
    pval = chi2cdf (chisqstat, df);
    if (nargout > 2)
      ci = cat (dim, sumsq ./ chi2inv (1 - alpha, df), Inf (size (pval)));
    endif
  elseif (strcmpi (tail, "left"))
    pval = chi2cdf (chisqstat, df);
    if (nargout > 2)
      ci = cat (dim, zeros (size (pval)), sumsq ./ chi2inv (alpha, df));
    endif
  endif
  ## Determine the test outcome
  h = double (pval < alpha);
  h(isnan (pval)) = NaN;
  ## Create stats output structure (if requested)
  if (nargout > 3)
    stats = struct ("chisqstat", chisqstat, "df", df);
  endif

endfunction

## Test input validation
%!error<vartest: too few input arguments.> vartest ();
%!error<vartest: invalid value for variance.> vartest ([1, 2, 3, 4], -0.5);
%!error<vartest: invalid value for alpha.> ...
%! vartest ([1, 2, 3, 4], 1, "alpha", 0);
%!error<vartest: invalid value for alpha.> ...
%! vartest ([1, 2, 3, 4], 1, "alpha", 1.2);
%!error<vartest: invalid value for alpha.> ...
%! vartest ([1, 2, 3, 4], 1, "alpha", "val");
%!error<vartest: invalid value for tail.>  ...
%! vartest ([1, 2, 3, 4], 1, "tail", "val");
%!error<vartest: invalid value for tail.>  ...
%! vartest ([1, 2, 3, 4], 1, "alpha", 0.01, "tail", "val");
%!error<vartest: invalid value for operating dimension.> ...
%! vartest ([1, 2, 3, 4], 1, "dim", 3);
%!error<vartest: invalid value for operating dimension.> ...
%! vartest ([1, 2, 3, 4], 1, "alpha", 0.01, "tail", "both", "dim", 3);
%!error<vartest: invalid name for optional arguments.> ...
%! vartest ([1, 2, 3, 4], 1, 2, "alpha", 0.01, "tail", "both", "badoption", 3);
%!test
%! load carsmall
%! [h, pval, ci] = vartest (MPG, 7^2);
%! assert (h, 1);
%! assert (pval, 0.04335086742174443, 1e-14);
%! assert (ci, [49.397; 88.039], 1e-3);
