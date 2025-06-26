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
## @deftypefn  {statistics} {@var{h} =} ztest (@var{x}, @var{m}, @var{sigma})
## @deftypefnx {statistics} {@var{h} =} ztest (@var{x}, @var{m}, @var{sigma}, @var{Name}, @var{Value})
## @deftypefnx {statistics} {[@var{h}, @var{pval}] =} ztest (@dots{})
## @deftypefnx {statistics} {[@var{h}, @var{pval}, @var{ci}] =} ztest (@dots{})
## @deftypefnx {statistics} {[@var{h}, @var{pval}, @var{ci}, @var{zvalue}] =} ztest (@dots{})
##
## One-sample Z-test.
##
## @code{@var{h} = ztest (@var{x}, @var{v})} performs a Z-test of the hypothesis
## that the data in the vector @var{x} come from a normal distribution with mean
## @var{m}, against the alternative that @var{x} comes from a normal
## distribution with a different mean @var{m}.  The result is @var{h} = 0 if the
## null hypothesis ("mean is M") cannot be rejected at the 5% significance
## level, or @var{h} = 1 if the null hypothesis can be rejected at the 5% level.
##
## @var{x} may also be a matrix or an N-D array.  For matrices, @code{ztest}
## performs separate tests along each column of @var{x}, and returns a vector of
## results.  For N-D arrays, @code{ztest} works along the first non-singleton
## dimension of @var{x}.  @var{m} and @var{sigma} must be a scalars.
##
## @code{ztest} treats NaNs as missing values, and ignores them.
##
## @code{[@var{h}, @var{pval}] = ztest (@dots{})} returns the p-value.  That
## is the probability of observing the given result, or one more extreme, by
## chance if the null hypothesis true.
##
## @code{[@var{h}, @var{pval}, @var{ci}] = ztest (@dots{})} returns a
## 100 * (1 - @var{alpha})% confidence interval for the true mean.
##
## @code{[@var{h}, @var{pval}, @var{ci}, @var{zvalue}] = ztest (@dots{})}
## returns the value of the test statistic.
##
## @code{[@dots{}] = ztest (@dots{}, @var{Name}, @var{Value}, @dots{})}
## specifies one or more of the following @var{Name}/@var{Value} pairs:
##
## @multitable @columnfractions 0.05 0.2 0.75
## @headitem @tab @var{Name} @tab @var{Value}
## @item @tab "alpha" @tab the significance level. Default is 0.05.
##
## @item @tab "dim" @tab dimension to work along a matrix or an N-D array.
##
## @item @tab "tail" @tab a string specifying the alternative hypothesis:
## @end multitable
## @multitable @columnfractions 0.1 0.15 0.75
## @item @tab "both" @tab "mean is not @var{m}" (two-tailed, default)
## @item @tab "left" @tab "mean is less than @var{m}" (left-tailed)
## @item @tab "right" @tab "mean is greater than @var{m}" (right-tailed)
## @end multitable
##
## @seealso{ttest, vartest, signtest, kstest}
## @end deftypefn

function [h, pval, ci, zvalue] = ztest (x, m, sigma, varargin)

  ## Validate input arguments
  if (nargin < 3)
    error ("ztest: too few input arguments.");
  endif
  if (! isscalar (m) || ! isnumeric(m) || ! isreal(m))
    error ("ztest: invalid value for mean.");
  endif
  if (! isscalar (sigma) || ! isnumeric(sigma) || ! isreal(sigma) || sigma < 0)
    error ("ztest: invalid value for standard deviation.");
  endif
  ## Add defaults
  alpha = 0.05;
  tail = "both";
  dim = [];
  if (nargin > 3)
    for idx = 4:2:nargin
      name = varargin{idx-3};
      value = varargin{idx-2};
      switch (lower (name))
        case "alpha"
          alpha = value;
          if (! isscalar (alpha) || ! isnumeric (alpha) || ...
                alpha <= 0 || alpha >= 1)
            error ("ztest: invalid VALUE for alpha.");
          endif
        case "tail"
          tail = value;
          if (! any (strcmpi (tail, {"both", "left", "right"})))
            error ("ztest: invalid VALUE for tail.");
          endif
        case "dim"
          dim = value;
          if (! isscalar (dim) || ! ismember (dim, 1:ndims (x)))
            error ("ztest: invalid VALUE for operating dimension.");
          endif
        otherwise
          error ("ztest: invalid NAME for optional arguments.");
      endswitch
    endfor
  endif
  ## Figure out which dimension mean will work along
  if (isempty (dim))
    dim = find (size (x) != 1, 1);
  endif
  ## Replace all NaNs with zeros
  is_nan = isnan (x);
  ## Find sample size for each group (if more than one)
  if (any (is_nan(:)))
    sz = sum (! is_nan, dim);
  else
    sz = size (x, dim);
  endif
  ## Calculate mean, strandard error and z-value for each group
  x_mean = sum (x(! is_nan), dim) ./ max (1, sz);
  stderr = sigma ./ sqrt (sz);
  zvalue = (x_mean - m) ./ stderr;
  ## Calculate p-value for the test and confidence intervals (if requested)
  if (strcmpi (tail, "both"))
    pval = 2 * normcdf (- abs (zvalue), 0, 1);
    if (nargout > 2)
      crit = norminv (1 - alpha / 2, 0, 1) .* stderr;
      ci = cat (dim, x_mean - crit, x_mean + crit);
    endif
  elseif (strcmpi (tail, "right"))
    p = normcdf (- zvalue,0,1);
    if (nargout > 2)
      crit = norminv (1 - alpha, 0, 1) .* stderr;
      ci = cat (dim, x_mean - crit, Inf (size (p)));
    endif
  elseif (strcmpi (tail, "left"))
    p = normcdf (zvalue, 0, 1);
    if (nargout > 2)
      crit = norminv (1 - alpha, 0, 1) .* stderr;
      ci = cat (dim, - Inf (size (p)), x_mean + crit);
    endif
  endif
  ## Determine the test outcome
  h = double (pval < alpha);
  h(isnan (pval)) = NaN;

endfunction

## Test input validation
%!error<ztest: too few input arguments.> ztest ();
%!error<ztest: invalid value for standard deviation.> ...
%! ztest ([1, 2, 3, 4], 2, -0.5);
%!error<ztest: invalid VALUE for alpha.> ...
%! ztest ([1, 2, 3, 4], 1, 2, "alpha", 0);
%!error<ztest: invalid VALUE for alpha.> ...
%! ztest ([1, 2, 3, 4], 1, 2, "alpha", 1.2);
%!error<ztest: invalid VALUE for alpha.> ...
%! ztest ([1, 2, 3, 4], 1, 2, "alpha", "val");
%!error<ztest: invalid VALUE for tail.>  ...
%! ztest ([1, 2, 3, 4], 1, 2, "tail", "val");
%!error<ztest: invalid VALUE for tail.>  ...
%! ztest ([1, 2, 3, 4], 1, 2, "alpha", 0.01, "tail", "val");
%!error<ztest: invalid VALUE for operating dimension.> ...
%! ztest ([1, 2, 3, 4], 1, 2, "dim", 3);
%!error<ztest: invalid VALUE for operating dimension.> ...
%! ztest ([1, 2, 3, 4], 1, 2, "alpha", 0.01, "tail", "both", "dim", 3);
%!error<ztest: invalid NAME for optional arguments.> ...
%! ztest ([1, 2, 3, 4], 1, 2, "alpha", 0.01, "tail", "both", "badoption", 3);
## Test results
%!test
%! load carsmall
%! [h, pval, ci] = ztest (MPG, mean (MPG, "omitnan"), std (MPG, "omitnan"));
%! assert (h, 0);
%! assert (pval, 1, 1e-14);
%! assert (ci, [22.094; 25.343], 1e-3);
%!test
%! load carsmall
%! [h, pval, ci] = ztest (MPG, 26, 8);
%! assert (h, 1);
%! assert (pval, 0.00568359158544743, 1e-14);
%! assert (ci, [22.101; 25.335], 1e-3);
%!test
%! load carsmall
%! [h, pval, ci] = ztest (MPG, 26, 4);
%! assert (h, 1);
%! assert (pval, 3.184168011941316e-08, 1e-14);
%! assert (ci, [22.909; 24.527], 1e-3);
