## Copyright (C) 1996-2017 Kurt Hornik
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
## @deftypefn  {statistics} {@var{h} =} ztest2 (@var{x1}, @var{n1}, @var{x2}, @var{n2})
## @deftypefnx {statistics} {@var{h} =} ztest2 (@var{x1}, @var{n1}, @var{x2}, @var{n2}, @var{name}, @var{value})
## @deftypefnx {statistics} {[@var{h}, @var{pval}] =} ztest2 (@dots{})
## @deftypefnx {statistics} {[@var{h}, @var{pval}, @var{zvalue}] =} ztest2 (@dots{})
##
## Two proportion Z-test.
##
## If @var{x1} and @var{n1} are the counts of successes and trials in one
## sample, and @var{x2} and @var{n2} those in a second one, test the null
## hypothesis that the success probabilities @math{p1} and @math{p2} are the
## same.  The result is @var{h} = 0 if the null hypothesis cannot be rejected at
## the 5% significance level, or @var{h} = 1 if the null hypothesis can be
## rejected at the 5% level.
##
## Under the null, the test statistic @var{zvalue} approximately follows a
## standard normal distribution.
##
## The size of @var{h}, @var{pval}, and @var{zvalue} is the common size of
## @var{x}, @var{n1}, @var{x2}, and @var{n2}, which must be scalars or of common
## size.  A scalar input functions as a constant matrix of the same size as the
## other inputs.
##
## @code{[@var{h}, @var{pval}] = ztest2 (@dots{})} returns the p-value.  That
## is the probability of observing the given result, or one more extreme, by
## chance if the null hypothesisis true.
##
## @code{[@var{h}, @var{pval}, @var{zvalue}] = ztest2 (@dots{})} returns the
## value of the test statistic.
##
## @code{[@dots{}] = ztest2 (@dots{}, @var{name}, @var{value}, @dots{})}
## specifies one or more of the following name/value pairs:
##
## @multitable @columnfractions 0.05 0.2 0.75
## @headitem @tab Name @tab Value
## @item @tab "alpha" @tab the significance level. Default is 0.05.
##
## @item @tab "tail" @tab a string specifying the alternative hypothesis:
## @end multitable
## @multitable @columnfractions 0.1 0.25 0.65
## @item @tab "both" @tab @math{p1} is not @math{p2} (two-tailed, default)
## @item @tab "left" @tab @math{p1} is less than @math{p2} (left-tailed)
## @item @tab "right" @tab @math{p1} is greater than @math{p2} (right-tailed)
## @end multitable
##
## @seealso{chi2test, fishertest}
## @end deftypefn

function [h, pval, zvalue] = ztest2 (x1, n1, x2, n2, varargin)

  if (nargin < 4)
    print_usage ();
  endif

  if (! isscalar (x1) || ! isscalar (n1) || ! isscalar (x2) || ! isscalar (n2))
    [retval, x1, n1, x2, n2] = common_size (x1, n1, x2, n2);
    if (retval > 0)
      error ("ztest2: X1, N1, X2, and N2 must be of common size or scalars.");
    endif
  endif

  if (iscomplex (x1) || iscomplex (n1) || iscomplex (x2) || iscomplex(n2))
    error ("ztest2: X1, N1, X2, and N2 must not be complex.");
  endif

  ## Add defaults and parse optional arguments
  alpha = 0.05;
  tail = "both";
  if (nargin > 4)
    params = numel (varargin);
    if ((params / 2) != fix (params / 2))
      error ("ztest2: optional arguments must be in Name-Value pairs.")
    endif
    for idx = 1:2:params
      name = varargin{idx};
      value = varargin{idx+1};
      switch (lower (name))
        case "alpha"
          alpha = value;
          if (! isscalar (alpha) || ! isnumeric (alpha) || ...
                alpha <= 0 || alpha >= 1)
            error ("ztest2: invalid value for alpha.");
          endif
        case "tail"
          tail = value;
          if (! any (strcmpi (tail, {"both", "left", "right"})))
            error ("ztest2: invalid value for tail.");
          endif
        otherwise
          error ("ztest2: invalid name for optional arguments.");
      endswitch
    endfor
  endif

  p1 = x1 ./ n1;
  p2 = x2 ./ n2;
  pc = (x1 + x2) ./ (n1 + n2);

  zvalue  = (p1 - p2) ./ sqrt (pc .* (1 - pc) .* (1 ./ n1 + 1 ./ n2));

  cdf = stdnormal_cdf (zvalue);

  if (strcmpi (tail, "both"))
    pval = 2 * min (cdf, 1 - cdf);
  elseif (strcmpi (tail, "right"))
    pval = 1 - cdf;
  elseif (strcmpi (tail, "left"))
    pval = cdf;
  endif

  ## Determine the test outcome
  h = double (pval < alpha);
  h(isnan (pval)) = NaN;

endfunction

## Test input validation
%!error ztest2 ();
%!error ztest2 (1);
%!error ztest2 (1, 2);
%!error ztest2 (1, 2, 3);
%!error<ztest2: invalid value for alpha.> ...
%! ztest2 (1, 2, 3, 4, "alpha", 0);
%!error<ztest2: invalid value for alpha.> ...
%! ztest2 (1, 2, 3, 4, "alpha", 1.2);
%!error<ztest2: invalid value for alpha.> ...
%! ztest2 (1, 2, 3, 4, "alpha", "val");
%!error<ztest2: invalid value for tail.>  ...
%! ztest2 (1, 2, 3, 4, "tail", "val");
%!error<ztest2: invalid value for tail.>  ...
%! ztest2 (1, 2, 3, 4, "alpha", 0.01, "tail", "val");
%!error<ztest: invalid name for optional arguments.> ...
%! ztest (1, 2, 3, 4, "alpha", 0.01, "tail", "both", "badoption", 3);
