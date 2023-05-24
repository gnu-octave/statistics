## Copyright (C) 2012 Rik Wehbring
## Copyright (C) 1995-2016 Kurt Hornik
## Copyright (C) 2022-2023 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn {statistics} {@var{x} =} tinv (@var{p}, @var{df})
##
## Inverse of the Student's T cumulative distribution function (iCDF).
##
## For each element of @var{p}, compute the quantile (the inverse of the CDF)
## at @var{p} of the Student's T distribution with @var{df} degrees of freedom.
##
## The size of @var{p} is the common size of @var{x} and @var{df}. A scalar
## input functions as a constant matrix of the same size as the other input.
##
## This function is analogous to looking in a table for the t-value of a
## single-tailed distribution.  For very large @var{df} (>10000), the inverse of
## the standard normal distribution is used.
##
## @seealso{tcdf, tpdf, trnd, tstat}
## @end deftypefn

function x = tinv (p, df)

  if (nargin != 2)
    print_usage ();
  endif

  if (! isscalar (p) || ! isscalar (df))
    [retval, p, df] = common_size (p, df);
    if (retval > 0)
      error ("tinv: P and DF must be of common size or scalars.");
    endif
  endif

  if (iscomplex (p) || iscomplex (df))
    error ("tinv: P and DF must not be complex.");
  endif

  if (isa (p, "single") || isa (df, "single"))
    x = NaN (size (p), "single");
  else
    x = NaN (size (p));
  endif

  k = (p == 0) & (df > 0);
  x(k) = -Inf;

  k = (p == 1) & (df > 0);
  x(k) = Inf;

  if (isscalar (df))
    k = (p > 0) & (p < 1);
    if ((df > 0) && (df < 10000))
      x(k) = (sign (p(k) - 1/2)
                .* sqrt (df * (1 ./ betainv (2*min (p(k), 1 - p(k)),
                                            df/2, 1/2) - 1)));
    elseif (df >= 10000)
      ## For large df, use the quantiles of the standard normal
      x(k) = -sqrt (2) * erfcinv (2 * p(k));
    endif
  else
    k = (p > 0) & (p < 1) & (df > 0) & (df < 10000);
    x(k) = (sign (p(k) - 1/2)
              .* sqrt (df(k) .* (1 ./ betainv (2*min (p(k), 1 - p(k)),
                                              df(k)/2, 1/2) - 1)));

    ## For large df, use the quantiles of the standard normal
    k = (p > 0) & (p < 1) & (df >= 10000);
    x(k) = -sqrt (2) * erfcinv (2 * p(k));
  endif

endfunction


%!shared p
%! p = [-1 0 0.5 1 2];
%!assert (tinv (p, ones (1,5)), [NaN -Inf 0 Inf NaN])
%!assert (tinv (p, 1), [NaN -Inf 0 Inf NaN], eps)
%!assert (tinv (p, [1 0 NaN 1 1]), [NaN NaN NaN Inf NaN], eps)
%!assert (tinv ([p(1:2) NaN p(4:5)], 1), [NaN -Inf NaN Inf NaN])

## Test class of input preserved
%!assert (tinv ([p, NaN], 1), [NaN -Inf 0 Inf NaN NaN], eps)
%!assert (tinv (single ([p, NaN]), 1), single ([NaN -Inf 0 Inf NaN NaN]), eps ("single"))
%!assert (tinv ([p, NaN], single (1)), single ([NaN -Inf 0 Inf NaN NaN]), eps ("single"))

## Test input validation
%!error tinv ()
%!error tinv (1)
%!error tinv (1,2,3)
%!error tinv (ones (3), ones (2))
%!error tinv (ones (2), ones (3))
%!error tinv (i, 2)
%!error tinv (2, i)
