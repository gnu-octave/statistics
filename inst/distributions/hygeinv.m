## Copyright (C) 2012 Rik Wehbring
## Copyright (C) 1997-2016 Kurt Hornik
## Copyright (C) 2022 Nicholas R. Jankowski
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
## @deftypefn  {statistics} @var{p} = hygeinv (@var{x}, @var{t}, @var{m}, @var{n})
##
## Inverse of the hypergeometric cumulative distribution function (iCDF).
##
## For each element of @var{x}, compute the quantile (the inverse of the CDF)
## at @var{x} of the hypergeometric distribution with parameters @var{t},
## @var{m}, and @var{n}.  The size of @var{p} is the common size of the input
## parameters.  A scalar input functions as a constant matrix of the same size
## as the other inputs.
##
## This is the probability of obtaining @var{x} marked items when randomly
## drawing a sample of size @var{n} without replacement from a population of
## total size @var{t} containing @var{m} marked items.  The parameters @var{t},
## @var{m}, and @var{n} must be positive integers with @var{m} and @var{n} not
## greater than @var{t}.
##
## @seealso{hygeinv, hygepdf, hygernd, hygestat}
## @end deftypefn

function p = hygeinv (x, t, m, n)

  if (nargin != 4)
    print_usage ();
  endif

  if (! isscalar (x) || ! isscalar (t) || ! isscalar (m) || ! isscalar (n))
    [retval, x, t, m, n] = common_size (x, t, m, n);
    if (retval > 0)
      error ("hygeinv: X, T, M, and N must be of common size or scalars.");
    endif
  endif

  if (iscomplex (x) || iscomplex (t) || iscomplex (m) || iscomplex (n))
    error ("hygeinv: X, T, M, and N must not be complex.");
  endif

  if (isa (x, "single") || isa (t, "single")
      || isa (m, "single") || isa (n, "single"))
    p = NaN (size (x), "single");
  else
    p = NaN (size (x));
  endif

  ok = ((t >= 0) & (m >= 0) & (n > 0) & (m <= t) & (n <= t) &
        (t == fix (t)) & (m == fix (m)) & (n == fix (n)));

  if (isscalar (t))
    if (ok)
      p = discrete_inv (x, 0 : n, hygepdf (0 : n, t, m, n));
      p(x == 0) = 0;  # Hack to return correct value for start of distribution
    endif
  else
    k = (x == 0);
    p (ok & k) = 0; # set any x=0 to 0 if not already set to output NaN
    k = (x == 1);
    p (ok & k) = n(ok & k);
    ok &= (x>0 & x<1); #remove 0's and x's outside (0,1), leave unfilled as NaN

    if (any (ok(:)))
      n = n(ok);
      v = 0 : max (n(:));

      ## manually perform discrete_inv to enable vectorizing with array input
      p_tmp = cumsum (hygepdf (v, t(ok), m(ok), n, "vectorexpand"), 2);
      sz_p = size (p_tmp);
      end_locs = sub2ind (sz_p, [1 : numel(n)]', n(:) + 1);

      ## manual row-wise vectorization of lookup, which returns index of element
      ## less than or equal to test value, zero if test value less than lowest
      ## number, and max index if greater than highest number. operated on
      ## flipped p_tmp, adjusting for different vector lengths in array rows.
      p_tmp = (p_tmp ./ p_tmp(end_locs))(:, end:-1:1) - x(ok)(:);
      p_tmp(p_tmp>=0) = NaN;
      [p_match, p_match_idx] = max (p_tmp, [], 2);
      p_match_idx(isnan(p_match)) = v(end) + 2;

      p(ok) = v(v(end) - p_match_idx + 3);

    endif
  endif
endfunction


%!shared x
%! x = [-1 0 0.5 1 2];
%!assert (hygeinv (x, 4*ones (1,5), 2*ones (1,5), 2*ones (1,5)), [NaN 0 1 2 NaN])
%!assert (hygeinv (x, 4*ones (1,5), 2, 2), [NaN 0 1 2 NaN])
%!assert (hygeinv (x, 4, 2*ones (1,5), 2), [NaN 0 1 2 NaN])
%!assert (hygeinv (x, 4, 2, 2*ones (1,5)), [NaN 0 1 2 NaN])
%!assert (hygeinv (x, 4*[1 -1 NaN 1.1 1], 2, 2), [NaN NaN NaN NaN NaN])
%!assert (hygeinv (x, 4, 2*[1 -1 NaN 1.1 1], 2), [NaN NaN NaN NaN NaN])
%!assert (hygeinv (x, 4, 5, 2), [NaN NaN NaN NaN NaN])
%!assert (hygeinv (x, 4, 2, 2*[1 -1 NaN 1.1 1]), [NaN NaN NaN NaN NaN])
%!assert (hygeinv (x, 4, 2, 5), [NaN NaN NaN NaN NaN])
%!assert (hygeinv ([x(1:2) NaN x(4:5)], 4, 2, 2), [NaN 0 NaN 2 NaN])

## Test class of input preserved
%!assert (hygeinv ([x, NaN], 4, 2, 2), [NaN 0 1 2 NaN NaN])
%!assert (hygeinv (single ([x, NaN]), 4, 2, 2), single ([NaN 0 1 2 NaN NaN]))
%!assert (hygeinv ([x, NaN], single (4), 2, 2), single ([NaN 0 1 2 NaN NaN]))
%!assert (hygeinv ([x, NaN], 4, single (2), 2), single ([NaN 0 1 2 NaN NaN]))
%!assert (hygeinv ([x, NaN], 4, 2, single (2)), single ([NaN 0 1 2 NaN NaN]))

## Test input validation
%!error hygeinv ()
%!error hygeinv (1)
%!error hygeinv (1,2)
%!error hygeinv (1,2,3)
%!error hygeinv (1,2,3,4,5)
%!error hygeinv (ones (2), ones (3), 1, 1)
%!error hygeinv (1, ones (2), ones (3), 1)
%!error hygeinv (1, 1, ones (2), ones (3))
%!error hygeinv (i, 2, 2, 2)
%!error hygeinv (2, i, 2, 2)
%!error hygeinv (2, 2, i, 2)
%!error hygeinv (2, 2, 2, i)
