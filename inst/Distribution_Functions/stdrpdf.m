## Copyright (C) 2026 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
## more details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn {statistics} {@var{y} =} stdrpdf (@var{x}, @var{k}, @var{df})
##
## Studentized range probability density function (PDF).
##
## For each element of @var{x}, compute the probability density function (PDF)
## of the studentized range distribution for @var{k} groups and @var{df}
## degrees of freedom.  The size of @var{y} is the common size of @var{x},
## @var{k} and @var{df}.  A scalar input functions as a constant matrix of the
## same size as the other inputs.
##
## @var{k} must be an integer of at least 2 and @var{df} positive, @code{Inf}
## included; otherwise @var{y} is @code{NaN}.  The density is zero for
## negative @var{x} and, for @var{k} greater than 2, at zero.
##
## MATLAB has no public counterpart of this function.
##
## Further information about the studentized range distribution can be found
## at @url{https://en.wikipedia.org/wiki/Studentized_range_distribution}
##
## Input arguments must be @qcode{double} or @qcode{single}; integer, logical,
## and character arrays are rejected.
##
## @seealso{stdrcdf, stdrinv, stdrrnd, stdrstat}
## @end deftypefn

function y = stdrpdf (x, k, df)

  ## Check for valid number of input arguments
  if (nargin < 3)
    error ("stdrpdf: function called with too few input arguments.");
  endif

  ## Check for common size of X, K, and DF
  if (! isscalar (x) || ! isscalar (k) || ! isscalar (df))
    [err, x, k, df] = common_size (x, k, df);
    if (err > 0)
      error ("stdrpdf: X, K, and DF must be of common size or scalars.");
    endif
  endif

  ## Check for X, K, and DF being double or single
  if (! (isfloat (x) && isfloat (k) && isfloat (df)))
    error ("stdrpdf: X, K, and DF must be double or single.");
  endif

  ## Check for X, K, and DF being reals
  if (iscomplex (x) || iscomplex (k) || iscomplex (df))
    error ("stdrpdf: X, K, and DF must not be complex.");
  endif

  ## Check for class type
  if (isa (x, 'single') || isa (k, 'single') || isa (df, 'single'))
    cls = 'single';
  else
    cls = 'double';
  endif
  sz = size (x);
  x = double (x(:));
  k = double (k(:)) .* ones (numel (x), 1);
  df = double (df(:)) .* ones (numel (x), 1);

  ## Invalid parameters or X give NaN
  valid = (k >= 2) & (k == fix (k)) & isfinite (k) & (df > 0);
  y = zeros (numel (x), 1);
  y(! valid | isnan (x)) = NaN;

  ## At zero only two groups have a density: their range is sqrt (2) times
  ## the absolute value of a t variable
  z = valid & (x == 0) & (k == 2);
  y(z) = sqrt (2) * tpdf (0, df(z));

  ## Evaluate the positive finite values, one parameter pair at a time
  todo = valid & (x > 0) & isfinite (x);
  if (any (todo))
    pairs = unique ([k(todo), df(todo)], 'rows');
    for i = 1:rows (pairs)
      idx = todo & (k == pairs(i,1)) & (df == pairs(i,2));
      y(idx) = __stdr__ (x(idx), pairs(i,1), pairs(i,2), 'pdf');
    endfor
  endif
  y = cast (reshape (y, sz), cls);

endfunction

%!demo
%! ## Plot various PDFs from the studentized range distribution
%! x = 0:0.01:8;
%! y1 = stdrpdf (x, 2, 5);
%! y2 = stdrpdf (x, 3, 5);
%! y3 = stdrpdf (x, 5, 5);
%! y4 = stdrpdf (x, 5, Inf);
%! plot (x, y1, '-b', x, y2, '-g', x, y3, '-r', x, y4, '-m')
%! grid on
%! legend ({'k = 2, df = 5', 'k = 3, df = 5', 'k = 5, df = 5', ...
%!          'k = 5, df = \infty'}, 'location', 'northeast')
%! title ('Studentized range PDF')
%! xlabel ('values in x')
%! ylabel ('density')

## Two groups are sqrt (2) times the absolute value of a t variable
%!shared x
%! x = [0.5, 2, 4, 8];
%!assert_equal (stdrpdf (x, 2, 5), sqrt (2) * tpdf (x / sqrt (2), 5), -1e-13)
%!assert_equal (stdrpdf (x, 2, Inf), exp (-x .^ 2 / 4) / sqrt (pi), -1e-13)
%!assert_equal (stdrpdf (0, 2, 5), sqrt (2) * tpdf (0, 5), -1e-13)
%!assert_equal (stdrpdf (0, 3, 5), 0)

## The density is the derivative of the distribution function
%!test
%! h = 1e-4;
%! d = (stdrcdf (3 + h, 5, 10) - stdrcdf (3 - h, 5, 10)) / (2 * h);
%! assert_equal (stdrpdf (3, 5, 10), d, -1e-7);

## Edge values and invalid parameters
%!assert_equal (stdrpdf ([-1, Inf, NaN], 3, 10), [0, 0, NaN])
%!assert_equal (stdrpdf (2, [1, 2.5, NaN, Inf], 10), NaN (1, 4))
%!assert_equal (stdrpdf (2, 3, [0, -1, NaN]), NaN (1, 3))

## Test class of input preserved
%!assert_equal (class (stdrpdf (single (2), 3, 10)), 'single')
%!assert_equal (class (stdrpdf (2, 3, single (10))), 'single')

## Test input validation
%!error<stdrpdf: function called with too few input arguments.> stdrpdf ()
%!error<stdrpdf: function called with too few input arguments.> stdrpdf (1, 2)
%!error<stdrpdf: X, K, and DF must be of common size or scalars.> ...
%! stdrpdf (ones (3), ones (2), 3)
%!error<stdrpdf: X, K, and DF must be double or single.> ...
%! stdrpdf (int32 (2), 3, 10)
%!error<stdrpdf: X, K, and DF must be double or single.> ...
%! stdrpdf (true, 3, 10)
%!error<stdrpdf: X, K, and DF must be double or single.> ...
%! stdrpdf ('a', 3, 10)
%!error<stdrpdf: X, K, and DF must not be complex.> stdrpdf (i, 3, 10)
%!error<stdrpdf: X, K, and DF must not be complex.> stdrpdf (2, i, 10)
%!error<stdrpdf: X, K, and DF must not be complex.> stdrpdf (2, 3, i)
