## Copyright (C) 2016-2017 Lachlan Andrew
## Copyright (C) 2012-2016 Rik Wehbring
## Copyright (C) 1995-2012 Kurt Hornik
## Copyright (C) 2023 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {statistics} {@var{x} =} binoinv (@var{p}, @var{n}, @var{ps})
##
## Inverse of the Binomial cumulative distribution function (iCDF).
##
## For each element of @var{p}, compute the quantile (the inverse of the CDF) of
## the binomial distribution with parameters @var{n} and @var{ps}, where @var{n}
## is the number of trials and @var{ps} is the probability of success.  The size
## of @var{x} is the common size of @var{p}, @var{n}, and @var{ps}.  A scalar
## input functions as a constant matrix of the same size as the other inputs.
##
## Further information about the binomial distribution can be found at
## @url{https://en.wikipedia.org/wiki/Binomial_distribution}
##
## Input arguments must be @qcode{double} or @qcode{single}; integer, logical,
## and character arrays are rejected.  MATLAB accepts a character array and
## evaluates it at the character codes, which Octave deliberately does not,
## since a character array is an integer type and integers are refused too.
##
## @seealso{binocdf, binopdf, binornd, binofit, binolike, binostat, binotest}
## @end deftypefn

function x = binoinv (p, n, ps)

  ## Check for valid number of input arguments
  if (nargin < 3)
    error ("binoinv: function called with too few input arguments.");
  endif

  ## Check for common size of P, N, and PS
  if (! isscalar (n) || ! isscalar (ps))
    [retval, p, n, ps] = common_size (p, n, ps);
    if (retval > 0)
      error ("binoinv: P, N, and PS must be of common size or scalars.");
    endif
  endif

  ## Check for P, N, and PS being double or single
  if (! (isfloat (p) && isfloat (n) && isfloat (ps)))
    error ("binoinv: P, N, and PS must be double or single.");
  endif

  ## Check for P, N, and PS being reals
  if (iscomplex (p) || iscomplex (n) || iscomplex (ps))
    error ("binoinv: P, N, and PS must not be complex.");
  endif

  ## Check for class type
  if (isa (p, 'single') || isa (n, 'single') || isa (ps, 'single'));
    x = zeros (size (p), 'single');
  else
    x = zeros (size (p));
  endif

  k = (! (p >= 0) | ! (p <= 1) | ! (n >= 0) | (n != fix (n)) | ! (ps >= 0) | ...
       ! (ps <= 1));
  x(k) = NaN;

  k = find ((p >= 0) & (p <= 1) & (n >= 0) & (n == fix (n) ...
                     & (ps >= 0) & (ps <= 1)));
  if (! isempty (k))
    pk = p(k)(:);
    if (isscalar (n) && isscalar (ps))
      [xk, unfinished] = scalar_binoinv (pk, n, ps);
      if (! isempty (unfinished))
        xk(unfinished) = bin_search_binoinv (pk(unfinished), n, ps);
      endif
      xk = tie_correct (xk, pk, n, ps);
      ## P of 1 attains N.  The CDF saturates to 1 well below it, reading 1 at
      ## 206 for 500 trials at 0.25, so no search can locate it.
      xk(pk == 1) = n;
    else
      nk = n(k)(:);
      psk = ps(k)(:);
      [xk, unfinished] = vector_binoinv (pk, nk, psk);
      if (! isempty (unfinished))
        xk(unfinished) = bin_search_binoinv (pk(unfinished), nk(unfinished), ...
                                             psk(unfinished));
      endif
      xk = tie_correct (xk, pk, nk, psk);
      xk(pk == 1) = nk(pk == 1);
    endif
    x(k) = xk;
  endif

endfunction

## Core algorithm to calculate the inverse binomial, for n and ps real scalars
## and x a column vector, and for which the output is not NaN or Inf.
## Compute CDF in batches of doubling size until CDF > p, or answer > 500
## Return the locations of unfinished cases in k.
function [m, k] = scalar_binoinv (p, n, ps)

  k = 1:length (p);
  m = zeros (size (p));
  prev_limit = 0;
  limit = 10;
  cdf = 0;
  v = 0;
  do
    cdf = binocdf (prev_limit:limit-1, n, ps);
    r = bsxfun (@le, p(k), cdf);
    [v, m(k)] = max (r, [], 2);     # find first instance of p <= cdf
    m(k) += prev_limit - 1;
    k = k(v == 0);

    prev_limit = limit;
    limit += limit;
  until (isempty (k) || limit >= 1000)

endfunction

## Core algorithm to calculate the inverse binomial, for n, ps, and x column
## vectors, and for which the output is not NaN or Inf.
## Compute CDF in batches of doubling size until CDF > p, or answer > 500
## Return the locations of unfinished cases in k.
## Calculates CDF by summing PDF, which is faster than calls to binocdf.
function [m, k] = vector_binoinv (p, n, ps)

  k = 1:length (p);
  m = zeros (size (p));
  prev_limit = 0;
  limit = 10;
  cdf = 0;
  v = 0;
  do
    xx = repmat (prev_limit:limit-1, [length(k), 1]);
    nn = kron (ones (1, limit-prev_limit), n(k));
    pp = kron (ones (1, limit-prev_limit), ps(k));
    pdf = binopdf (xx, nn, pp);
    pdf(:,1) += cdf(v==0, end);
    cdf = cumsum (pdf, 2);
    r = bsxfun (@le, p(k), cdf);
    [v, m(k)] = max (r, [], 2);     # find first instance of p <= cdf
    m(k) += prev_limit - 1;
    k = k(v == 0);

    prev_limit = limit;
    limit += min (limit, max (1e4/numel (k), 10));  # limit memory use
  until (isempty (k) || limit >= 1000)

endfunction

## Vectorized binary search.
## Can handle vectors n and ps, and is faster than the scalar case when the
## answer is large.
## Could be optimized to call binocdf only for a subset of the p at each stage,
## but care must be taken to handle both scalar and vector n, ps.  Bookkeeping
## may cost more than the extra computations.
function m = bin_search_binoinv (p, n, ps)

  ## binocdf returns a column, so a row P would broadcast into a full matrix.
  p = p(:);
  k = 1:length (p);
  lower = zeros (size (p));
  limit = 500;              # lower bound on point at which prev phase finished
  while (any (k) && limit < 1e100)
    cdf = binocdf (limit, n, ps);
    k = (p > cdf);
    lower(k) = limit;
    limit += limit;
  endwhile
  upper = max (2*lower, 1);
  k = find (lower != limit/2);       # elements for which above loop finished
  for i = 1:ceil (log2 (max (lower)))
    mid = (upper + lower)/2;
    cdf = binocdf (floor (mid(:)), n, ps);
    r = (p <= cdf);
    upper(r)  = mid(r);
    lower(! r) = mid(! r);
  endfor
  m = ceil (lower);
  m(p > binocdf (m(:), n, ps)) += 1;  # fix off-by-one errors from binary search

endfunction

## Step the answer back onto M-1 wherever P lies within the error of the CDF
## there.  The CDF carries a relative error that grows with N, so an exactly
## attained probability reads short and sends the search past its answer:
## binocdf (2, 5, 0.5) is 0.49999999999999989 where the exact value, 16/32,
## is representable, and the shortfall reaches 2355 eps by 2001 trials.  The
## slack is capped at half the probability of M so that it can never cross a
## real step of the distribution and round P to the wrong side of one.
function m = tie_correct (m, p, n, ps)

  j = find (m > 0);
  if (isempty (j))
    return;
  endif
  mj = m(j);
  if (isscalar (n))
    nj = n;
  else
    nj = n(j);
  endif
  if (isscalar (ps))
    psj = ps;
  else
    psj = ps(j);
  endif
  lo = binocdf (mj - 1, nj, psj);
  tol = min (32 * nj .* eps (lo), 0.5 * binopdf (mj, nj, psj));
  m(j(p(j) <= lo + tol)) -= 1;

endfunction

%!demo
%! ## Plot various iCDFs from the binomial distribution
%! p = 0.001:0.001:0.999;
%! x1 = binoinv (p, 20, 0.5);
%! x2 = binoinv (p, 20, 0.7);
%! x3 = binoinv (p, 40, 0.5);
%! plot (p, x1, '-b', p, x2, '-g', p, x3, '-r')
%! grid on
%! legend ({'n = 20, ps = 0.5', 'n = 20, ps = 0.7', ...
%!          'n = 40, ps = 0.5'}, 'location', 'southeast')
%! title ('Binomial iCDF')
%! xlabel ('probability')
%! ylabel ('values in x (number of successes)')

## Test output
%!shared p
%! p = [-1 0 0.5 1 2];
%!assert_equal (binoinv (p, 2*ones (1,5), 0.5*ones (1,5)), [NaN 0 1 2 NaN])
%!assert_equal (binoinv (p, 2, 0.5*ones (1,5)), [NaN 0 1 2 NaN])
%!assert_equal (binoinv (p, 2*ones (1,5), 0.5), [NaN 0 1 2 NaN])
%!assert_equal (binoinv (p, 2*[0 -1 NaN 1.1 1], 0.5), [NaN NaN NaN NaN NaN])
%!assert_equal (binoinv (p, 2, 0.5*[0 -1 NaN 3 1]), [NaN NaN NaN NaN NaN])
%!assert_equal (binoinv ([p(1:2) NaN p(4:5)], 2, 0.5), [NaN 0 NaN 2 NaN])

## Test class of input preserved
%!assert_equal (binoinv ([p, NaN], 2, 0.5), [NaN 0 1 2 NaN NaN])
%!assert_equal (binoinv (single ([p, NaN]), 2, 0.5), single ([NaN 0 1 2 NaN NaN]))
%!assert_equal (binoinv ([p, NaN], single (2), 0.5), single ([NaN 0 1 2 NaN NaN]))
%!assert_equal (binoinv ([p, NaN], 2, single (0.5)), single ([NaN 0 1 2 NaN NaN]))

## Test accuracy against a round trip through the CDF
%!shared x
%! x = magic (3) + 1;
%!assert_equal (binoinv (binocdf (1:10, 11, 0.1), 11, 0.1), 1:10)
%!assert_equal (binoinv (binocdf (1:10, 2*(1:10), 0.1), 2*(1:10), 0.1), 1:10)
%!assert_equal (binoinv (binocdf (x, 2*x, 1./x), 2*x, 1./x), x)

## A symmetric binomial with an odd number of trials attains 0.5 exactly at
## its median, which the CDF reads short of by an error growing with N.
%!assert_equal (binoinv (0.5, 5, 0.5), 2)
%!assert_equal (binoinv (0.5, 7, 0.5), 3)
%!assert_equal (binoinv (0.5, 501, 0.5), 250)
%!assert_equal (binoinv (0.5, 2001, 0.5), 1000)
%!assert_equal (binoinv (0.5, 100001, 0.5), 50000)

## P of 1 attains N even where the CDF saturates to 1 far below it.
%!assert_equal (binoinv (1, 500, 0.25), 500)
%!assert_equal (binoinv (1, [50, 500], [0.5, 0.1]), [50, 500])

## Answers above 500 take the binary search, where the CDF is a column and P
## keeps the orientation it was given.
%!assert_equal (binoinv ([0.5, 0.9], 2001, 0.5), [1000, 1029])
%!assert_equal (binoinv ([0.5; 0.9], 2001, 0.5), [1000; 1029])
%!assert_equal (binoinv ([NaN, 0.5], 2001, 0.5), [NaN, 1000])
%!assert_equal (binoinv ([2, 0.5], [2001, 2001], [0.5, 0.5]), [NaN, 1000])

## Test input validation
%!error<binoinv: function called with too few input arguments.> binoinv ()
%!error<binoinv: function called with too few input arguments.> binoinv (1)
%!error<binoinv: function called with too few input arguments.> binoinv (1,2)
%!error<binoinv: function called with too many inputs> binoinv (1,2,3,4)
%!error<binoinv: P, N, and PS must be of common size or scalars.> ...
%! binoinv (ones (3), ones (2), ones (2))
%!error<binoinv: P, N, and PS must be of common size or scalars.> ...
%! binoinv (ones (2), ones (3), ones (2))
%!error<binoinv: P, N, and PS must be of common size or scalars.> ...
%! binoinv (ones (2), ones (2), ones (3))
%!error<binoinv: P, N, and PS must be double or single.> binoinv (int32 (2), 2, 2)
%!error<binoinv: P, N, and PS must be double or single.> binoinv (true, 2, 2)
%!error<binoinv: P, N, and PS must be double or single.> binoinv ('a', 2, 2)
%!error<binoinv: P, N, and PS must not be complex.> binoinv (i, 2, 2)
%!error<binoinv: P, N, and PS must not be complex.> binoinv (2, i, 2)
%!error<binoinv: P, N, and PS must not be complex.> binoinv (2, 2, i)
