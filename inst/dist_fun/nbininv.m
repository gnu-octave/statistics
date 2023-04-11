## Copyright (C) 1995-2012 Kurt Hornik
## Copyright (C) 2012-2016 Rik Wehbring
## Copyright (C) 2016-2017 Lachlan Andrew
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
## @deftypefn  {statistics} {@var{x} =} nbininv (@var{p}, @var{n}, @var{ps})
##
## Inverse of the negative binomial cumulative distribution function (iCDF).
##
## For each element of @var{p}, compute the quantile (the inverse of the CDF)
## at @var{p} of the negative binomial distribution with parameters @var{n} and
## @var{ps}.  The size of @var{x} is the common size of @var{p}, @var{n}, and
## @var{ps}.  A scalar input functions as a constant matrix of the same size as
## the other inputs.
##
## When @var{n} is integer this is the Pascal distribution.
## When @var{n} is extended to real numbers this is the Polya distribution.
##
## The number of failures in a Bernoulli experiment with success probability
## @var{ps} before the @var{n}-th success follows this distribution.
##
## @seealso{nbininv, nbinpdf, nbinrnd, nbinstat}
## @end deftypefn

function x = nbininv (p, n, ps)

  ## Check for valid number of input arguments
  if (nargin != 3)
    print_usage ();
  endif

  ## Check for common size of P, N, and PS
  if (! isscalar (p) || ! isscalar (n) || ! isscalar (ps))
    [retval, p, n, ps] = common_size (p, n, ps);
    if (retval > 0)
      error ("nbininv: P, N, and PS must be of common size or scalars.");
    endif
  endif

  ## Check for P, N, and PS being reals
  if (iscomplex (p) || iscomplex (n) || iscomplex (ps))
    error ("nbininv: P, N, and PS must not be complex.");
  endif

  ## Check for appropriate class
  if (isa (p, "single") || isa (n, "single") || isa (ps, "single"))
    x = zeros (size (p), "single");
  else
    x = zeros (size (p));
  endif

  k = (isnan (p) | (p < 0) | (p > 1) | isnan (n) | (n < 1) | (n == Inf)
       | isnan (ps) | (ps < 0) | (ps > 1));
  x(k) = NaN;

  k = (p == 1) & (n > 0) & (n < Inf) & (ps >= 0) & (ps <= 1);
  x(k) = Inf;

  k = find ((p >= 0) & (p < 1) & (n > 0) & (n < Inf)
            & (ps > 0) & (ps <= 1));
  if (! isempty (k))
    p = p(k);
    m = zeros (size (k));
    if (isscalar (n) && isscalar (ps))
      [m, unfinished] = scalar_nbininv (p(:), n, ps);
      m(unfinished) = bin_search_nbininv (p(unfinished), n, ps);
    else
      m = bin_search_nbininv (p, n(k), ps(k));
    endif
    x(k) = m;
  endif

endfunction


## Core algorithm to calculate the inverse negative binomial, for n and ps real
## scalars and y a column vector, and for which the output is not NaN or Inf.
## Compute CDF in batches of doubling size until CDF > p, or answer > 500.
## Return the locations of unfinished cases in k.
function [m, k] = scalar_nbininv (p, n, ps)

  k = 1:length (p);
  m = zeros (size (p));
  prev_limit = 0;
  limit = 10;
  do
    cdf = nbincdf (prev_limit:limit, n, ps);
    r = bsxfun (@le, p(k), cdf);
    [v, m(k)] = max (r, [], 2);     # find first instance of p <= cdf
    m(k) += prev_limit - 1;
    k = k(v == 0);

    prev_limit = limit;
    limit += limit;
  until (isempty (k) || limit >= 1000)

endfunction

## Vectorized binary search.
## Can handle vectors n and ps, and is faster than the scalar case when the
## answer is large.
## Could be optimized to call nbincdf only for a subset of the p at each stage,
## but care must be taken to handle both scalar and vector n,ps.  Bookkeeping
## may cost more than the extra computations.
function m = bin_search_nbininv (p, n, ps)

  k = 1:length (p);
  lower = zeros (size (p));
  limit = 1;
  while (any (k) && limit < 1e100)
    cdf = nbincdf (limit, n, ps);
    k = (p > cdf);
    lower(k) = limit;
    limit += limit;
  endwhile
  upper = max (2*lower, 1);
  k = find (lower != limit/2);    # elements for which above loop finished
  for i = 1:ceil (log2 (max (lower)))
    mid = (upper + lower)/2;
    cdf = nbincdf (floor (mid), n, ps);
    r = (p <= cdf);
    upper(r)  = mid(r);
    lower(! r) = mid(! r);
  endfor
  m = ceil (lower);
  m(p > nbincdf (m, n, ps)) += 1;  # fix off-by-one errors from binary search

endfunction


%!shared p
%! p = [-1 0 3/4 1 2];
%!assert (nbininv (p, ones (1,5), 0.5*ones (1,5)), [NaN 0 1 Inf NaN])
%!assert (nbininv (p, 1, 0.5*ones (1,5)), [NaN 0 1 Inf NaN])
%!assert (nbininv (p, ones (1,5), 0.5), [NaN 0 1 Inf NaN])
%!assert (nbininv (p, [1 0 NaN Inf 1], 0.5), [NaN NaN NaN NaN NaN])
%!assert (nbininv (p, [1 0 1.5 Inf 1], 0.5), [NaN NaN 2 NaN NaN])
%!assert (nbininv (p, 1, 0.5*[1 -Inf NaN Inf 1]), [NaN NaN NaN NaN NaN])
%!assert (nbininv ([p(1:2) NaN p(4:5)], 1, 0.5), [NaN 0 NaN Inf NaN])

## Test class of input preserved
%!assert (nbininv ([p, NaN], 1, 0.5), [NaN 0 1 Inf NaN NaN])
%!assert (nbininv (single ([p, NaN]), 1, 0.5), single ([NaN 0 1 Inf NaN NaN]))
%!assert (nbininv ([p, NaN], single (1), 0.5), single ([NaN 0 1 Inf NaN NaN]))
%!assert (nbininv ([p, NaN], 1, single (0.5)), single ([NaN 0 1 Inf NaN NaN]))

## Test accuracy, to within +/- 1 since it is a discrete distribution
%!shared y, tol
%! y = magic (3) + 1;
%! tol = 1;
%!assert (nbininv (nbincdf (1:10, 3, 0.1), 3, 0.1), 1:10, tol)
%!assert (nbininv (nbincdf (1:10, 3./(1:10), 0.1), 3./(1:10), 0.1), 1:10, tol)
%!assert (nbininv (nbincdf (y, 3./y, 1./y), 3./y, 1./y), y, tol)

## Test input validation
%!error nbininv ()
%!error nbininv (1)
%!error nbininv (1,2)
%!error nbininv (1,2,3,4)
%!error nbininv (ones (3), ones (2), ones (2))
%!error nbininv (ones (2), ones (3), ones (2))
%!error nbininv (ones (2), ones (2), ones (3))
%!error nbininv (i, 2, 2)
%!error nbininv (2, i, 2)
%!error nbininv (2, 2, i)
