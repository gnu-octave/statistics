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
## @deftypefn  {statistics} {@var{x} =} nbininv (@var{p}, @var{r}, @var{ps})
##
## Inverse of the negative binomial cumulative distribution function (iCDF).
##
## For each element of @var{p}, compute the quantile (the inverse of the CDF) of
## the negative binomial distribution with parameters @var{r} and @var{ps},
## where @var{r} is the number of successes until the experiment is stopped and
## @var{ps} is the probability of success in each experiment, given the
## probability in @var{p}.  The size of @var{x} is the common size of @var{p},
## @var{r}, and @var{ps}.  A scalar input functions as a constant matrix of the
## same size as the other inputs.
##
## When @var{r} is an integer, the negative binomial distribution is also known
## as the Pascal distribution and it models the number of failures in @var{x}
## before a specified number of successes is reached in a series of independent,
## identical trials.  Its parameters are the probability of success in a single
## trial, @var{ps}, and the number of successes, @var{r}.  A special case of the
## negative binomial distribution, when @qcode{@var{r} = 1}, is the geometric
## distribution, which models the number of failures before the first success.
##
## @var{r} can also have non-integer positive values, in which form the negative
## binomial distribution, also known as the Polya distribution, has no
## interpretation in terms of repeated trials, but, like the Poisson
## distribution, it is useful in modeling count data.  The negative binomial
## distribution is more general than the Poisson distribution because it has a
## variance that is greater than its mean, making it suitable for count data
## that do not meet the assumptions of the Poisson distribution.  In the limit,
## as @var{r} increases to infinity, the negative binomial distribution
## approaches the Poisson distribution.
##
## Further information about the negative binomial distribution can be found at
## @url{https://en.wikipedia.org/wiki/Negative_binomial_distribution}
##
## @seealso{nbininv, nbinpdf, nbinrnd, nbinfit, nbinlike, nbinstat}
## @end deftypefn

function x = nbininv (p, r, ps)

  ## Check for valid number of input arguments
  if (nargin < 3)
    error ("nbininv: function called with too few input arguments.");
  endif

  ## Check for common size of P, R, and PS
  if (! isscalar (p) || ! isscalar (r) || ! isscalar (ps))
    [retval, p, r, ps] = common_size (p, r, ps);
    if (retval > 0)
      error ("nbininv: P, R, and PS must be of common size or scalars.");
    endif
  endif

  ## Check for P, R, and PS being reals
  if (iscomplex (p) || iscomplex (r) || iscomplex (ps))
    error ("nbininv: P, R, and PS must not be complex.");
  endif

  ## Check for class type
  if (isa (p, "single") || isa (r, "single") || isa (ps, "single"))
    x = zeros (size (p), "single");
  else
    x = zeros (size (p));
  endif

  k = (isnan (p) | (p < 0) | (p > 1) | isnan (r) | (r < 1) | (r == Inf)
       | isnan (ps) | (ps < 0) | (ps > 1));
  x(k) = NaN;

  k = (p == 1) & (r > 0) & (r < Inf) & (ps >= 0) & (ps <= 1);
  x(k) = Inf;

  k = find ((p >= 0) & (p < 1) & (r > 0) & (r < Inf)
            & (ps > 0) & (ps <= 1));
  if (! isempty (k))
    p = p(k);
    m = zeros (size (k));
    if (isscalar (r) && isscalar (ps))
      [m, unfinished] = scalar_nbininv (p(:), r, ps);
      m(unfinished) = bin_search_nbininv (p(unfinished), r, ps);
    else
      m = bin_search_nbininv (p, r(k), ps(k));
    endif
    x(k) = m;
  endif

endfunction


## Core algorithm to calculate the inverse negative binomial, for r and ps real
## scalars and y a column vector, and for which the output is not NaN or Inf.
## Compute CDF in batches of doubling size until CDF > p, or answer > 500.
## Return the locations of unfinished cases in k.
function [m, k] = scalar_nbininv (p, r, ps)

  k = 1:length (p);
  m = zeros (size (p));
  prev_limit = 0;
  limit = 10;
  do
    cdf = nbincdf (prev_limit:limit, r, ps);
    rr = bsxfun (@le, p(k), cdf);
    [v, m(k)] = max (rr, [], 2);     # find first instance of p <= cdf
    m(k) += prev_limit - 1;
    k = k(v == 0);

    prev_limit = limit;
    limit += limit;
  until (isempty (k) || limit >= 1000)

endfunction

## Vectorized binary search.
## Can handle vectors r and ps, and is faster than the scalar case when the
## answer is large.
## Could be optimized to call nbincdf only for a subset of the p at each stage,
## but care must be taken to handle both scalar and vector r,ps.  Bookkeeping
## may cost more than the extra computations.
function m = bin_search_nbininv (p, r, ps)

  k = 1:length (p);
  lower = zeros (size (p));
  limit = 1;
  while (any (k) && limit < 1e100)
    cdf = nbincdf (limit, r, ps);
    k = (p > cdf);
    lower(k) = limit;
    limit += limit;
  endwhile
  upper = max (2*lower, 1);
  k = find (lower != limit/2);    # elements for which above loop finished
  for i = 1:ceil (log2 (max (lower)))
    mid = (upper + lower)/2;
    cdf = nbincdf (floor (mid), r, ps);
    rr = (p <= cdf);
    upper(rr)  = mid(rr);
    lower(! rr) = mid(! rr);
  endfor
  m = ceil (lower);
  m(p > nbincdf (m, r, ps)) += 1;  # fix off-by-one errors from binary search

endfunction

%!demo
%! ## Plot various iCDFs from the negative binomial distribution
%! p = 0.001:0.001:0.999;
%! x1 = nbininv (p, 2, 0.15);
%! x2 = nbininv (p, 5, 0.2);
%! x3 = nbininv (p, 4, 0.4);
%! x4 = nbininv (p, 10, 0.3);
%! plot (p, x1, "-r", p, x2, "-g", p, x3, "-k", p, x4, "-m")
%! grid on
%! ylim ([0, 40])
%! legend ({"r = 2, ps = 0.15", "r = 5, ps = 0.2", "r = 4, p = 0.4", ...
%!          "r = 10, ps = 0.3"}, "location", "northwest")
%! title ("Negative binomial iCDF")
%! xlabel ("probability")
%! ylabel ("values in x (number of failures)")

## Test output
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
%!error<nbininv: function called with too few input arguments.> nbininv ()
%!error<nbininv: function called with too few input arguments.> nbininv (1)
%!error<nbininv: function called with too few input arguments.> nbininv (1, 2)
%!error<nbininv: P, R, and PS must be of common size or scalars.> ...
%! nbininv (ones (3), ones (2), ones (2))
%!error<nbininv: P, R, and PS must be of common size or scalars.> ...
%! nbininv (ones (2), ones (3), ones (2))
%!error<nbininv: P, R, and PS must be of common size or scalars.> ...
%! nbininv (ones (2), ones (2), ones (3))
%!error<nbininv: P, R, and PS must not be complex.> nbininv (i, 2, 2)
%!error<nbininv: P, R, and PS must not be complex.> nbininv (2, i, 2)
%!error<nbininv: P, R, and PS must not be complex.> nbininv (2, 2, i)
