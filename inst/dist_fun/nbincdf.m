## Copyright (C) 2012 Rik Wehbring
## Copyright (C) 1995-2016 Kurt Hornik
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
## @deftypefn  {statistics} {@var{p} =} nbincdf (@var{x}, @var{r}, @var{ps})
## @deftypefnx {statistics} {@var{p} =} nbincdf (@var{x}, @var{r}, @var{ps}, @qcode{"upper"})
##
## Negative binomial cumulative distribution function (CDF).
##
## For each element of @var{x}, compute the cumulative distribution function
## (CDF) of the negative binomial distribution with parameters @var{r} and
## @var{ps}, where @var{r} is the number of successes until the experiment is
## stopped and @var{ps} is the probability of success in each experiment, given
## the number of failures in @var{x}.  The size of @var{p} is the common size of
## @var{x}, @var{r}, and @var{ps}.  A scalar input functions as a constant
## matrix of the same size as the other inputs.
##
## The algorithm uses the cumulative sums of the binomial masses.
##
## @code{@var{p} = nbincdf (@var{x}, @var{r}, @var{ps}, "upper")} computes the
## upper tail probability of the negative binomial distribution with parameters
## @var{r} and @var{ps}, at the values in @var{x}.
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

function p = nbincdf (x, r, ps, uflag)

  ## Check for valid number of input arguments
  if (nargin < 3)
    error ("nbincdf: function called with too few input arguments.");
  endif

  ## Check for "upper" flag
  if (nargin == 4 && strcmpi (uflag, "upper"))
    uflag = true;
  elseif (nargin == 4  && ! strcmpi (uflag, "upper"))
    error ("nbincdf: invalid argument for upper tail.");
  else
    uflag = false;
  endif

  ## Check for R and PS being scalars
  scalarNPS = (isscalar(r) & isscalar(ps));

  ## Check for common size of X, R, and PS
  if (! isscalar (x) || ! isscalar (r) || ! isscalar (ps))
    [retval, x, r, ps] = common_size (x, r, ps);
    if (retval > 0)
      error ("nbincdf: X, R, and PS must be of common size or scalars.");
    endif
  endif

  ## Check for X, R, and PS being reals
  if (iscomplex (x) || iscomplex (r) || iscomplex (ps))
    error ("nbincdf: X, R, and PS must not be complex.");
  endif

  ## Check for class type
  if (isa (x, "single") || isa (r, "single") || isa (ps, "single"))
    p = zeros (size (x), "single");
  else
    p = zeros (size (x));
  endif

  ## Force NaN for out of range or missing parameters and missing data NaN
  is_nan = (isnan (x) | isnan (r) | (r <= 0) | (r == Inf) | (ps < 0) | (ps > 1));
  p(is_nan) = NaN;

  ## Compute P for X >= 0
  xf = floor (x);
  k = find (xf >= 0 & ! is_nan);

  ## Return 1 for positive infinite values of X, unless "upper" is given: p = 0
  k1 = find (isinf (xf(k)));
  if (any (k1))
    if (uflag)
      p(k(k1)) = 0;
    else
      p(k(k1)) = 1;
    endif
    k(k1) = [];
  endif

  ## Return 1 when X < 0 and "upper" is given
  k1 = (x < 0 & ! is_nan);
  if (any (k1))
    if (uflag)
      p(k1) = 1;
    endif
  endif

  ## Accumulate probabilities up to the maximum value in X
  if (any (k))
    if (uflag)
      p(k) = betainc (ps(k), r(k), xf(k) + 1, "upper");
    else
      max_val = max (xf(k));
      if (scalarNPS)
        tmp = cumsum (nbinpdf (0:max_val, r(1), ps(1)));
        p(k) = tmp(xf(k) + 1);
      else
        idx = (0:max_val)';
        compare = idx(:, ones (size (k)));
        index = xf(k);
        index = index(:);
        index = index(:, ones (size (idx)))';
        n_big = r(k);
        n_big = n_big(:);
        n_big = n_big(:, ones (size (idx)))';
        ps_big = ps(k);
        ps_big = ps_big(:);
        ps_big = ps_big(:, ones (size (idx)))';
        p0 = nbinpdf (compare, n_big, ps_big);
        indicator = find (compare > index);
        p0(indicator) = zeros (size (indicator));
        p(k) = sum(p0,1);
      endif
    endif
  endif

  ## Prevent round-off errors
  p(p > 1) = 1;

endfunction

%!demo
%! ## Plot various CDFs from the negative binomial distribution
%! x = 0:50;
%! p1 = nbincdf (x, 2, 0.15);
%! p2 = nbincdf (x, 5, 0.2);
%! p3 = nbincdf (x, 4, 0.4);
%! p4 = nbincdf (x, 10, 0.3);
%! plot (x, p1, "*r", x, p2, "*g", x, p3, "*k", x, p4, "*m")
%! grid on
%! xlim ([0, 40])
%! legend ({"r = 2, ps = 0.15", "r = 5, ps = 0.2", "r = 4, p = 0.4", ...
%!          "r = 10, ps = 0.3"}, "location", "southeast")
%! title ("Negative binomial CDF")
%! xlabel ("values in x (number of failures)")
%! ylabel ("probability")

## Test output
%!shared x, y
%! x = [-1 0 1 2 Inf];
%! y = [0 1/2 3/4 7/8 1];
%!assert (nbincdf (x, ones (1,5), 0.5*ones (1,5)), y)
%!assert (nbincdf (x, 1, 0.5*ones (1,5)), y)
%!assert (nbincdf (x, ones (1,5), 0.5), y)
%!assert (nbincdf (x, ones (1,5), 0.5, "upper"), 1 - y, eps)
%!assert (nbincdf ([x(1:3) 0 x(5)], [0 1 NaN 1.5 Inf], 0.5), ...
%! [NaN 1/2 NaN nbinpdf(0,1.5,0.5) NaN], eps)
%!assert (nbincdf (x, 1, 0.5*[-1 NaN 4 1 1]), [NaN NaN NaN y(4:5)])
%!assert (nbincdf ([x(1:2) NaN x(4:5)], 1, 0.5), [y(1:2) NaN y(4:5)])

## Test class of input preserved
%!assert (nbincdf ([x, NaN], 1, 0.5), [y, NaN])
%!assert (nbincdf (single ([x, NaN]), 1, 0.5), single ([y, NaN]))
%!assert (nbincdf ([x, NaN], single (1), 0.5), single ([y, NaN]))
%!assert (nbincdf ([x, NaN], 1, single (0.5)), single ([y, NaN]))

## Test input validation
%!error<nbincdf: function called with too few input arguments.> nbincdf ()
%!error<nbincdf: function called with too few input arguments.> nbincdf (1)
%!error<nbincdf: function called with too few input arguments.> nbincdf (1, 2)
%!error<nbincdf: invalid argument for upper tail.> nbincdf (1, 2, 3, 4)
%!error<nbincdf: invalid argument for upper tail.> nbincdf (1, 2, 3, "some")
%!error<nbincdf: X, R, and PS must be of common size or scalars.> ...
%! nbincdf (ones (3), ones (2), ones (2))
%!error<nbincdf: X, R, and PS must be of common size or scalars.> ...
%! nbincdf (ones (2), ones (3), ones (2))
%!error<nbincdf: X, R, and PS must be of common size or scalars.> ...
%! nbincdf (ones (2), ones (2), ones (3))
%!error<nbincdf: X, R, and PS must not be complex.> nbincdf (i, 2, 2)
%!error<nbincdf: X, R, and PS must not be complex.> nbincdf (2, i, 2)
%!error<nbincdf: X, R, and PS must not be complex.> nbincdf (2, 2, i)
