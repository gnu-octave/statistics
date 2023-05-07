## Copyright (C) 2012 Rik Wehbring
## Copyright (C) 1995-2016 Kurt Hornik
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
## @deftypefn  {statistics} {@var{p} =} nbincdf (@var{x}, @var{n}, @var{ps})
## @deftypefnx {statistics} {@var{p} =} nbincdf (@var{x}, @var{n}, @var{ps}, @qcode{"upper"})
##
## Negative binomial cumulative distribution function (CDF).
##
## For each element of @var{x}, compute the cumulative distribution function
## (CDF) at @var{x} of the negative binomial distribution with parameters
## @var{n} and @var{ps}.  The size of @var{p} is the common size of @var{x},
## @var{n}, and @var{ps}.  A scalar input functions as a constant matrix of the
## same size as the other inputs.
##
## The algorithm uses the cumulative sums of the binomial masses.
##
## @code{@var{p} = nbincdf (@var{x}, @var{n}, @var{ps}, "upper")} computes the
## upper tail probability of the lognormal distribution.
##
## When @var{n} is integer this is the Pascal distribution.
## When @var{n} is extended to real numbers this is the Polya distribution.
##
## The number of failures in a Bernoulli experiment with success probability
## @var{ps} before the @var{n}-th success follows this distribution.
##
## @seealso{nbininv, nbinpdf, nbinrnd, nbinstat}
## @end deftypefn

function p = nbincdf (x, n, ps, uflag)

  ## Check for valid number of input arguments
  if (nargin < 3 || nargin > 4)
    error ("nbincdf: invalid number of input arguments.");
  endif

  ## Check for "upper" flag
  if (nargin == 4 && strcmpi (uflag, "upper"))
    uflag = true;
  elseif (nargin == 4  && ! strcmpi (uflag, "upper"))
    error ("nbincdf: invalid argument for upper tail.");
  else
    uflag = false;
  endif

  ## Check for N and PS being scalars
  scalarNPS = (isscalar(n) & isscalar(ps));

  ## Check for common size of X, N, and PS
  if (! isscalar (x) || ! isscalar (n) || ! isscalar (ps))
    [retval, x, n, ps] = common_size (x, n, ps);
    if (retval > 0)
      error ("nbincdf: X, N, and PS must be of common size or scalars.");
    endif
  endif

  ## Check for X, N, and PS being reals
  if (iscomplex (x) || iscomplex (n) || iscomplex (ps))
    error ("nbincdf: X, N, and PS must not be complex.");
  endif

  ## Check for appropriate class
  if (isa (x, "single") || isa (n, "single") || isa (ps, "single"))
    p = zeros (size (x), "single");
  else
    p = zeros (size (x));
  endif

  ## Force NaN for out of range or missing parameters and missing data NaN
  is_nan = (isnan (x) | isnan (n) | (n <= 0) | (n == Inf) | (ps < 0) | (ps > 1));
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
      p(k) = betainc (ps(k), n(k), xf(k) + 1, "upper");
    else
      max_val = max (xf(k));
      if (scalarNPS)
        tmp = cumsum (nbinpdf (0:max_val, n(1), ps(1)));
        p(k) = tmp(xf(k) + 1);
      else
        idx = (0:max_val)';
        compare = idx(:, ones (size (k)));
        index = xf(k);
        index = index(:);
        index = index(:, ones (size (idx)))';
        n_big = n(k);
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
%!error nbincdf ()
%!error nbincdf (1)
%!error nbincdf (1, 2)
%!error nbincdf (1, 2, 3, 4, 5)
%!error<nbincdf: invalid argument for upper tail.> nbincdf (1, 2, 3, 4)
%!error<nbincdf: invalid argument for upper tail.> nbincdf (1, 2, 3, "some")
%!error nbincdf (ones (3), ones (2), ones (2))
%!error nbincdf (ones (2), ones (3), ones (2))
%!error nbincdf (ones (2), ones (2), ones (3))
%!error nbincdf (i, 2, 2)
%!error nbincdf (2, i, 2)
%!error nbincdf (2, 2, i)
