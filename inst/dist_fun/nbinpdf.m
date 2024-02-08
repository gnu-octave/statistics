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
## @deftypefn  {statistics} {@var{y} =} nbinpdf (@var{x}, @var{r}, @var{ps})
##
## Negative binomial probability density function (PDF).
##
## For each element of @var{x}, compute the probability density function (PDF)
## at @var{x} of the negative binomial distribution with parameters @var{r} and
## @var{ps}, where @var{r} is the number of successes until the experiment is
## stopped and @var{ps} is the probability of success in each experiment, given
## the number of failures in @var{x}.  The size of @var{y} is the common size of
## @var{x}, @var{r}, and @var{ps}.  A scalar input functions as a constant
## matrix of the same size as the other inputs.
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
## @seealso{nbininv, nbininv, nbinrnd, nbinfit, nbinlike, nbinstat}
## @end deftypefn

function y = nbinpdf (x, r, ps)

  ## Check for valid number of input arguments
  if (nargin < 3)
    error ("nbinpdf: function called with too few input arguments.");
  endif

  ## Check for common size of X, R, and PS
  if (! isscalar (x) || ! isscalar (r) || ! isscalar (ps))
    [retval, x, r, ps] = common_size (x, r, ps);
    if (retval > 0)
      error ("nbinpdf: X, R, and PS must be of common size or scalars.");
    endif
  endif

  ## Check for X, R, and PS being reals
  if (iscomplex (x) || iscomplex (r) || iscomplex (ps))
    error ("nbinpdf: X, R, and PS must not be complex.");
  endif

  ## Check for class type
  if (isa (x, "single") || isa (r, "single") || isa (ps, "single"))
    y = NaN (size (x), "single");
  else
    y = NaN (size (x));
  endif

  ok = (x < Inf) & (x == fix (x)) & (r > 0) & (r < Inf) & (ps >= 0) & (ps <= 1);

  k = (x < 0) & ok;
  y(k) = 0;

  k = (x >= 0) & ok;
  if (isscalar (r) && isscalar (ps))
    y(k) = bincoeff (-r, x(k)) .* (ps ^ r) .* ((ps - 1) .^ x(k));
  else
    y(k) = bincoeff (-r(k), x(k)) .* (ps(k) .^ r(k)) .* ((ps(k) - 1) .^ x(k));
  endif


endfunction

%!demo
%! ## Plot various PDFs from the negative binomial distribution
%! x = 0:40;
%! y1 = nbinpdf (x, 2, 0.15);
%! y2 = nbinpdf (x, 5, 0.2);
%! y3 = nbinpdf (x, 4, 0.4);
%! y4 = nbinpdf (x, 10, 0.3);
%! plot (x, y1, "*r", x, y2, "*g", x, y3, "*k", x, y4, "*m")
%! grid on
%! xlim ([0, 40])
%! ylim ([0, 0.12])
%! legend ({"r = 2, ps = 0.15", "r = 5, ps = 0.2", "r = 4, p = 0.4", ...
%!          "r = 10, ps = 0.3"}, "location", "northeast")
%! title ("Negative binomial PDF")
%! xlabel ("values in x (number of failures)")
%! ylabel ("density")

## Test output
%!shared x, y
%! x = [-1 0 1 2 Inf];
%! y = [0 1/2 1/4 1/8 NaN];
%!assert (nbinpdf (x, ones (1,5), 0.5*ones (1,5)), y)
%!assert (nbinpdf (x, 1, 0.5*ones (1,5)), y)
%!assert (nbinpdf (x, ones (1,5), 0.5), y)
%!assert (nbinpdf (x, [0 1 NaN 1.5 Inf], 0.5), [NaN 1/2 NaN 1.875*0.5^1.5/4 NaN], eps)
%!assert (nbinpdf (x, 1, 0.5*[-1 NaN 4 1 1]), [NaN NaN NaN y(4:5)])
%!assert (nbinpdf ([x, NaN], 1, 0.5), [y, NaN])

## Test class of input preserved
%!assert (nbinpdf (single ([x, NaN]), 1, 0.5), single ([y, NaN]))
%!assert (nbinpdf ([x, NaN], single (1), 0.5), single ([y, NaN]))
%!assert (nbinpdf ([x, NaN], 1, single (0.5)), single ([y, NaN]))

## Test input validation
%!error<nbinpdf: function called with too few input arguments.> nbinpdf ()
%!error<nbinpdf: function called with too few input arguments.> nbinpdf (1)
%!error<nbinpdf: function called with too few input arguments.> nbinpdf (1, 2)
%!error<nbinpdf: X, R, and PS must be of common size or scalars.> ...
%! nbinpdf (ones (3), ones (2), ones (2))
%!error<nbinpdf: X, R, and PS must be of common size or scalars.> ...
%! nbinpdf (ones (2), ones (3), ones (2))
%!error<nbinpdf: X, R, and PS must be of common size or scalars.> ...
%! nbinpdf (ones (2), ones (2), ones (3))
%!error<nbinpdf: X, R, and PS must not be complex.> nbinpdf (i, 2, 2)
%!error<nbinpdf: X, R, and PS must not be complex.> nbinpdf (2, i, 2)
%!error<nbinpdf: X, R, and PS must not be complex.> nbinpdf (2, 2, i)
