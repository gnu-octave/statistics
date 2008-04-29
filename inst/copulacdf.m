## Copyright (C) 2008  Arno Onken
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn {Function File} {@var{p} =} copulacdf (@var{family}, @var{x}, @var{theta})
## Compute the cumulative distribution function of a copula family.
##
## @subheading Arguments
##
## @itemize @bullet
## @item
## @var{family} is the copula family name. Currently, @var{family} can
## be @code{'Clayton'} for the Clayton family, @code{'Gumbel'} for the
## Gumbel-Hougaard family, or @code{'Frank'} for the Frank family.
##
## @item
## @var{x} is the support where each row corresponds to an observation.
##
## @item
## @var{theta} is the parameter of the copula. The elements of
## @var{theta} must be greater than or equal to @code{-1} in the
## bivariate case and non-negative for greater dimensions for the
## Clayton family, greater than or equal to @code{1} for the
## Gumbel-Hougaard family, and arbitrary in the bivariate case and
## non-negative for greater dimensions for the Frank family. @var{theta}
## must be a column vector with the same number of rows as @var{x} or be
## scalar.
## @end itemize
##
## @subheading Return values
##
## @itemize @bullet
## @item
## @var{p} is the cumulative distribution of the copula at each row of
## @var{x} and corresponding parameter @var{theta}.
## @end itemize
##
## @subheading Examples
##
## @example
## @group
## x = [0.2:0.2:0.6; 0.2:0.2:0.6];
## theta = [1; 2];
## p = copulacdf ("Clayton", x, theta)
## @end group
##
## @group
## p = copulacdf ("Gumbel", x, 2)
## @end group
## @end example
##
## @subheading References
##
## @enumerate
## @item
## Roger B. Nelsen. @cite{An Introduction to Copulas}. Springer,
## New York, second edition, 2006.
## @end enumerate
## @end deftypefn

## Author: Arno Onken <asnelt@asnelt.org>
## Description: CDF of a copula family

function p = copulacdf (family, x, theta)

  # Check arguments
  if (nargin != 3)
    print_usage ();
  endif

  if (! ischar (family))
    error ("copulacdf: family must be one of 'Clayton', 'Gumbel', and 'Frank'");
  endif

  if (! isempty (x) && ! ismatrix (x))
    error ("copulacdf: x must be a numeric matrix");
  endif

  [n, d] = size (x);

  if (! isvector (theta) || (! isscalar (theta) && size (theta, 1) != n))
    error ("copulacdf: theta must be a column vector with the same number of rows as x or be scalar");
  endif

  if (n == 0)
    # Input is empty
    p = zeros (0, 1);
  else
    if (n > 1 && isscalar (theta))
      theta = repmat (theta, n, 1);
    endif

    # Truncate input to unit hypercube
    x(x < 0) = 0;
    x(x > 1) = 1;

    # Compute the cumulative distribution function according to family
    lowerarg = lower (family);

    if (strcmp (lowerarg, "clayton"))
      # The Clayton family
      p = exp (-log (sum (x .^ (repmat (-theta, 1, d)), 2) - d + 1) ./ theta);
      p(p < 0) = 0;

      k = find (theta == 0);
      if (any (k))
        # Product copula at columns k
        p(k) = prod (x(k, :), 2);
      endif

      # Check theta
      if (d > 2)
        k = find (! (theta >= 0) | ! (theta < inf));
      else
        # In the bivariate case, theta can be negative
        k = find (! (theta >= -1) | ! (theta < inf));
      endif
      if (any (k))
        p(k) = NaN;
      endif
    elseif (strcmp (lowerarg, "gumbel"))
      # The Gumbel-Hougaard family
      p = exp (-(sum ((-log (x)) .^ repmat (theta, 1, d), 2)) .^ (1 ./ theta));
      p(! (theta >= 1) | ! (theta < inf)) = NaN;
    elseif (strcmp (lowerarg, "frank"))
      # The Frank family
      p = -log (1 + (prod (expm1 (repmat (-theta, 1, d) .* x), 2)) ./ (expm1 (-theta) .^ (d - 1))) ./ theta;

      k = find (theta == 0);
      if (any (k))
        # Product copula at columns k
        p(k) = prod (x(k, :), 2);
      endif

      # Check theta
      if (d > 2)
        k = find (! (theta >= 0) | ! (theta < inf));
      else
        # In the bivariate case, theta can be negative
        k = find (! (theta > -inf) | ! (theta < inf));
      endif
      if (any (k))
        p(k) = NaN;
      endif
    else
      error ("copulacdf: unknown copula family '%s'", family);
    endif
  endif

endfunction

%!test
%! x = [0.2:0.2:0.6; 0.2:0.2:0.6];
%! theta = [1; 2];
%! p = copulacdf ("Clayton", x, theta);
%! expected_p = [0.1395; 0.1767];
%! assert (p, expected_p, 0.001);

%!test
%! x = [0.2:0.2:0.6; 0.2:0.2:0.6];
%! p = copulacdf ("Gumbel", x, 2);
%! expected_p = [0.1464; 0.1464];
%! assert (p, expected_p, 0.001);

%!test
%! x = [0.2:0.2:0.6; 0.2:0.2:0.6];
%! theta = [1; 2];
%! p = copulacdf ("Frank", x, theta);
%! expected_p = [0.0699; 0.0930];
%! assert (p, expected_p, 0.001);
