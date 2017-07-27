## Copyright (C) 2006, 2007 Arno Onken <asnelt@asnelt.org>
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
## @deftypefn {Function File} {[@var{m}, @var{v}] =} betastat (@var{a}, @var{b})
## Compute mean and variance of the beta distribution.
##
## @subheading Arguments
##
## @itemize @bullet
## @item
## @var{a} is the first parameter of the beta distribution. @var{a} must be
## positive
##
## @item
## @var{b} is the second parameter of the beta distribution. @var{b} must be
## positive
## @end itemize
## @var{a} and @var{b} must be of common size or one of them must be scalar
##
## @subheading Return values
##
## @itemize @bullet
## @item
## @var{m} is the mean of the beta distribution
##
## @item
## @var{v} is the variance of the beta distribution
## @end itemize
##
## @subheading Examples
##
## @example
## @group
## a = 1:6;
## b = 1:0.2:2;
## [m, v] = betastat (a, b)
## @end group
##
## @group
## [m, v] = betastat (a, 1.5)
## @end group
## @end example
##
## @subheading References
##
## @enumerate
## @item
## Wendy L. Martinez and Angel R. Martinez. @cite{Computational Statistics
## Handbook with MATLAB}. Appendix E, pages 547-557, Chapman & Hall/CRC,
## 2001.
##
## @item
## Athanasios Papoulis. @cite{Probability, Random Variables, and Stochastic
## Processes}. McGraw-Hill, New York, second edition, 1984.
## @end enumerate
## @end deftypefn

## Author: Arno Onken <asnelt@asnelt.org>
## Description: Moments of the beta distribution

function [m, v] = betastat (a, b)

  if (nargin != 2)
    print_usage ();
  elseif (! isscalar (a) && ! isscalar (b) && ! size_equal (a, b))
    error ("betastat: a and b must be of common size or scalar");
  endif

  k = find (! (a > 0 & b > 0));

  # Calculate moments
  a_b = a + b;
  m = a ./ (a_b);
  m(k) = NaN;

  if (nargout > 1)
    v = (a .* b) ./ ((a_b .^ 2) .* (a_b + 1));
    v(k) = NaN;
  endif

endfunction

%!test
%! a = -2:6;
%! b = 0.4:0.2:2;
%! [m, v] = betastat (a, b);
%! expected_m = [NaN NaN NaN 1/2 2/3.2 3/4.4 4/5.6 5/6.8 6/8];
%! expected_v = [NaN NaN NaN 0.0833, 0.0558, 0.0402, 0.0309, 0.0250, 0.0208];
%! assert (m, expected_m, eps*100);
%! assert (v, expected_v, 0.001);

%!test
%! a = -2:1:6;
%! [m, v] = betastat (a, 1.5);
%! expected_m = [NaN NaN NaN 1/2.5 2/3.5 3/4.5 4/5.5 5/6.5 6/7.5];
%! expected_v = [NaN NaN NaN 0.0686, 0.0544, 0.0404, 0.0305, 0.0237, 0.0188];
%! assert (m, expected_m);
%! assert (v, expected_v, 0.001);

%!test
%! a = [14  Inf   10  NaN  10];
%! b = [12    9  NaN  Inf  12];
%! [m, v] = betastat (a, b);
%! expected_m = [14/26 NaN NaN NaN 10/22];
%! expected_v = [168/18252 NaN NaN NaN 120/11132];
%! assert (m, expected_m);
%! assert (v, expected_v);

%!assert (nthargout (1:2, @betastat, 5, []), {[], []})
%!assert (nthargout (1:2, @betastat, [], 5), {[], []})
%!assert (nthargout (1:2, @betastat, "", 5), {[], []})
%!assert (nthargout (1:2, @betastat, true, 5), {1/6, 5/252})

%!assert (size (betastat (rand (10, 5, 4), rand (10, 5, 4))), [10 5 4])
%!assert (size (betastat (rand (10, 5, 4), 7)), [10 5 4])

