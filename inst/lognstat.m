## Copyright (C) 2006 Arno Onken
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

## -*- texinfo -*-
## @deftypefn {Function File} {[@var{m}, @var{v}] =} lognstat (@var{mu}, @var{sigma})
## Returns mean and variance of the lognormal distribution
##
## Arguments are
##
## @itemize
## @item
## @var{mu} is the first parameter of the lognormal distribution
## @item
## @var{sigma} is the second parameter of the lognormal distribution.
## @var{sigma} must be positive or zero
## @end itemize
## @var{mu} and @var{sigma} must be of common size or one of them must be
## scalar
##
## Return values are
##
## @itemize
## @item
## @var{m} is the mean of the lognormal distribution
## @item
## @var{v} is the variance of the lognormal distribution
## @end itemize
##
## Examples:
##
## @example
## mu = 0:0.2:1;
## sigma = 0.2:0.2:1.2;
## [m, v] = lognstat (mu, sigma)
##
## [m, v] = lognstat (0, sigma)
## @end example
##
## References:
##
## @itemize
## @item
## @cite{Matlab 7.0 documentation (pdf)}
## @item
## @uref{http://en.wikipedia.org/wiki/Lognormal_distribution}
## @end itemize
##
## @end deftypefn

## Author: Arno Onken <whyly@gmx.net>

function [m, v] = lognstat (mu, sigma)

  # Check arguments
  if (nargin != 2)
    usage ("[m, v] = lognstat (mu, sigma)");
  endif

  if (! isempty (mu) && ! ismatrix (mu))
    error ("lognstat: mu must be a numeric matrix");
  endif
  if (! isempty (sigma) && ! ismatrix (sigma))
    error ("lognstat: sigma must be a numeric matrix");
  endif

  if (! isscalar (mu) || ! isscalar (sigma))
    [retval, mu, sigma] = common_size (mu, sigma);
    if (retval > 0)
      error ("lognstat: mu and sigma must be of common size or scalar");
    endif
  endif

  # Calculate moments
  m = exp (mu + (sigma .^ 2) ./ 2);
  v = (exp (sigma .^ 2) - 1) .* exp (2 .* mu + sigma .^ 2);

  # Continue argument check
  k = find (! (sigma >= 0) | ! (sigma < Inf));
  if (any (k))
    m (k) = NaN;
    v (k) = NaN;
  endif

endfunction

%!test
%! mu = 0:0.2:1;
%! sigma = 0.2:0.2:1.2;
%! [m, v] = lognstat (mu, sigma);
%! expected_m = [1.0202, 1.3231, 1.7860, 2.5093,  3.6693,   5.5845];
%! expected_v = [0.0425, 0.3038, 1.3823, 5.6447, 23.1345, 100.4437];
%! assert (m, expected_m, 0.001);
%! assert (v, expected_v, 0.001);

%!test
%! sigma = 0.2:0.2:1.2;
%! [m, v] = lognstat (0, sigma);
%! expected_m = [1.0202, 1.0833, 1.1972, 1.3771, 1.6487,  2.0544];
%! expected_v = [0.0425, 0.2036, 0.6211, 1.7002, 4.6708, 13.5936];
%! assert (m, expected_m, 0.001);
%! assert (v, expected_v, 0.001);
