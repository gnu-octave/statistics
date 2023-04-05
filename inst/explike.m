## Copyright (C) 2021 Nir Krakauer <nkrakauer@ccny.cuny.edu>
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
## @deftypefn  {statistics} {[@var{nlogL}, @var{avar}] =} explike (@var{mu}, @var{x})
##
## Negative log-likelihood for the exponential distribution.
##
## @subheading Arguments
##
## @itemize @bullet
## @item
## @var{mu} is a scalar containing the mean of the exponential distribution
## (@math{@var{mu} = 1 / λ}, where λ is the rate, or inverse scale parameter).
## @item
## @var{x} is the vector of given values.
##
## @end itemize
##
## @subheading Return values
##
## @itemize @bullet
## @item
## @var{nlogL} is the negative log-likelihood.
## @item
## @var{avar} is the inverse of the Fisher information matrix.
## (The Fisher information matrix is the second derivative of the negative log
## likelihood with respect to the parameter value.)
##
## @end itemize
##
## @seealso{expcdf, expinv, exppdf, exprnd, expfit, expstat}
## @end deftypefn

function [nlogL, avar] = explike (mu, x)

  ## Check arguments
  if (nargin != 2)
    print_usage;
  endif

  beta = mu(1);
  n = numel (x);
  sx = sum (x(:));
  sxb = sx/beta;

  ## Calculate negative log likelihood
  nlogL = sxb + n*log(beta);

  ## Optionally calculate the inverse (reciprocal) of the second derivative
  ## of the negative log likelihood with respect to parameter
  if (nargout > 1)
    avar = (beta ^ 2) ./ (2 * sxb - n);
  endif

endfunction


%!test
%! x = 12;
%! beta = 5;
%! [L, V] = explike (beta, x);
%! expected_L = 4.0094;
%! expected_V = 6.5789;
%! assert (L, expected_L, 0.001);
%! assert (V, expected_V, 0.001);

%!test
%! x = 1:5;
%! beta = 2;
%! [L, V] = explike (beta, x);
%! expected_L = 10.9657;
%! expected_V = 0.4;
%! assert (L, expected_L, 0.001);
%! assert (V, expected_V, 0.001);
%!error explike ();
