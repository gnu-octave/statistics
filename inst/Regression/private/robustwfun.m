## Copyright (C) 2026 Andreas Bertsatos <abertsatos@biol.uoa.gr>
## Copyright (C) 2026 Avanish Salunke <avanishsalunke16@gmail.com>
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
## @deftypefn {Private Function} {[@var{w}, @var{psi}, @var{psip}] =} robustwfun (@var{z}, @var{wfun})
##
## Weight, influence function and its derivative for a robust weight function.
##
## @var{z} is the vector of residuals already scaled by the robust scale and the
## tuning constant, and @var{wfun} is one of @qcode{"andrews"},
## @qcode{"bisquare"}, @qcode{"cauchy"}, @qcode{"fair"}, @qcode{"huber"},
## @qcode{"logistic"}, @qcode{"ols"}, @qcode{"talwar"} or @qcode{"welsch"}, or a
## function handle returning the weights.  @var{w} is the weight, @var{psi} is
## @code{@var{z} .* @var{w}} and @var{psip} is its derivative, computed
## analytically for a named function and by central difference for a handle.
##
## @qcode{"andrews"} and @qcode{"logistic"} take their limiting value of 1 at
## @code{@var{z} == 0}, where the quotient defining them is @math{0/0}.
##
## @seealso{robustfit, nlinfit, robusttune}
## @end deftypefn

function [w, psi, psip] = robustwfun (z, wfun)

  if (is_function_handle (wfun))
    w = wfun (z);
    if (nargout > 1)
      d = 1e-6;
      psi = z .* w;
      psip = ((z + d) .* wfun (z + d) - (z - d) .* wfun (z - d)) / (2 * d);
    endif
    return;
  endif
  a = abs (z);
  switch (wfun)
    case "andrews"
      in = a < pi;
      w = (sin (z) ./ (z + (z == 0))) .* in;
      w(z == 0) = 1;
      if (nargout > 1); psi = sin (z) .* in;  psip = cos (z) .* in;  endif
    case "bisquare"
      in = a < 1;
      w = (1 - z .^ 2) .^ 2 .* in;
      if (nargout > 1)
        psi = z .* w;  psip = (1 - z .^ 2) .* (1 - 5 * z .^ 2) .* in;
      endif
    case "cauchy"
      w = 1 ./ (1 + z .^ 2);
      if (nargout > 1)
        psi = z .* w;  psip = (1 - z .^ 2) ./ (1 + z .^ 2) .^ 2;
      endif
    case "fair"
      w = 1 ./ (1 + a);
      if (nargout > 1); psi = z .* w;  psip = 1 ./ (1 + a) .^ 2;  endif
    case "huber"
      w = 1 ./ max (1, a);
      if (nargout > 1); psi = z .* w;  psip = double (a < 1);  endif
    case "logistic"
      th = tanh (z);
      w = th ./ (z + (z == 0));
      w(z == 0) = 1;
      if (nargout > 1); psi = th;  psip = 1 - th .^ 2;  endif
    case "ols"
      w = ones (size (z));
      if (nargout > 1); psi = z;  psip = ones (size (z));  endif
    case "talwar"
      in = a < 1;
      w = double (in);
      if (nargout > 1); psi = z .* in;  psip = double (in);  endif
    case "welsch"
      e = exp (- z .^ 2);
      w = e;
      if (nargout > 1); psi = z .* e;  psip = (1 - 2 * z .^ 2) .* e;  endif
  endswitch

endfunction
