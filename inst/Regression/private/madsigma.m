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
## @deftypefn  {Private Function} {@var{s} =} madsigma (@var{r}, @var{p})
## @deftypefnx {Private Function} {@var{s} =} madsigma (@var{r}, @var{p}, @var{y})
##
## Robust scale from the median absolute deviation of the residuals @var{r},
## dropping the @math{@var{p}-1} smallest to account for the fitted parameters.
##
## Given the response @var{y}, the scale is floored at @code{1e-6 * std (y)},
## which is what an iteration needs and what the reported scale must not have.
##
## @seealso{robustfit, nlinfit}
## @end deftypefn

function s = madsigma (r, p, y)

  rs = sort (abs (r));
  s = median (rs(max (1, p):end)) / 0.6745;

  ## The scale that weights an iteration is floored relative to the response,
  ## so a fit that is already exact does not grade its weights out of rounding
  ## noise: on y = 3x through the origin the residuals are 1e-15 and the
  ## unfloored scale is 1e-15 with them, leaving the ratios order one and the
  ## weights graded where MATLAB returns ones.  Measured on R2024a, on a
  ## fixture where std (Y), norm (Y) and max (abs (Y)) differ by six orders of
  ## magnitude: the floor is 1e-6 * std (Y).  The scale reported back to a
  ## caller is the unfloored one, so the covariance step asks for two
  ## arguments and only an iteration asks for three.
  if (nargin > 2)
    tiny = 1e-6 * std (y);
    if (tiny == 0)
      tiny = 1;
    endif
    s = max (s, tiny);
  endif

endfunction
