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
## @deftypefn {Private Function} {@var{s} =} madsigma (@var{r}, @var{p})
##
## Robust scale from the median absolute deviation of the residuals @var{r},
## dropping the @math{@var{p}-1} smallest to account for the fitted parameters.
##
## @seealso{robustfit, nlinfit}
## @end deftypefn

function s = madsigma (r, p)

  rs = sort (abs (r));
  s = median (rs(max (1, p):end)) / 0.6745;

endfunction
