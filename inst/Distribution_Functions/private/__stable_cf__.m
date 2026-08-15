## Copyright (C) 2026 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn {Function File} {@var{phi} =} __stable_cf__ (@var{t}, @var{alpha}, @var{beta})
##
## Characteristic function of the standard stable distribution.
##
## Returns the characteristic function @var{phi} evaluated at the points @var{t}
## of the standard stable distribution @code{S(@var{alpha}, @var{beta}, 1, 0)}
## in the Nolan @qcode{S0} parameterization.  This is a private helper for
## @code{stblpdf} and @code{stblcdf}, which recover the density and the
## cumulative probability by numerical inversion.
##
## @end deftypefn

function phi = __stable_cf__ (t, alpha, beta)

  at = abs (t);
  st = sign (t);
  if (abs (alpha - 1) < eps)
    ## The |t|*log|t| term vanishes at the origin
    lt = at .* log (at);
    lt(at == 0) = 0;
    phi = exp (-at - 1i .* beta .* (2 ./ pi) .* st .* lt);
  else
    phi = exp (-at .^ alpha - 1i .* beta .* tan (pi .* alpha ./ 2) ...
               .* st .* (at - at .^ alpha));
  endif

endfunction
