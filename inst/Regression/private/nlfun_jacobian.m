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
## @deftypefn {statistics} {@var{J} =} nlfun_jacobian (@var{modelfun}, @var{beta}, @var{X}, @var{derivstep})
##
## Numeric Jacobian of @var{modelfun} with respect to @var{beta} by forward
## differences, evaluated at @var{X}.  Internal helper shared by @code{nlinfit}
## and @code{nlpredci}; not intended to be called directly.
##
## @end deftypefn

function J = nlfun_jacobian (modelfun, beta, X, derivstep)

  beta = beta(:);
  f0   = modelfun (beta, X);
  f0   = f0(:);
  n    = numel (f0);
  p    = numel (beta);
  J    = zeros (n, p);
  for j = 1:p
    h = derivstep * max (abs (beta(j)), 1);
    if (h == 0)
      h = derivstep;
    endif
    bj     = beta;
    bj(j)  = bj(j) + h;
    fj     = modelfun (bj, X);
    J(:,j) = (fj(:) - f0) / h;
  endfor

endfunction
