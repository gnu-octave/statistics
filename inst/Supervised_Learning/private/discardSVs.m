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
## @deftypefn {Private Function} {@var{m} =} discardSVs (@var{m})
## Collapse a linear LIBSVM model onto the one vector that decides it.
##
## A linear kernel makes the decision function @math{sum_i c_i <s_i, x> - rho},
## which is @math{<w, x> - rho} for @math{w = sum_i c_i s_i}.  One support
## vector holding @math{w}, with a coefficient of one, therefore decides
## exactly what all of them did, and the rest need not be kept.
##
## This is what lets the support vectors be discarded without touching a
## single scoring path: @code{svmpredict} goes on being the engine, over a
## model of one vector instead of many.  Emptying the properties alone would
## save nothing, the engine's own copy being where they are kept.
##
## The model must be linear.  Callers check that, so that the error names the
## method the user called.  Collapsing an already collapsed model returns it
## unchanged, one vector of @math{w} summing to @math{w}.
## @end deftypefn

function m = discardSVs (m)

  w = m.sv_coef' * m.SVs;
  m.SVs = sparse (w);
  m.sv_coef = 1;
  m.totalSV = 1;
  m.nSV = [1; 0];
  m.sv_indices = 1;

endfunction
