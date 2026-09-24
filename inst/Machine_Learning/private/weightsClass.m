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
## FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
## more details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn {Private Function} {@var{errmsg} =} weightsClass (@var{W})
## Check that observation weights are single or double.
##
## @var{errmsg} is the body of the message the caller should raise, or empty
## when @var{W} is empty or a real single or double array.  Logical and integer
## weights are refused, as MATLAB refuses them.  The other checks on the
## weights stay with the caller.
## @end deftypefn

function errmsg = weightsClass (W)

  errmsg = "";
  if (! isempty (W) && ! (isfloat (W) && isreal (W)))
    errmsg = "'Weights' must be a real vector of class single or double.";
  endif

endfunction
