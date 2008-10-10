## Copyright (C) 2001 Paul Kienzle
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
## along with this program; If not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn {Function File} harmmean (@var{x})
## @deftypefnx{Function File} harmmean (@var{x}, @var{dim})
## Compute the harmonic mean.
##
## This function does the same as @code{mean (x, "h")}.
##
## @seealso{mean}
## @end deftypefn

function a = harmmean(x, dim)
  if (nargin == 1)
    a = mean(s, "h");
  else
    a = mean(x, "h", dim);
  endif
endfunction
