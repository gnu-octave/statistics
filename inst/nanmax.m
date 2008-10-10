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
## @deftypefn {Function File} {[@var{v}, @var{idx}] =} nanmax (@var{X})
## @deftypefnx{Function File} {[@var{v}, @var{idx}] =} nanmax (@var{X}, @var{Y})
## Find the maximal element while ignoring NaN values.
##
## @code{nanmax} is identical to the @code{max} function except that NaN values
## are ignored.  If all values in a column are NaN, the maximum is 
## returned as NaN rather than []. 
##
## @seealso{max, nansum, nanmin, nanmean, nanmedian}
## @end deftypefn

function [v, idx] = nanmax (X, Y, DIM) 
  if nargin < 1 || nargin > 3
    usage ("[v, idx] = nanmax(X [, Y, [DIM]])");
  elseif nargin == 1 || (nargin == 2 && isempty(Y))
    nanvals = isnan(X);
    X(nanvals) = -Inf;
    v = max (X);
    v(all(nanvals)) = NaN;
  elseif (nargin == 3 && isempty(Y))
    nanvals = isnan(X);
    X(nanvals) = -Inf;
    v = max (X,[],DIM);
    v(all(nanvals)) = NaN;
  else
    Xnan = isnan(X);
    Ynan = isnan(Y);
    X(Xnan) = -Inf;
    Y(Ynan) = -Inf;
    if (nargin == 3)
      v = max(X,Y,DIM);
    else
      v = max(X,Y);
    endif
    v(Xnan & Ynan) = NaN;
  endif
endfunction
