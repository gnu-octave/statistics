# Copyright (C) 2008   Sylvain Pelissier   <sylvain.pelissier@gmail.com>
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
## @deftypefn {Function File} {[@var{y}] =} nanvar(@var{X} [, @var{opt} [, @var{dim}]])
## nanstd is identical to the var function except that NaN values are ignored.
## @seealso{nanmean,nanstd,nanmax,nanmin}
## @end deftypefn

function y = nanvar(x,w,dim)
  if nargin < 1
    print_usage ();
  else
    if ((nargin < 2) || isempty(w))
      w = 0;
    endif
    
    if nargin < 3
      dim = min(find(size(x)>1));
      if isempty(dim)
	dim=1;
      endif
    endif
    
    y = nanstd(x,w,dim).^2;
  endif
endfunction

## Tests
%!shared x
%! x = [1 2 nan 3 4 5];
%!assert (nanvar (x), var (x(! isnan (x))), 10*eps)
