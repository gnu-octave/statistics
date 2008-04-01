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
## @deftypefn {Function File}  {@var{a} =} prctile (@var{x}, @var{p})
## @deftypefnx {Function File}  {@var{a} =} prctile (@var{x}, @var{p}, @var{dim})
## Computes the value associated with the @var{p}-th percentile of
## @var{x}. If @var{x} is a matrix, computes @var{p} for each column of
## @var{x}. If @var{p} is a vector, the returned value is a matrix with
## one row for each element of @var{p} and one column for each column of
## @var{x}. If @var{x} is a vector, then @var{y} is a vector with the
## same size as @var{p}. If the argument @var{dim} is given, operate
## along this dimesnion rather than the columns of the matrix.
##
## The first and last values are pegged at 0 percent and 100 percent
## respectively, and the rest of the values are uniformly spaced between
## them, with linear interpolation between the points.  This is
## consistent with the definition of quantile given in the R statistics
## package, but inconsistent with that of the statistics toolbox from
## Matlab.
## @end deftypefn

function a = prctile(x, p, dim)

  if (nargin != 2 && nargin != 3)
    print_usage();
  endif

  nd = ndims (x);
  sz = size (x);

  if (nargin == 2)
    ## Find the first non-singleton dimension.
    dim  = 1;
    while (dim < nd + 1 && sz(dim) == 1)
      dim = dim + 1;
    endwhile
    if (dim > nd)
      dim = 1;
    endif
  else
    if (! (isscalar (dim) && dim == round (dim)) && dim > 0
	&& dim < (nd + 1))
      error ("prctile: dim must be an integer and valid dimension");
    endif
  endif
  perm = 1:nd;
  perm(1) = dim;
  perm(dim) = 1;
  x = permute (x, perm);
  y = sort (x);
  trim = 1 + (size (y, 1) - 1) * p(:) * 0.01;
  delta = (trim - floor (trim)) * ones (1, size (y, 2));
  a = y (floor (trim), :) .* (1 - delta) + y (ceil (trim), :) .* delta;
  if (isvector (a))
    a = reshape (a, size (p));
  else
    a = ipermute (a, perm);
  endif

endfunction

%!assert(prctile([0,1],[25,75]),[.25,.75],eps)
%!assert(prctile([0,0;1,2],[25,75]),[.25,.5;.75,1.5],eps)
%!assert(prctile([0,1;0,2],[25,75],2),[.25,.75;.5,1.5],eps)
