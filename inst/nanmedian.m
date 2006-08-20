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
## along with this program; if not, write to the Free Software
## Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

## v = nanmedian(X [, dim]);
## nanmedian is identical to the median function except that NaN values are
## ignored.  If all values are NaN, the median is returned as NaN. 
## [Is this behaviour compatible?]
##
## See also: nanmin, nanmax, nansum, nanmean
function v = nanmedian (X, varargin)
  if nargin < 1 || nargin > 2
    usage ("v = nanmedian(X [, dim])");
  endif
  if nargin < 2
    dim = min(find(size(X)>1));
    if isempty(dim), dim=1; endif;
  else
    dim = varargin{:};
  endif

  sz = size (X);
  if (prod (sz) > 1)
    try dfi = do_fortran_indexing;
    catch dfi = 0;
    end
    try wfi = warn_fortran_indexing;
    catch wfi = 0;
    end
    unwind_protect
      do_fortran_indexing = 1;
      warn_fortran_indexing = 0;
      ## Find lengths of datasets after excluding NaNs; valid datasets
      ## are those that are not empty after you remove all the NaNs
      n = sz(dim) - sum (isnan(X),varargin{:});

      ## When n is equal to zero, force it to one, so that median
      ## picks up a NaN value below
      n (n==0) = 1;

      ## Sort the datasets, with the NaN going to the end of the data
      X = sort (X, varargin{:});

      ## Determine the offset for each column in single index mode
      colidx = reshape((0:(prod(sz) / sz(dim) - 1)), size(n)); 
      colidx = floor(colidx / prod(sz(1:dim-1))) * prod(sz(1:dim)) + ...
	  mod(colidx,prod(sz(1:dim-1)));
      stride = prod(sz(1:dim-1));

      ## Average the two central values of the sorted list to compute
      ## the median, but only do so for valid rows.  If the dataset
      ## is odd length, the single central value will be used twice.
      ## E.g., 
      ##   for n==5, ceil(2.5+0.5) is 3 and floor(2.5+0.5) is also 3
      ##   for n==6, ceil(3.0+0.5) is 4 and floor(3.0+0.5) is 3
      ## correction made for stride of data "stride*ceil(2.5-0.5)+1"
      v = (X(colidx + stride*ceil(n./2-0.5) + 1)  + ...
	   X(colidx + stride*floor(n./2-0.5) + 1)) ./ 2;
    unwind_protect_cleanup
      do_fortran_indexing = dfi;
      warn_fortran_indexing = wfi;
    end_unwind_protect
  else
    error ("nanmedian: invalid matrix argument");
  endif
endfunction
