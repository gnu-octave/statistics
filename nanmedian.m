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
function v = nanmedian (X, dim)
  if nargin < 1 || nargin > 2
    usage ("v = nanmean(X [, dim])");
  else
    if nargin == 1
      if size(X,1) == 1
	dim = 2; 
      else
        dim = 1;
      endif
    endif
    if (dim == 2) X = X.'; endif
    dfi = do_fortran_indexing;
    unwind_protect
      do_fortran_indexing = 1;
      ## Find lengths of datasets after excluding NaNs; valid datasets
      ## are those that are not empty after you remove all the NaNs
      n = size(X,1) - sum (isnan(X));
      valid = find(n!=0);

      ## Extract all non-empty datasets and sort, replacing NaN with Inf
      ## so that the invalid elements go toward the ends of the columns
      X (isnan(X)) = Inf;
      X = sort ( X (:, valid) );

      ## Determine the offset for each remaining column in single index mode
      colidx = (0:size(X,2)-1)*size(X,1);

      ## Assume the median for all datasets will be NaNs
      v = NaN*ones(size(n));

      ## Average the two central values of the sorted list to compute
      ## the median, but only do so for valid rows.  If the dataset
      ## is odd length, the single central value will be used twice.
      ## E.g., 
      ##   for n==5, ceil(2.5+0.4) is 3 and floor(2.5+0.6) is also 3
      ##   for n==6, ceil(3.0+0.4) is 4 and floor(3.0+0.6) is 3
      v(valid) = ( X (colidx + floor(n(valid)./2+0.6)) ...
		 + X (colidx + ceil(n(valid)./2+0.4)) ) ./ 2;
    unwind_protect_cleanup
      do_fortran_indexing = dfi;
    end_unwind_protect
    if (dim == 2) v = v.'; endif
  endif
endfunction
