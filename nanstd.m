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

## v = nanstd(X [, dim]);
## nanstd is identical to the std function except that NaN values are
## ignored.  If all values are NaN, the std is returned as NaN. If there
## is only a single non-NaN value, the std is returned as 0. 
## [Is this behaviour compatible?]
##
## See also: nanmin, nanmax, nansum, nanmedian, nanmean
function v = nanstd (X, dim)
  if nargin < 1
    usage ("v = nanstd(X [, dim])");
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
    wdz = warn_divide_by_zero;
    unwind_protect
      do_fortran_indexing = 1;
      warn_divide_by_zero = 0;

      ## determine the number of non-missing points in each data set
      n = sum (!isnan(X));

      ## replace missing data with zero and compute the mean
      X(isnan(X)) = 0;
      meanX = sum (X) ./ n;

      ## subtract the mean from the data and compute the sum squared
      v = sumsq (X - ones(size(X,1), 1) * meanX);

      ## because the missing data was set to zero each missing data
      ## point will contribute (-meanX)^2 to sumsq, so remove these
      v = v - (meanX .^ 2) .* (size(X,1) - n);

      ## compute the standard deviation from the corrected sumsq
      v = sqrt ( v ./ (n - 1) );

      ## set special values of std for n=0 and n=1
      ## v(n == 0) = NaN;  # meanX = 0/0 -> NaN above, so not necessary
      v(n == 1) = 0;
    unwind_protect_cleanup
      do_fortran_indexing = dfi;
      warn_divide_by_zero = wdz;
    end_unwind_protect
    if (dim == 2) v = v.'; endif
  endif
endfunction
