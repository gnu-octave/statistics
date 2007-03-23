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
## Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA

## v = nanstd(X [, opt [, dim]]);
## nanstd is identical to the std function except that NaN values are
## ignored.  If all values are NaN, the std is returned as NaN. If there
## is only a single non-NaN value, the std is returned as 0. 
##
##	0:  normalizes with N-1 [default]
##		provides the square root of best unbiased estimator of the variance
##	1:  normalizes with N,
##		this provides the square root of the second moment around the mean
##
## See also: nanmin, nanmax, nansum, nanmedian, nanmean
function v = nanstd (X, opt, varargin)
  if nargin < 1
    usage ("v = nanstd(X [, opt [, dim]])");
  else
    if nargin < 3
      dim = min(find(size(X)>1));
      if isempty(dim), dim=1; endif;
    else
      dim = varargin{1};
    endif
    if ((nargin < 2) || isempty(opt))
      opt = 0;
    endif

    ## determine the number of non-missing points in each data set
    n = sum (!isnan(X), varargin{:});
    
    ## replace missing data with zero and compute the mean
    X(isnan(X)) = 0;
    meanX = sum (X, varargin{:}) ./ n;
    
    ## subtract the mean from the data and compute the sum squared
    sz = ones(1,length(size(X)));
    sz(dim) = size(X,dim);
    v = sumsq (X - repmat(meanX,sz), varargin{:});
    
    ## because the missing data was set to zero each missing data
    ## point will contribute (-meanX)^2 to sumsq, so remove these
    v = v - (meanX .^ 2) .* (size(X,dim) - n);
    
    if (opt == 0)
      ## compute the standard deviation from the corrected sumsq using
      ## max(n-1,1) in the denominator so that the std for a single point is 0
      v = sqrt ( v ./ max(n - 1, 1) );
    elseif (opt == 1)
      ## compute the standard deviation from the corrected sumsq
      v = sqrt ( v ./ n );
    else
      error ("std: unrecognized normalization type");
    endif

  endif
endfunction
