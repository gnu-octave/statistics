## Copyright (C) 2001 Paul Kienzle <pkienzle@users.sf.net>
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
## @deftypefn {Function File} {@var{v} =} nanmean (@var{X})
## @deftypefnx{Function File} {@var{v} =} nanmean (@var{X}, @var{dim})
## Compute the mean value while ignoring NaN values.
##
## This is a legacy function.  It is best to use the @code{mean} function with
## the "omitnan" flag. 
##
## @seealso{mean, nanmin, nanmax, nansum, nanmedian}
## @end deftypefn

function v = nanmean (X, varargin) 
  if nargin < 1
    print_usage;
  else
    n = sum (!isnan(X), varargin{:});
    n(n == 0) = NaN;
    X(isnan(X)) = 0;
    v = sum (X, varargin{:}) ./ n;
  endif
endfunction

%!test
%! x = [1 2 nan 3 4 5];
%! assert (nanmean (x), mean (x(! isnan (x)), "omitnan"), 10*eps)