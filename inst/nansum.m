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
## @deftypefn  {Built-in Function} {} nansum (@var{x})
## @deftypefnx {Built-in Function} {} nansum (@var{x}, @var{dim})
## @deftypefnx {Built-in Function} {} nansum (@dots{}, @qcode{"native"})
## @deftypefnx {Built-in Function} {} nansum (@dots{}, @qcode{"double"})
## @deftypefnx {Built-in Function} {} nansum (@dots{}, @qcode{"extra"})
## Compute the sum while ignoring NaN values.
##
## @code{nansum} is identical to the @code{sum} function except that NaN
## values are treated as 0 and so ignored.  If all values are NaN, the sum is
## returned as 0.
##
## See help text of @code{sum} for details on the options.
##
## @seealso{sum, nanmin, nanmax, nanmean, nanmedian}
## @end deftypefn

function v = nansum (X, varargin)
  if (nargin < 1)
    print_usage ();
  else
    X(isnan (X)) = 0;
    v = sum (X, varargin{:});
  endif
endfunction

%!assert (nansum ([2 4 NaN 7]), 13)
%!assert (nansum ([2 4 NaN Inf]), Inf)

%!assert (nansum ([1 NaN 3; NaN 5 6; 7 8 NaN]), [8 13 9])
%!assert (nansum ([1 NaN 3; NaN 5 6; 7 8 NaN], 2), [4; 11; 15])
%!assert (nansum (single ([1 NaN 3; NaN 5 6; 7 8 NaN])), single ([8 13 9]))
%!assert (nansum (single ([1 NaN 3; NaN 5 6; 7 8 NaN]), "double"), [8 13 9])

%!assert (nansum (uint8 ([2 4 1 7])), 14)
%!assert (nansum (uint8 ([2 4 1 7]), "native"), uint8 (14))
%!assert (nansum (uint8 ([2 4 1 7])), 14)
