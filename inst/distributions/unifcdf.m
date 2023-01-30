## Copyright (C) 2012 Rik Wehbring
## Copyright (C) 1995-2016 Kurt Hornik
## Copyright (C) 2023 Andreas Bertsatos <abertsatos@biol.uoa.gr>
##
## This file is part of the statistics package for GNU Octave.
##
## This program is free software: you can redistribute it and/or
## modify it under the terms of the GNU General Public License as
## published by the Free Software Foundation, either version 3 of the
## License, or (at your option) any later version.
##
## This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; see the file COPYING.  If not, see
## <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn  {statistics} @var{p} = unifcdf (@var{x})
## @deftypefnx {statistics} @var{p} = unifcdf (@var{x}, @var{a})
## @deftypefnx {statistics} @var{p} = unifcdf (@var{x}, @var{a}, @var{b})
## @deftypefnx {statistics} @var{p} = unifcdf (@dots{}, "upper")
##
## Uniform cumulative distribution function (CDF).
##
## For each element of @var{x}, compute the cumulative distribution function
## (CDF) at @var{x} of the uniform distribution on the interval [@var{a},
## @var{b}].  The size of @var{p} is the common size of the input arguments.
## A scalar input functions as a constant matrix of the same size as the other
## inputs.
##
## Default values are @var{a} = 0, @var{b} = 1.
##
## @code{[@dots{}] = unifcdf (@dots{}, "upper")} computes the upper tail
## probability of the lognormal distribution.
##
## @seealso{unifinv, unifpdf, unifrnd, unifstat}
## @end deftypefn

function p = unifcdf (x, varargin)

  ## Check for valid number of input arguments
  if (nargin < 1 || nargin > 4)
    error ("unifcdf: invalid number of input arguments.");
  endif

  ## Check for "upper" flag
  if (nargin > 1 && strcmpi (varargin{end}, "upper"))
    uflag = true;
    varargin(end) = [];
  elseif (nargin > 1 && ischar (varargin{end}) && ...
                      ! strcmpi (varargin{end}, "upper"))
    error ("unifcdf: invalid argument for upper tail.");
  elseif (nargin > 3 && ! strcmpi (varargin{end}, "upper"))
    error ("unifcdf: invalid argument for upper tail.");
  else
    uflag = false;
  endif

  ## Get extra arguments (if they exist) or add defaults
  if (numel (varargin) > 0)
    a = varargin{1};
  else
    a = 0;
  endif
  if (numel (varargin) > 1)
    b = varargin{2};
  else
    b = 1;
  endif

  ## Check for common size of X, A, and B
  if (! isscalar (x) || ! isscalar (a) || ! isscalar (b))
    [retval, x, a, b] = common_size (x, a, b);
    if (retval > 0)
      error ("unifcdf: X, A, and B must be of common size or scalars.");
    endif
  endif

  ## Check for X and SIGMA being reals
  if (iscomplex (x) || iscomplex (a) || iscomplex (b))
    error ("unifcdf: X, A, and B must not be complex.");
  endif

  ## Check for appropriate class
  if (isa (x, "single") || isa (a, "single") || isa (b, "single"))
    p = zeros (size (x), "single");
  else
    p = zeros (size (x));
  endif

  ## Calculate Rayleigh CDF for valid parameter and data range
  k = find(x > a & x < b & a < b);
  if (uflag)
    p(x <= a & a < b) = 1;
    p(x >= b & a < b) = 0;
    if any(k)
      p(k) = (b(k)- x(k)) ./ (b(k) - a(k));
    endif
  else
    p(x <= a & a < b) = 0;
    p(x >= b & a < b) = 1;
    if any(k)
      p(k) = (x(k) - a(k)) ./ (b(k) - a(k));
    endif
  endif

  ## Continue argument check
  p(a >= b) = NaN;
  p(isnan(x) | isnan(a) | isnan(b)) = NaN;

endfunction


%!shared x,y
%! x = [-1 0 0.5 1 2] + 1;
%! y = [0 0 0.5 1 1];
%!assert (unifcdf (x, ones (1,5), 2*ones (1,5)), y)
%!assert (unifcdf (x, ones (1,5), 2*ones (1,5), "upper"), 1 - y)
%!assert (unifcdf (x, 1, 2*ones (1,5)), y)
%!assert (unifcdf (x, 1, 2*ones (1,5), "upper"), 1 - y)
%!assert (unifcdf (x, ones (1,5), 2), y)
%!assert (unifcdf (x, ones (1,5), 2, "upper"), 1 - y)
%!assert (unifcdf (x, [2 1 NaN 1 1], 2), [NaN 0 NaN 1 1])
%!assert (unifcdf (x, [2 1 NaN 1 1], 2, "upper"), 1 - [NaN 0 NaN 1 1])
%!assert (unifcdf (x, 1, 2*[0 1 NaN 1 1]), [NaN 0 NaN 1 1])
%!assert (unifcdf (x, 1, 2*[0 1 NaN 1 1], "upper"), 1 - [NaN 0 NaN 1 1])
%!assert (unifcdf ([x(1:2) NaN x(4:5)], 1, 2), [y(1:2) NaN y(4:5)])
%!assert (unifcdf ([x(1:2) NaN x(4:5)], 1, 2, "upper"), 1 - [y(1:2) NaN y(4:5)])

## Test class of input preserved
%!assert (unifcdf ([x, NaN], 1, 2), [y, NaN])
%!assert (unifcdf (single ([x, NaN]), 1, 2), single ([y, NaN]))
%!assert (unifcdf ([x, NaN], single (1), 2), single ([y, NaN]))
%!assert (unifcdf ([x, NaN], 1, single (2)), single ([y, NaN]))

## Test input validation
%!error unifcdf ()
%!error unifcdf (1, 2, 3, 4)
%!error unifcdf (1, 2, 3,"upper", 4)
%!error unifcdf (ones (3), ones (2), ones (2))
%!error unifcdf (ones (2), ones (3), ones (2))
%!error unifcdf (ones (2), ones (2), ones (3))
%!error unifcdf (i, 2, 2)
%!error unifcdf (2, i, 2)
%!error unifcdf (2, 2, i)
