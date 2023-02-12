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
## @deftypefn  {statistics} {@var{p} =} betacdf (@var{x}, @var{a}, @var{b})
## @deftypefnx {statistics} {@var{p} =} betacdf (@var{x}, @var{a}, @var{b}, "upper")
##
## Beta cumulative distribution function (CDF).
##
## For each element of @var{x}, compute the cumulative distribution function
## at @var{x} of the Beta distribution with parameters @var{a} and @var{b}.  The
## size of @var{p} is the common size of @var{x}, @var{a} and @var{b}.  A scalar
## input functions as a constant matrix of the same size as the other inputs.
##
## @code{@var{p} = betacdf (@var{x}, @var{a}, @var{b}, "upper")} computes the
## upper tail probability of the Beta distribution with parameters @var{a} and
## @var{b} at the values in @var{x}.
##
## @seealso{betainv, betapdf, betarnd, betastat}
## @end deftypefn

function p = betacdf (x, a, b, varargin)

  ## Check for valid number of input arguments
  if (nargin < 3 || nargin > 4)
    error ("betacdf: invalid number of input arguments.");
  endif

  ## Check for valid "upper" flag
  if (nargin > 3 && ! strcmpi (varargin{1}, "upper"))
    error ("betacdf: invalid argument for upper tail.");
  endif

  ## Check for common size of X, A and B
  if (! isscalar (x) || ! isscalar (a) || ! isscalar (b))
    [err, x, a, b] = common_size (x, a, b);
    if (err > 0)
      error ("betacdf: X, A, and B must be of common size or scalars.");
    endif
  endif

  ## Check for X,A and B being reals
  if (iscomplex (x) || iscomplex (a) || iscomplex (b))
    error ("betacdf: X, A, and B must not be complex.");
  endif

  ## Check for class type
  if (isa (x, "single") || isa (a, "single") || isa (b, "single"))
    is_type = "single";
  else
    is_type = "double";
  endif

  ## Find valid values in parameters and data
  okPARAM = (0 < a & a < Inf) & (0 < b & b < Inf);
  okDATA = (okPARAM & (0 <= x & x <= 1));
  all_OK = all (okDATA(:));

  ## Force NaNs for out of range parameters.
  ## Fill in edges cases when X is outside 0 or 1.
  if (! all_OK)
    p = NaN (size (okDATA), is_type);
    if (nargin > 3 && strcmpi (varargin{1}, "upper"))
      p(okPARAM & x <= 0) = 1;
      p(okPARAM & x >= 1) = 0;
    else
      p(okPARAM & x < 0) = 0;
      p(okPARAM & x > 1) = 1;
    endif
    ## Remove the out of range/edge cases. Return, if there's nothing left.
    if (any (okDATA(:)))
      if (numel (x) > 1)
        x = x(okDATA);
      endif
      if (numel (a) > 1)
        a = a(okDATA);
      endif
      if (numel (b) > 1)
        b = b(okDATA);
      endif
    else
      return;
    endif
  endif

  ## Call betainc for the actual work
  pk = betainc (x, a, b, varargin{:});

  ## Relocate the values to the correct places if necessary.
  if all_OK
    p = pk;
  else
    p(okDATA) = pk;
  endif

endfunction

## Test output
%!shared x, y, x1, x2
%! x = [-1 0 0.5 1 2];
%! y = [0 0 0.75 1 1];
%!assert (betacdf (x, ones (1,5), 2*ones (1,5)), y)
%!assert (betacdf (x, 1, 2*ones (1,5)), y)
%!assert (betacdf (x, ones (1,5), 2), y)
%!assert (betacdf (x, [0 1 NaN 1 1], 2), [NaN 0 NaN 1 1])
%!assert (betacdf (x, 1, 2*[0 1 NaN 1 1]), [NaN 0 NaN 1 1])
%!assert (betacdf ([x(1:2) NaN x(4:5)], 1, 2), [y(1:2) NaN y(4:5)])
%! x1 = [0.1:0.2:0.9];
%!assert (betacdf (x1, 2, 2), [0.028, 0.216, 0.5, 0.784, 0.972], 1e-14);
%!assert (betacdf (x1, 2, 2, "upper"), 1 - [0.028, 0.216, 0.5, 0.784, 0.972],...
%!        1e-14);
%! x2 = [1, 2, 3];
%!assert (betacdf (0.5, x2, x2), [0.5, 0.5, 0.5], 1e-14);

## Test class of input preserved
%!assert (betacdf ([x, NaN], 1, 2), [y, NaN])
%!assert (betacdf (single ([x, NaN]), 1, 2), single ([y, NaN]))
%!assert (betacdf ([x, NaN], single (1), 2), single ([y, NaN]))
%!assert (betacdf ([x, NaN], 1, single (2)), single ([y, NaN]))

## Test input validation
%!error<betacdf: invalid number of input arguments.> betacdf ()
%!error<betacdf: invalid number of input arguments.> betacdf (1)
%!error<betacdf: invalid number of input arguments.> betacdf (1,2)
%!error<betacdf: invalid number of input arguments.> betacdf (1,2,3,4,5)
%!error<betacdf: invalid argument for upper tail.> betacdf (1,2,3,"tail")
%!error<betacdf: X, A, and B must be of common size or scalars.> ...
%! betacdf (ones (3), ones (2), ones (2))
%!error<betacdf: X, A, and B must be of common size or scalars.> ...
%! betacdf (ones (2), ones (3), ones (2))
%!error<betacdf: X, A, and B must be of common size or scalars.> ...
%! betacdf (ones (2), ones (2), ones (3))
%!error<betacdf: X, A, and B must not be complex.> betacdf (i, 2, 2)
%!error<betacdf: X, A, and B must not be complex.> betacdf (2, i, 2)
%!error<betacdf: X, A, and B must not be complex.> betacdf (2, 2, i)
