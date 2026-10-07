## Copyright (C) 2006 Frederick (Rick) A Niles <niles@rickniles.com>
## Copyright (C) 2023 Andreas Bertsatos <abertsatos@biol.uoa.gr>
##
## This file is part of the statistics package for GNU Octave.
##
## This program is free software; you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation; either version 3 of the License, or (at your option) any later
## version.
##
## This program is distributed in the hope that it will be useful, but WITHOUT
## ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
## FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
## more details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn  {statistics} {@var{y} =} jsupdf (@var{x})
## @deftypefnx {statistics} {@var{y} =} jsupdf (@var{x}, @var{alpha1})
## @deftypefnx {statistics} {@var{y} =} jsupdf (@var{x}, @var{alpha1}, @var{alpha2})
##
## Johnson SU probability density function (PDF).
##
## For each element of @var{x}, compute the probability density function (PDF)
## at @var{x} of the Johnson SU distribution with shape parameters @var{alpha1}
## and @var{alpha2}, which is
## @code{@var{alpha2} * normpdf (@var{alpha1} + @var{alpha2} * asinh (@var{x}))
## / sqrt (@var{x}^2 + 1)}.  The size of @var{y} is the common size of the input
## arguments @var{x}, @var{alpha1}, and @var{alpha2}.  A scalar input functions
## as a constant matrix of the same size as the other inputs.
##
## Default values are @var{alpha1} = 1, @var{alpha2} = 1.  @var{alpha2} must be
## positive; where it is not, the result is @qcode{NaN}.
##
## Input arguments must be @qcode{double} or @qcode{single}; integer, logical,
## and character arrays are rejected.  MATLAB accepts a character array and
## evaluates it at the character codes, which Octave deliberately does not,
## since a character array is an integer type and integers are refused too.
##
## @seealso{jsucdf}
## @end deftypefn

function y = jsupdf (x, alpha1, alpha2)

  if (nargin < 1 || nargin > 3)
    print_usage;
  endif

  if (nargin == 1)
    alpha1 = 1;
    alpha2 = 1;
  elseif (nargin == 2)
    alpha2 = 1;
  endif

  ## Check for X, ALPHA1, and ALPHA2 being double or single
  if (! (isfloat (x) && isfloat (alpha1) && isfloat (alpha2)))
    error ("jsupdf: X, ALPHA1, and ALPHA2 must be double or single.");
  endif

  if (! isscalar (x) || ! isscalar (alpha1) || ! isscalar (alpha2))
    [retval, x, alpha1, alpha2] = common_size (x, alpha1, alpha2);
    if (retval > 0)
      error (strcat ("jsupdf: X, ALPHA1, and ALPHA2 must be of common", ...
                     " size or scalars."));
    endif
  endif

  y = alpha2 ./ hypot (x, 1) .* normpdf (alpha1 + alpha2 .* asinh (x));
  y((alpha2 <= 0) & true (size (y))) = NaN;

endfunction

%!assert_equal (jsupdf (0), normpdf (1), -1e-15)
%!assert_equal (jsupdf (1, 0.5, 2), ...
%!              sqrt (2) * normpdf (0.5 + 2 * log (1 + sqrt (2))), -1e-14)
%!assert_equal (jsupdf ([-Inf, NaN, Inf]), [0, NaN, 0])
%!assert_equal (integral (@(x) jsupdf (x, 0.5, 2), -Inf, 1), ...
%!              jsucdf (1, 0.5, 2), -1e-9)
%!assert_equal (jsupdf (single (0)), single (normpdf (1)), -eps ('single'))
%!assert_equal (jsupdf (1, 0, 0), NaN)
%!assert_equal (jsupdf ([1, 2], 1, -1), [NaN, NaN])
%!assert_equal (jsupdf (1, 1, [1, 0]), [jsupdf(1), NaN])

%!error<jsupdf: X, ALPHA1, and ALPHA2 must be double or single.> jsupdf (int32 (2), 1, 1)
%!error<jsupdf: X, ALPHA1, and ALPHA2 must be double or single.> jsupdf (true, 1, 1)
%!error<jsupdf: X, ALPHA1, and ALPHA2 must be double or single.> jsupdf ('a', 1, 1)
%!error jsupdf ()
%!error jsupdf (1, 2, 3, 4)
%!error<jsupdf: X, ALPHA1, and ALPHA2 must be of common size or scalars.> ...
%! jsupdf (1, ones (2), ones (3))

