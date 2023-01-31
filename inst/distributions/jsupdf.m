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
## FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
## details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn  {statistics} @var{y} = jsupdf (@var{x})
## @deftypefnx {statistics} @var{y} = jsupdf (@var{x}, @var{alpha1})
## @deftypefnx {statistics} @var{y} = jsupdf (@var{x}, @var{alpha1}, @var{alpha2})
##
## Johnson SU probability density function (PDF).
##
## For each element of @var{x}, compute the probability density function (PDF)
## at @var{x} of the Johnson SU distribution with shape parameters @var{alpha1}
## and @var{alpha2}.  The size of @var{p} is the common size of the input
## arguments @var{x}, @var{alpha1}, and @var{alpha2}.  A scalar input functions
## as a constant matrix of the same size as the other
##
## Default values are @var{alpha1} = 1, @var{alpha2} = 1.
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

  if (! isscalar (x) || ! isscalar (alpha1) || ! isscalar(alpha2))
    [retval, x, alpha1, alpha2] = common_size (x, alpha1, alpha2);
    if (retval > 0)
      error (strcat (["jsupdf: X, ALPHA1, and ALPHA2 must be of common"], ...
                     [" size or scalars."]));
    endif
  endif

  one = ones (size (x));
  sr = sqrt (x .* x + one);
  y = (alpha2 ./ sr) .* ...
      stdnormal_pdf (alpha1 .* one + alpha2 .* log (x + sr));

endfunction

%!error jsupdf ()
%!error jsupdf (1, 2, 3, 4)
%!error<jsupdf: X, ALPHA1, and ALPHA2 must be of common size or scalars.> ...
%! jsupdf (1, ones (2), ones (3))

