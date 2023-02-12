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
## @deftypefn  {statistics} {@var{p} =} jsucdf (@var{x})
## @deftypefnx {statistics} {@var{p} =} jsucdf (@var{x}, @var{alpha1})
## @deftypefnx {statistics} {@var{p} =} jsucdf (@var{x}, @var{alpha1}, @var{alpha2})
##
## Johnson SU cumulative distribution function (CDF).
##
## For each element of @var{x}, return the cumulative distribution functions
## (CDF) at @var{x} of the Johnson SU distribution with shape parameters
## @var{alpha1} and @var{alpha2}.  The size of @var{p} is the common size of the
## input arguments @var{x}, @var{alpha1}, and @var{alpha2}.  A scalar input
## functions as a constant matrix of the same size as the other
##
## Default values are @var{alpha1} = 1, @var{alpha2} = 1.
##
## @seealso{jsupdf}
## @end deftypefn

function p = jsucdf (x, alpha1, alpha2)

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
      error (strcat (["jsucdf: X, ALPHA1, and ALPHA2 must be of common"], ...
                     [" size or scalars."]));
    endif
  endif

  one = ones (size (x));
  p = stdnormal_cdf (alpha1 .* one + alpha2 .* log (x + sqrt (x .* x + one)));

endfunction

%!error jsucdf ()
%!error jsucdf (1, 2, 3, 4)
%!error<jsucdf: X, ALPHA1, and ALPHA2 must be of common size or scalars.> ...
%! jsucdf (1, ones (2), ones (3))
