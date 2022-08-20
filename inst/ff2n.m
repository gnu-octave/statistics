## Copyright (C) 2022 Andreas Bertsatos <abertsatos@biol.uoa.gr>
## based on public domain work by Paul Kienzle <pkienzle@users.sf.net>
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
## @deftypefn {Function File} ff2n (@var{n})
## Full-factor design with n binary terms.
##
## @seealso {fullfact}
## @end deftypefn

function A = ff2n (n)
  if (nargin != 1)
    error ("ff2n: wrong number of input arguments.");
  endif
  if (floor (n) != n || numel (n) != 1 || n < 1 ...
                     || ! isfinite (n) || ! isreal (n))
    error ("ff2n: @var{N} must be a positive integer scalar.");
  endif
  A = fullfact (2 * ones (1, n)) - 1;
endfunction

%!error ff2n ();
%!error ff2n (2, 5);
%!error ff2n (2.5);
%!error ff2n (0);
%!error ff2n (-3);
%!error ff2n (3+2i);
%!error ff2n (Inf);
%!error ff2n (NaN);
%!test
%! A = ff2n (3);
%! assert (A, fullfact (3));
%!test
%! A = ff2n (8);
%! assert (A, fullfact (8));

