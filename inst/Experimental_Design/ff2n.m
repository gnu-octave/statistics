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
## FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
## more details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn  {statistics} {@var{dFF2} =} ff2n (@var{n})
##
## Two-level full factorial design.
##
## @code{@var{dFF2} = ff2n (@var{n})} gives factor settings dFF2 for a two-level
## full factorial design with n factors.  @var{dFF2} is m-by-n, where m is the
## number of treatments in the full-factorial design.  Each row of @var{dFF2}
## corresponds to a single treatment.  Each column contains the settings for a
## single factor, with values of 0 and 1 for the two levels.
##
## @var{n} must be a non-negative integer scalar.  @code{ff2n (0)} returns a
## 1-by-0 matrix, a single treatment with no factors.  An empty @var{n} is an
## error, whereas MATLAB returns an empty matrix.
##
## @seealso{fullfact}
## @end deftypefn

function A = ff2n (n)
  if (nargin != 1)
    error ("ff2n: wrong number of input arguments.");
  endif
  if (! (isscalar (n) && (isnumeric (n) || islogical (n)) && isreal (n)
         && isfinite (n) && n >= 0 && fix (n) == n))
    error ("ff2n: N must be a non-negative integer scalar.");
  endif
  if (n == 0)
    A = zeros (1, 0);
  else
    A = flip (fullfact (2 * ones (1, n)), 2) - 1;
  endif
endfunction

%!error ff2n ();
%!error ff2n (2, 5);
%!error <ff2n: N must be a non-negative integer scalar.> ff2n ([])
%!error <ff2n: N must be a non-negative integer scalar.> ff2n ([1, 2])
%!error <ff2n: N must be a non-negative integer scalar.> ff2n ('a')
%!error <ff2n: N must be a non-negative integer scalar.> ff2n (2.5)
%!error <ff2n: N must be a non-negative integer scalar.> ff2n (-3)
%!error <ff2n: N must be a non-negative integer scalar.> ff2n (3+2i)
%!error <ff2n: N must be a non-negative integer scalar.> ff2n (Inf)
%!error <ff2n: N must be a non-negative integer scalar.> ff2n (NaN)
%!assert_equal (ff2n (0), zeros (1, 0))
%!assert_equal (ff2n (true), [0; 1])
%!assert_equal (ff2n (int8 (2)), [0, 0; 0, 1; 1, 0; 1, 1])
%!test
%! A = ff2n (3);
%! assert_equal (A, [0, 0, 0; 0, 0, 1; 0, 1, 0; 0, 1, 1; ...
%!             1, 0, 0; 1, 0, 1; 1, 1, 0; 1, 1, 1]);
%!test
%! A = ff2n (2);
%! assert_equal (A, [0, 0; 0, 1; 1, 0; 1, 1]);

