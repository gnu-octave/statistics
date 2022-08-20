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
## @deftypefn {Function File} fullfact (@var{N})
## Full factorial design.
##
## If @var{N} is a scalar, return the full factorial design with @var{N} binary
## choices, 0 and 1.
##
## If @var{N} is a vector, return the full factorial design with ordinal choices
## 1 through @var{n_i} for each factor @var{i}.
##
## Values in @var{N} must be positive integers.
##
## @seealso {ff2n}
## @end deftypefn

function A = fullfact(n)
  if (nargin != 1)
    error ("fullfact: wrong number of input arguments.");
  endif
  if length(n) == 1
    if (floor (n) != n || n < 1 || ! isfinite (n) || ! isreal (n))
      error ("fullfact: @var{N} must be a positive integer.");
    endif
    ## Combinatorial design with n binary choices
    A = fullfact(2*ones(1,n))-1;
  else
    if (any (floor (n) != n) || any (n < 1) || any (! isfinite (n)) ...
                             || any (! isreal (n)))
      error ("fullfact: values in @var{N} must be positive integers.");
    endif
    ## Combinatorial design with n(i) ordinal choices
    A = [1:n(end)]';
    for i=length(n)-1:-1:1
      A = [kron([1:n(i)]',ones(rows(A),1)), repmat(A,n(i),1)];
    end
  end
endfunction

%!demo
%! ## Full factorial design with 3 binary variables
%! fullfact (3)

%!demo
%! ## Full factorial design with 3 ordinal variables
%! fullfact ([2, 3, 4])

%!error fullfact ();
%!error fullfact (2, 5);
%!error fullfact (2.5);
%!error fullfact (0);
%!error fullfact (-3);
%!error fullfact (3+2i);
%!error fullfact (Inf);
%!error fullfact (NaN);
%!error fullfact ([1, 2, -3]);
%!error fullfact ([0, 1, 2]);
%!error fullfact ([1, 2, NaN]);
%!error fullfact ([1, 2, Inf]);
%!test
%! A = fullfact (2);
%! assert (A, [0, 0; 0, 1; 1, 0; 1, 1]);
%!test
%! A = fullfact ([1, 2]);
%! assert (A, [1, 1; 1, 2]);
%!test
%! A = fullfact ([1, 2, 4]);
%! A_out = [1, 1, 1; 1, 1, 2; 1, 1, 3; 1, 1, 4; ...
%!          1, 2, 1; 1, 2, 2; 1, 2, 3; 1, 2, 4];
%! assert (A, A_out);
