## Copyright (C) 2022-2025 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {statistics} {@var{A} =} fullfact (@var{levels})
##
## Full factorial design.
##
## @code{@var{A} =} fullfact (@var{levels}) returns a numeric matrix @var{A}
## with the treatments of a full factorial design specified by @var{levels},
## which must be a numeric vector of real positive integer values with each
## value specifying the number of levels of each individual factor.
##
## Each row of @var{A} corresponds to a single treatment and each column to a
## single factor.  For binary full factorial design, use @code{ff2n}.
##
## @seealso {ff2n}
## @end deftypefn

function A = fullfact (levels)
  if (nargin != 1)
    error ("fullfact: one input argument is required.");
  endif
  if (! (isvector (levels) && isnumeric (levels) && isfinite (levels)))
    error ("fullfact: input argument must be a finite real numeric vector.");
  endif
  if (any (fix (levels) != levels) || any (levels < 1) || ! all (isreal (levels)))
    error ("fullfact: factor levels must be real positive integers.");
  endif
  rows = prod (levels);
  cols = numel (levels);
  A = zeros (rows, cols);
  n_seqs = rows;
  for i = 1:cols
    factor = [1:levels(i)]';
    n_reps = rows / n_seqs;
    factor = repelem (factor, n_reps, 1);
    n_seqs = n_seqs / levels(i);
    factor = repmat (factor, n_seqs, 1);
    A(:,i) = factor;
  endfor
endfunction

%!demo
%! ## Full factorial design with 3 ordinal variables
%! fullfact ([2, 3, 4])

%!error<fullfact: one input argument is required.> fullfact ();
%!error<fullfact: input argument must be a finite real numeric vector.> ...
%! fullfact (Inf);
%!error<fullfact: input argument must be a finite real numeric vector.> ...
%! fullfact (NaN);
%!error<fullfact: input argument must be a finite real numeric vector.> ...
%! fullfact (ones (2));
%!error<fullfact: input argument must be a finite real numeric vector.> ...
%! fullfact ([1, 2, NaN]);
%!error<fullfact: input argument must be a finite real numeric vector.> ...
%! fullfact ([1, 2, Inf]);
%!error<fullfact: factor levels must be real positive integers.> fullfact (2.5);
%!error<fullfact: factor levels must be real positive integers.> fullfact (0);
%!error<fullfact: factor levels must be real positive integers.> fullfact (-3);
%!error<fullfact: factor levels must be real positive integers.> fullfact (3+2i);
%!error<fullfact: factor levels must be real positive integers.> fullfact ([1, 2, -3]);
%!error<fullfact: factor levels must be real positive integers.> fullfact ([0, 1, 2]);
%!test
%! A = fullfact (1);
%! assert (A, 1);
%!test
%! A = fullfact (2);
%! assert (A, [1; 2]);
%!test
%!test
%! A = fullfact (3);
%! assert (A, [1; 2; 3]);
%!test
%! A = fullfact ([1, 2, 4]);
%! A_out = [1, 1, 1; 1, 2, 1; 1, 1, 2; 1, 2, 2; ...
%!          1, 1, 3; 1, 2, 3; 1, 1, 4; 1, 2, 4];
%! assert (A, A_out);
%!test
%! A = fullfact ([2, 2]);
%! assert (A, [1, 1; 2, 1; 1, 2; 2, 2]);
%!test
%! A = fullfact ([2, 2, 4]);
%! A_out = [1, 1, 1; 2, 1, 1; 1, 2, 1; 2, 2, 1; ...
%!          1, 1, 2; 2, 1, 2; 1, 2, 2; 2, 2, 2; ...
%!          1, 1, 3; 2, 1, 3; 1, 2, 3; 2, 2, 3; ...
%!          1, 1, 4; 2, 1, 4; 1, 2, 4; 2, 2, 4];
%! assert (A, A_out);
%!test
%! A = fullfact ([3, 2, 4]);
%! A_out = [1, 1, 1; 2, 1, 1; 3, 1, 1; 1, 2, 1; 2, 2, 1; 3, 2, 1; ...
%!          1, 1, 2; 2, 1, 2; 3, 1, 2; 1, 2, 2; 2, 2, 2; 3, 2, 2; ...
%!          1, 1, 3; 2, 1, 3; 3, 1, 3; 1, 2, 3; 2, 2, 3; 3, 2, 3; ...
%!          1, 1, 4; 2, 1, 4; 3, 1, 4; 1, 2, 4; 2, 2, 4; 3, 2, 4];
%! assert (A, A_out);
%!test
%! A = fullfact ([4, 2]);
%! assert (A, [1, 1; 2, 1; 3, 1; 4, 1; 1, 2; 2, 2; 3, 2; 4, 2]);
