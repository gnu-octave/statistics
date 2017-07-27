## Copyright (C) 2015 Carnë Draug <carandraug@octave.org>
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn  {Function File} {@var{z} =} squareform (@var{y})
## @deftypefnx {Function File} {@var{y} =} squareform (@var{z})
## @deftypefnx {Function File} {@var{z} =} squareform (@var{y}, @qcode{"tovector"})
## @deftypefnx {Function File} {@var{y} =} squareform (@var{z}, @qcode{"tomatrix"})
## Interchange between distance matrix and distance vector formats.
##
## Converts between an hollow (diagonal filled with zeros), square, and
## symmetric matrix and a vector with of the lower triangular part.
##
## Its target application is the conversion of the vector returned by
## @code{pdist} into a distance matrix.  It performs the opposite operation
## if input is a matrix.
##
## If @var{x} is a vector, its number of elements must fit into the
## triangular part of a matrix (main diagonal excluded).  In other words,
## @code{numel (@var{x}) = @var{n} * (@var{n} - 1) / 2} for some integer
## @var{n}.  The resulting matrix will be @var{n} by @var{n}.
##
## If @var{x} is a distance matrix, it must be square and the diagonal entries
## of @var{x} must all be zeros.  @code{squareform} will generate a warning if
## @var{x} is not symmetric.
##
## The second argument is used to specify the output type in case there
## is a single element.  It will defaults to @qcode{"tomatrix"} otherwise.
##
## @seealso{pdist}
## @end deftypefn

## Author: Carnë Draug <carandraug@octave.org>

function y = squareform (x, method)

  if (nargin < 1 || nargin > 2)
    print_usage ();
  elseif (! isnumeric (x) || ! ismatrix (x))
    error ("squareform: Y or Z must be a numeric matrix or vector");
  endif

  if (nargin == 1)
    ## This is ambiguous when numel (x) == 1, but that's the whole reason
    ## why the "method" option exists.
    if (isvector (x))
      method = "tomatrix";
    else
      method = "tovector";
    endif
  endif

  switch (tolower (method))
    case "tovector"
      if (! issquare (x))
        error ("squareform: Z is not a square matrix");
      elseif (any (diag (x) != 0))
        error ("squareform: Z is not a hollow matrix, i.e., with diagonal entries all zero");
      elseif (! issymmetric(x))
        warning ("squareform:symmetric",
                 "squareform: Z is not a symmetric matrix");
      endif

      y = vec (tril (x, -1, "pack"), 2);

    case "tomatrix"
      ## the dimensions of y are the solution to the quadratic formula for:
      ## length (x) = (sy - 1) * (sy / 2)
      sy = (1 + sqrt (1 + 8 * numel (x))) / 2;
      if (fix (sy) != sy)
        error ("squareform: the numel of Y cannot form a square matrix");
      endif

      y = zeros (sy, class (x));
      y(tril (true (sy), -1)) = x;  # fill lower triangular part
      y += y.'; # and then the upper triangular part

    otherwise
      error ("squareform: invalid METHOD '%s'", method);
  endswitch

endfunction

%!shared v, m
%! v = 1:6;
%! m = [0 1 2 3;1 0 4 5;2 4 0 6;3 5 6 0];

## make sure that it can go both directions automatically
%!assert (squareform (v), m)
%!assert (squareform (squareform (v)), v)
%!assert (squareform (m), v)

## treat row and column vectors equally
%!assert (squareform (v'), m)

## handle 1 element input properly
%!assert (squareform (1), [0 1;1 0])
%!assert (squareform (1, "tomatrix"), [0 1; 1 0])
%!assert (squareform (0, "tovector"), zeros (1, 0))

%!warning <not a symmetric matrix> squareform ([0 1 2; 3 0 4; 5 6 0]);

## confirm that it respects input class
%!test
%! for c = {@single, @double, @uint8, @uint32, @uint64}
%!   f = c{1};
%!   assert (squareform (f (v)), f (m))
%!   assert (squareform (f (m)), f (v))
%! endfor

