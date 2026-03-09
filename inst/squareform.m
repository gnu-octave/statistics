## Copyright (C) 2015 Carnë Draug <carandraug@octave.org>
## Copyright (C) 2026 Avanish Salunke <avanishsalunke16@gmail.com>
##
## This file is part of the statistics package for GNU Octave.
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
## @deftypefn  {statistics} {@var{z} =} squareform (@var{y})
## @deftypefnx {statistics} {@var{y} =} squareform (@var{z})
## @deftypefnx {statistics} {@var{z} =} squareform (@var{y}, @qcode{"tovector"})
## @deftypefnx {statistics} {@var{y} =} squareform (@var{z}, @qcode{"tomatrix"})
##
## Interchange between distance matrix and distance vector formats.
##
## Converts between a hollow (diagonal filled with zeros), square, and
## symmetric matrix and a vector of the lower triangular part.
##
## Its target application is the conversion of the vector returned by
## @code{pdist} into a distance matrix.  It performs the opposite operation
## if input is a matrix.
##
## If @var{x} is a numeric or logical vector, its number of elements must fit into the
## triangular part of a matrix (main diagonal excluded).  In other words,
## @code{numel (@var{x}) = @var{n} * (@var{n} - 1) / 2} for some integer
## @var{n}.  The resulting matrix will be @var{n} by @var{n}.
##
## If @var{x} is a numeric or logical distance matrix, it must be square and the diagonal entries
## of @var{x} must all be zeros.  If @var{x} is not symmetric, only the lower
## triangular part is used.
##
## The second argument is used to specify the output type in case there
## is a single element.  It will default to @qcode{"tomatrix"} otherwise.
## The string is case-insensitive and partial matches are accepted as
## long as they are not ambiguous (e.g., @qcode{"tom"} or @qcode{"tov"}).
##
## @seealso{pdist}
## @end deftypefn

function y = squareform (x, method)

  if (nargin < 1 || nargin > 2)
    print_usage ();
  elseif (! (isnumeric (x) || islogical (x)) || ! ismatrix (x))
    error ("squareform: Y or Z must be a numeric or logical matrix or vector.");
  endif

  if (nargin == 1)
    ## This is ambiguous when numel (x) == 1, but that's the whole reason
    ## why the "method" option exists.
    if (isvector (x))
      method = "tomatrix";
    else
      method = "tovector";
    endif
  elseif (ischar (method) || isstring (method))
    method = tolower (method);
    if (length (method) >= 3)
      if (strncmp (method, "tovector", length (method)))
        method = "tovector";
      elseif (strncmp (method, "tomatrix", length (method)))
        method = "tomatrix";
      endif
    endif
  endif

  switch (method)
    case "tovector"
      if (! issquare (x))
        error ("squareform: Z is not a square matrix.");
      elseif (any (diag (x) != 0))
        error ("squareform: Z is not a hollow matrix.");
      endif

      y = vec (tril (x, -1, "pack"), 2);

    case "tomatrix"
      ## the dimensions of y are the solution to the quadratic formula for:
      ## length (x) = (sy - 1) * (sy / 2)
      sy = (1 + sqrt (1 + 8 * numel (x))) / 2;
      if (fix (sy) != sy)
        error ("squareform: the numel of Y cannot form a square matrix.");
      endif

      y = zeros (sy, class (x));
      y(tril (true (sy), -1)) = x;  # fill lower triangular part
      if (islogical (x))
        y = y | y.'; 
      else
        y += y.'; # and then the upper triangular part
      endif

    otherwise
      error ("squareform: invalid METHOD '%s'.", method);
  endswitch

endfunction

%!test
%! v = 1:6;
%! m = [0 1 2 3; 1 0 4 5; 2 4 0 6; 3 5 6 0];
%! assert (squareform (v), m);
%! assert (squareform (m), v);
%! assert (squareform (squareform (v)), v);

## test the row and column vectors equally
%!test
%! v = 1:6;
%! m = [0 1 2 3; 1 0 4 5; 2 4 0 6; 3 5 6 0];
%! assert (squareform (v'), m);

## handle 1 element input properly
%!test
%! assert (squareform (1), [0 1; 1 0]);
%! assert (squareform (1, "tomatrix"), [0 1; 1 0]);
%! assert (squareform (0, "tovector"), zeros (1, 0));

## test logical inputs 
%!test
%! v = [true, false, true];
%! m = [false, true, false; true, false, true; false, true, false];
%! assert (squareform (v), m);
%! assert (squareform (m), v);
%! assert (class (squareform (v)), "logical");
%! assert (class (squareform (m)), "logical");

## confirm that it respects input class
%!test
%! v = 1:6;
%! m = [0 1 2 3; 1 0 4 5; 2 4 0 6; 3 5 6 0];
%! classes = {@logical, @single, @double, @uint8, @uint32, @uint64};
%! for i = 1:numel (classes)
%!   f = classes{i};
%!   assert (squareform (f (v)), f (m));
%!   assert (squareform (f (m)), f (v));
%! endfor

## test partial string matching and case insensitivity
%!test
%! v = 1:6;
%! m = [0 1 2 3; 1 0 4 5; 2 4 0 6; 3 5 6 0];
%! assert (squareform (v, "tom"), m);
%! assert (squareform (m, "tov"), v);
%! assert (squareform (v, "TOMATRIX"), m);

## input validations
%!error <squareform: Y or Z must be a numeric or logical matrix or vector.> squareform ("string")
%!error <squareform: Y or Z must be a numeric or logical matrix or vector.> squareform ({1, 2, 3})
%!error <squareform: Z is not a square matrix.> squareform ([1, 2, 3; 4, 5, 6], "tovector")
%!error <squareform: Z is not a hollow matrix.> squareform (eye (3), "tovector")
%!error <squareform: the numel of Y cannot form a square matrix.> squareform ([1, 2, 3, 4], "tomatrix")
%!error <squareform: invalid METHOD 'to'> squareform ([1, 2, 3], "to")
%!error <squareform: invalid METHOD 'invalid'> squareform ([1, 2, 3], "invalid")

