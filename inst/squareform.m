## Copyright (C) 2014 JD Walsh. Based on code (C) 2006, 2008 Bill Denney
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
## @deftypefn  {Function File} {@var{y} =} squareform (@var{z})
## @deftypefnx {Function File} {@var{z} =} squareform (@var{y}, @qcode{"tovector"})
## @deftypefnx {Function File} {@var{y} =} squareform (@var{z}, @qcode{"tomatrix"})
## Converts a vector from the @code{pdist} function into a distance matrix or
## a distance matrix back to vector form.
##
## If @var{x} is a vector, it must have
## @code{length(@var{x}) = @var{n} * (@var{n} - 1) / 2} for some integer
## @var{n}. The resulting matrix will be @var{n} by @var{n}.
##
## If @var{x} is a distance matrix, it must be square and the diagonal entries
## of @var{x} must all be zeros. @code{squareform} will generate a warning if
## @var{x} is not symmetric.
##
## The second argument is used to specify the output type in case there
## is a single element.  It will defaults to @qcode{"tomatrix"} otherwise.
##
## @seealso{pdist}
## @end deftypefn

## Author: JD Walsh <walsh@math.gatech.edu>
## Created: 2014-11-09
## Description: Convert distance matrix from vector to square form and back
## Keywords: distance format

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
      if (~issquare (x))
        usage ('squareform: x is not a square matrix');
      elseif (any (diag (x) ~= 0))
        usage ('squareform: x is not a hollow matrix');
      elseif (~issymmetric(x))
        warning ('squareform:symmetric', ...
                 'squareform: x is not a symmetric matrix');
      endif

      sx = size (x, 1);
      y = zeros (1, (sx - 1) * sx / 2);
      idx = 1;
      for i = 2 : sx
        newidx = idx + sx - i;
        y(1, idx:newidx) = x(i:sx, i-1);
        idx = newidx + 1;
      endfor

    case "tomatrix"
      ## make sure that x is a column
      x = x(:);

      ## the dimensions of y are the solution to the quadratic formula for:
      ## length (x) = (sy - 1) * (sy / 2)
      sy = (1 + sqrt (1 + 8 * length (x))) / 2;
      if (floor (sy) ~= sy)
        usage ('squareform: incorrect vector size; see help');
      else
        y = zeros (sy);
        for i = 1 : sy-1
          step = sy - i;
          y((sy-step+1):sy, i) = x(1:step);
          x(1:step) = [];
        endfor
        y = y + y';
      endif
    otherwise
      error ("squareform: invalid METHOD '%s'", method);
  endswitch

endfunction

## make sure that it can go both directions automatically
%!assert(squareform(1:6), [0 1 2 3;1 0 4 5;2 4 0 6;3 5 6 0])
%!assert(squareform(squareform(1:6)),[1:6])
%!assert(squareform([0 1 2 3;1 0 4 5;2 4 0 6;3 5 6 0]), [1:6])

## make sure that the command arguments force the correct behavior
## squareform(1, "tovector") correctly throws an error: invalid distance matrix
%!assert(squareform(1), [0 1;1 0])
%!assert(squareform(1, "tomatrix"), [0 1;1 0])
%!assert(squareform(0, "tovector"), zeros(1,0))

## default to "tomatrix" for 1 element input
%!assert (squareform (1), [0 1; 1 0])
