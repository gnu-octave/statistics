########################################################################
##
## Copyright (C) 1995-2021 The Octave Project Developers
##
## See the file COPYRIGHT.md in the top-level directory of this
## distribution or <https://octave.org/copyright/>.
##
## This file is part of Octave.
##
## Octave is free software: you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## Octave is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with Octave; see the file COPYING.  If not, see
## <https://www.gnu.org/licenses/>.
##
########################################################################

## -*- texinfo -*-
## @deftypefn {} {@var{R} =} rmmissing (@var{A})
## @deftypefnx {} {@var{R} =} rmmissing (@var{A}, @var{dim})
## @deftypefnx {} {@var{R} =} rmmissing (@dots{}, @var{Name}, @var{Value})
## @deftypefnx {} {[@var{R} @var{TF}] =} rmmissing (@dots{})
##
## Remove missing or incomplete data from an array.
##
## Given an input vector or matrix (2-D array) @var{A}, remove missing data
## from a vector or missing rows or columns from a matrix.  @var{A}
## can be a numeric array, char array, or an array of cell strings.
## @var{R} returns the array after removal of missing data.
##
## The values which represent missing data depend on the data type of @var{A}:
##
## @itemize
## @item
## @code{NaN}: @code{single}, @code{double}.
##
## @item
## @code{' '} (white space): @code{char}.
##
## @item
## @code{@{''@}}: string cells.
## @end itemize
##
## Choose to remove rows (default) or columns by setting optional input
## @var{dim}:
##
## @itemize
## @item
## @code{1}: rows.
##
## @item
## @code{2}: columns.
## @end itemize
##
## Note: data types with no default 'missing' value will always result in
## @code{R == A} and a TF output of @code{false(size(@var{A}))}.
##
## Additional optional parameters are set by @var{Name}-@var{Value} pairs.
## These are:
##
## @itemize
## @item
## @code{MinNumMissing}: minimum number of missing values to remove an entry,
## row or column, defined as a positive integer number.  E.g.: if
## @code{MinNumMissing} is set to @code{2}, remove the row of a numeric matrix
## only if it includes 2 or more NaN.
## @end itemize
##
## Optional return value @var{TF} is a logical array where @code{true} values
## represent removed entries, rows or columns from the original data @var{A}.
##
## @end deftypefn
##
## @seealso{ismissing, standardizeMissing}

function [R, TF] = rmmissing (A, varargin)

  if ((nargin < 1) || (nargin > 4))
     print_usage ();
  endif

  if ndims(A) > 2
    error ("rmmissing: input dimension cannot exceed 2");
  endif

  optDimensionI = 2; # default dimension: rows
  optMinNumMissingI = 1;

  ## parse options
  if (nargin > 1)
    if (isnumeric (varargin{1}))
      ## option "dim"
      switch (varargin{1})
        case 1
          optDimensionI = 2;
        case 2
          optDimensionI = 1;
        otherwise
          error ("rmmissing: 'dim' must be either 1 or 2");
      endswitch

      pair_index = 2;
    else
      [r, c] = size (A);

      ## first non singleton dimension, but only two dimensions considered
      if (r == 1 && c != 1)
        optDimensionI = 1;
      endif

      pair_index = 1;
    endif

    ## parse name-value parameters
    while (pair_index <= (nargin - 1))
      switch (lower (varargin{pair_index}))
        ## minimum number of missing values to remove entries;
        ## it must be a positive integer number
        case "minnummissing"
          if (! isnumeric (varargin{pair_index + 1}) ||
              ! isscalar (varargin{pair_index + 1}) ||
              floor (varargin{pair_index + 1}) != varargin{pair_index + 1} ||
              varargin{pair_index + 1} < 1)
            error (["rmmissing: 'MinNumMissing' requires a positive integer"...
                    " number as value"]);
          endif

          optMinNumMissingI = varargin{pair_index + 1};
        otherwise
          error ("rmmissing: unknown parameter name '%s'", ...
                 varargin{pair_index});
      endswitch

      pair_index += 2;
    endwhile
  endif

  ## main logic
  TF = ismissing (A);
  if (isvector (A))
    R = A(TF == 0);

  elseif (iscellstr(A) || ismatrix (A))
    ## matrix: ismissing returns an array, so it must be converted to a row or
    ## column vector according to the "dim" of choice
    if (optMinNumMissingI > 1)
      TF = sum (TF, optDimensionI);
      TF(TF < optMinNumMissingI) = 0;
      TF = logical (TF);
    else
      TF = any (TF, optDimensionI);
    endif

    if (optDimensionI == 2)
      ## remove the rows
      R = A((TF == 0), :);
    else
      ## remove the columns
      R = A(:, (TF == 0));
    endif
  else
    error ("rmmissing: unsupported data");
  endif
endfunction

%!assert (rmmissing ([1,NaN,3]), [1,3])
%!assert (rmmissing ('abcd f'), 'abcdf')
%!assert (rmmissing ({'xxx','','xyz'}), {'xxx','xyz'})
%!assert (rmmissing ({'xxx','';'xyz','yyy'}), {'xyz','yyy'})
%!assert (rmmissing ({'xxx','';'xyz','yyy'}, 2), {'xxx';'xyz'})
%!assert (rmmissing ([1,2;NaN,2]), [1,2])
%!assert (rmmissing ([1,2;NaN,2], 2), [2,2]')
%!assert (rmmissing ([1,2;NaN,4;NaN,NaN],"MinNumMissing", 2), [1,2;NaN,4])

## Test second output
%!test
%! x = [1:6];
%! x([2,4]) = NaN;
%! [~, idx] = rmmissing (x);
%! assert (idx, logical ([0, 1, 0, 1, 0, 0]));
%! assert (class(idx), 'logical');
%! x = reshape (x, [2, 3]);
%! [~, idx] = rmmissing (x);
%! assert (idx, logical ([0; 1]));
%! assert (class(idx), 'logical');
%! [~, idx] = rmmissing (x, 2);
%! assert (idx, logical ([1, 1, 0]));
%! assert (class(idx), 'logical');
%! [~, idx] = rmmissing (x, 1, "MinNumMissing", 2);
%! assert (idx, logical ([0; 1]));
%! assert (class(idx), 'logical');
%! [~, idx] = rmmissing (x, 2, "MinNumMissing", 2);
%! assert (idx, logical ([0, 0, 0]));
%! assert (class(idx), 'logical');

## Test data type handling
%!assert (rmmissing (single ([1 2 NaN; 3 4 5])), single ([3 4 5]))
%!assert (rmmissing (logical (ones (3))), logical (ones (3)))
%!assert (rmmissing (int32 (ones (3))), int32 (ones (3)))
%!assert (rmmissing (uint32 (ones (3))), uint32 (ones (3)))
%!assert (rmmissing ({1, 2, 3}), {1, 2, 3})
%!assert (rmmissing ([struct, struct, struct]), [struct, struct, struct])

## Test empty input handling
%!assert (rmmissing ([]), [])
%!assert (rmmissing (ones (1,0)), ones (1,0))
%!assert (rmmissing (ones (1,0), 1), ones (1,0))
%!assert (rmmissing (ones (1,0), 2), ones (1,0))
%!assert (rmmissing (ones (0,1)), ones (0,1))
%!assert (rmmissing (ones (0,1), 1), ones (0,1))
%!assert (rmmissing (ones (0,1), 2), ones (0,1))
%!error <input dimension> rmmissing (ones (0,1,2))

## Test input validation
%!error rmmissing ()
%!error <input dimension> rmmissing (ones(2,2,2))
%!error <must be either 1 or 2> rmmissing ([1 2; 3 4], 5)
%!error <unknown parameter name> rmmissing ([1 2; 3 4], "XXX", 1)
%!error <'MinNumMissing'> rmmissing ([1 2; 3 4], 2, "MinNumMissing", -2)
%!error <'MinNumMissing'> rmmissing ([1 2; 3 4], "MinNumMissing", 3.8)
%!error <'MinNumMissing'> rmmissing ([1 2; 3 4], "MinNumMissing", [1 2 3])
%!error <'MinNumMissing'> rmmissing ([1 2; 3 4], "MinNumMissing", 'xxx')

