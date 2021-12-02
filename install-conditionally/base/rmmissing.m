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
## Remove missing or incomplete data from an array or a matrix.
##
## Given an input array or matrix @var{A}, remove rows or columns with
## missing data from a matrix, or remove missing data from an array.  @var{A}
## can be a numeric or char matrix, a vector or an array of cell strings.
## @var{R} is the return matrix or array, after the removal of missing data.
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
## @seealso{ismissing, isnan}

function [R, TF] = rmmissing (A, varargin)

  if (nargin < 1)
     print_usage ();
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
  if (iscellstr (A) || isvector (A))
    TF = ismissing (A);

    R = A(find (TF == 0));
  elseif (ismatrix (A))
    ## matrix: ismissing returns a matrix, so it must be converted to a row or
    ## column vector according to the "dim" of choice
    TF = ismissing (A);

    if (optMinNumMissingI > 1)
      TF = sum (TF, optDimensionI);
      TF(find (TF < optMinNumMissingI)) = 0;
      TF = logical (TF);
    else
      TF = any (TF, optDimensionI);
    endif

    if (optDimensionI == 2)
      ## remove the rows
      R = A(find (TF == 0), :);
    else
      ## remove the columns 
      R = A(:, find (TF == 0));
    endif
  else
    error ("rmmissing: unsupported data");
  endif
endfunction


%!assert (rmmissing ([1,NaN,3]), [1,3])
%!assert (rmmissing ('abcd f'), 'abcdf')
%!assert (rmmissing ({'xxx','','xyz'}), {'xxx','xyz'})
%!assert (rmmissing ([1,2;NaN,2]), [1,2])
%!assert (rmmissing ([1,2;NaN,2], 2), [2,2]')
%!assert (rmmissing ([1,2;NaN,4;NaN,NaN],"MinNumMissing", 2), [1,2;NaN,4])

## Test input validation
%!error rmmissing ();
%!error rmmissing ({1, 2, 3});
%!error <must be either 1 or 2> rmmissing ([1 2; 3 4], 5);
%!error <unknown parameter name> rmmissing ([1 2; 3 4], "XXX", 1);
%!error <'MinNumMissing'> rmmissing ([1 2; 3 4], 2, "MinNumMissing", -2);
%!error <'MinNumMissing'> rmmissing ([1 2; 3 4], "MinNumMissing", 3.8);
%!error <'MinNumMissing'> rmmissing ([1 2; 3 4], "MinNumMissing", [1 2 3]);
%!error <'MinNumMissing'> rmmissing ([1 2; 3 4], "MinNumMissing", 'xxx');
