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
## @deftypefn {} {@var{TF} =} ismissing (@var{A})
## @deftypefnx {} {@var{TF} =} ismissing (@var{A}, @var{indicator})
##
## Find missing data in a numeric or string array.
##
## Given an input numeric data array, char array, or array of cell strings
## @var{A}, @code{ismissing} returns a logical array @var{TF} with
## the same dimensions as @var{A}, where @code{true} values match missing
## values in the input data.
##
## The optional input @var{indicator} is an array of values that represent
## missing values in the input data.  The values which represent missing data
## by default depend on the data type of @var{A}:
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
## Note: logical and numeric data types may be used in any combination
## for @var{A} and @var{indicator}. @var{A} and the indicator values will be
## compared as type double, and the output will have the same class as @var{A}.
## Data types other than those specified above have no defined 'missing' value.
## As such, the TF output for those inputs will always be
## @code{false(size(@var{A}))}. The exception to this is that @var{indicator}
## can be specified for logical and numeric inputs to designate values that
## will register as 'missing'.
##
## @end deftypefn
##
## @seealso{rmmissing, standardizeMissing}

function TF = ismissing (A, indicator)

  if (nargin < 1) || (nargin > 2)
    print_usage ();
  endif

  ## check "indicator"
  if (nargin != 2)
     indicator = [];
  endif

  ## if A is an array of cell strings and indicator just a string,
  ## convert indicator to a cell string with one element
  if (iscellstr (A) && ischar (indicator) && ! iscellstr (indicator))
    indicator = {indicator};
  endif

  if ((! isempty (indicator)) &&
      ((isnumeric (A) && ! (isnumeric (indicator) || islogical (indicator))) ||
       (ischar (A) && ! ischar (indicator)) ||
       (iscellstr (A) && ! (iscellstr (indicator)))))
    error ("ismissing: 'indicator' and 'A' must have the same data type");
  endif

  ## main logic
  if (isempty (indicator))

    if (isnumeric (A))
      ## numeric matrix: just find the NaNs
      ## integer types have no missing value, but isnan will return false
      TF = isnan (A);

    elseif (iscellstr (A))
      ## cell strings - find empty cells
      TF = cellfun ('isempty', A);

    elseif (ischar (A))
      ## char matrix: find the white spaces
      TF = isspace (A);

    else
      ##no missing type defined, return false
      TF = false (size (A));
    endif

  else
    ## indicator specified for missing data
    TF = false (size (A));
    if (isnumeric(A) || ischar (A) || islogical (A))
      for iter = 1 : numel (indicator)
        if (isnan (indicator(iter)))
          TF(isnan(A)) = true;
        else
          TF(A == indicator(iter)) = true;
        endif
      endfor
    elseif (iscellstr (A))
      for iter = 1 : numel (indicator)
        TF(strcmp (A, indicator(iter))) = true;
      endfor
    else
      error ("ismissing: indicators not supported for data type '%s'", ...
               class(A));
    endif
  endif
endfunction

%!assert (ismissing ([1,NaN,3]), [false,true,false])
%!assert (ismissing ('abcd f'), [false,false,false,false,true,false])
%!assert (ismissing ({'xxx','','xyz'}), [false,true,false])
%!assert (ismissing ({'x','','y'}), [false,true,false])
%!assert (ismissing ({'x','','y';'z','a',''}), logical([0,1,0;0,0,1]))
%!assert (ismissing ([1,2;NaN,2]), [false,false;true,false])
%!assert (ismissing ([1,2;NaN,2], 2), [false,true;false,true])
%!assert (ismissing ([1,2;NaN,2], [1 2]), [true,true;false,true])
%!assert (ismissing ([1,2;NaN,2], NaN), [false,false;true,false])

## test nD array data
%!assert (ismissing (cat(3,magic(2),magic(2))), logical (zeros (2,2,2)))
%!assert (ismissing (cat(3,magic(2),[1 2;3 NaN])), logical (cat(3,[0,0;0,0],[0,0;0,1])))
%!assert (ismissing ([1 2; 3 4], [5 1; 2 0]), logical([1 1; 0 0]))
%!assert (ismissing (cat(3,'f oo','ba r')), logical(cat(3,[0 1 0 0],[0 0 1 0])))
%!assert (ismissing (cat(3,{'foo'},{''},{'bar'})), logical(cat(3,0,1,0)))

## test data type handling
%!assert (ismissing (double (NaN)), true)
%!assert (ismissing (single (NaN)), true)
%!assert (ismissing (' '), true)
%!assert (ismissing ({''}), true)
%!assert (ismissing ({' '}), false)
%!assert (ismissing (double (eye(3)), single (1)), logical(eye(3)))
%!assert (ismissing (double (eye(3)), true), logical(eye(3)))
%!assert (ismissing (double (eye(3)), int32 (1)), logical(eye(3)))
%!assert (ismissing (single (eye(3)), true), logical(eye(3)))
%!assert (ismissing (single (eye(3)), double (1)), logical(eye(3)))
%!assert (ismissing (single(eye(3)), int32 (1)), logical(eye(3)))

## test data types without missing values
%!assert (ismissing ({'123', '', 123}), [false false false])
%!assert (ismissing (logical ([1 0 1])), [false false false])
%!assert (ismissing (int32 ([1 2 3])), [false false false])
%!assert (ismissing (uint32 ([1 2 3])), [false false false])
%!assert (ismissing ({1, 2, 3}), [false false false])
%!assert (ismissing ([struct struct struct]), [false false false])
%!assert (ismissing (logical (eye(3)), true), logical(eye(3)))
%!assert (ismissing (logical (eye(3)), double (1)), logical(eye(3)))
%!assert (ismissing (logical (eye(3)), single (1)), logical(eye(3)))
%!assert (ismissing (logical (eye(3)), int32 (1)), logical(eye(3)))
%!assert (ismissing (int32 (eye(3)), int32 (1)), logical(eye(3)))
%!assert (ismissing (int32 (eye(3)), true), logical(eye(3)))
%!assert (ismissing (int32 (eye(3)), double (1)), logical(eye(3)))
%!assert (ismissing (int32 (eye(3)), single (1)), logical(eye(3)))

## test empty input handling
%!assert (ismissing ([]), logical([]))
%!assert (ismissing (''), logical([]))
%!assert (ismissing (ones (0,1)), logical(ones(0,1)))
%!assert (ismissing (ones (1,0)), logical(ones(1,0)))
%!assert (ismissing (ones (1,2,0)), logical(ones(1,2,0)))

## Test input validation
%!error ismissing ()
%!error <'indicator' and 'A' must have the same> ismissing ([1 2; 3 4], "abc")
%!error <'indicator' and 'A' must have the same> ismissing ({"", "", ""}, 1)
%!error <'indicator' and 'A' must have the same> ismissing (1, struct)
%!error <indicators not supported for data type> ismissing (struct, 1)
