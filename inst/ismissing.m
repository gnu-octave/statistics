## Copyright (C) 2026 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {statistics} {@var{TF} =} ismissing (@var{A})
## @deftypefnx {statistics} {@var{TF} =} ismissing (@var{A}, @var{indicator})
##
## Find missing data in arrays.
##
## @code{@var{TF} = ismissing (@var{A})} returns a logical array, @var{TF}, with
## the same dimensions as @var{A}, where @code{true} values match the standard
## missing values in the input data according to their data type.
##
## Standard missing values and their corresponding data types are:
##
## @itemize
## @item @qcode{NaN} - for @qcode{double}, @qcode{single}, @qcode{duration}, and
## @qcode{calendarDuration} arrays.
## @item @qcode{NaT} - for @qcode{datetime} arrays.
## @item @qcode{<missing>} - for @qcode{string} arrays.
## @item @qcode{<undefined>} - for @qcode{categorical} arrays.
## @item @qcode{@{0x0 char@}} - for @qcode{cell} arrays of character vectors.
## @end itemize
##
## For any data types that do not support missing values, @code{ismissing}
## returns @code{@var{TF} = false (size (@var{A}))}.
##
## Note: the generic @code{ismissing} function from the statistics package only
## operates on core Octave datatypes and it explicitly identifies missing values
## in @qcode{double} and @qcode{single} arrays, as well as in @qcode{cell}
## arrays of character vectors.  All other data types are handled by the
## overloaded methods from their respective data class from the datatypes
## package.  Use @code{help class_name.ismissing} to find more information about
## the functional specialization of their respective class implementation.
##
## The optional input @var{indicator} can be a scalar or a vector, of the same
## type as the input data @var{A}, specifying alternative missing values in the
## input data.  When specifying @var{indicator} values, the standard missing
## values are ignored, unless explicitly stated in the @var{indicator}.
##
## Additional data type matches between @var{indicator} and @var{A} are:
##
## @itemize
## @item @qcode{double} indicators also match @qcode{single}, all integer types,
## and @qcode{logical} data in @var{A}.
##
## @item @qcode{string} and @qcode{char} indicators also match
## @qcode{categorical} data in @var{A}.
##
## @item @qcode{char} and @qcode{cellstr} indicators also match @qcode{string}
## data in @var{A}.
## @end itemize
##
## Note: the generic @code{ismissing} function from the statistics package only
## accepts @var{indicator} argument for numeric, @qcode{logical}, and
## @qcode{char} arrays, as well as for @qcode{cell} arrays of character vectors.
## For all other core Octave data types, @code{ismissing} produces an error.
## However, @var{indicator} is supported for data classes from the datatypes
## package through their respective class implementation of overloaded methods.
##
## @seealso{fillmissing, rmmissing, standardizeMissing}
## @end deftypefn

function TF = ismissing (A, indicator)

  if (nargin < 1) || (nargin > 2)
    print_usage ();
  endif

  ## Check INDICATOR
  if (nargin != 2)
    indicator = [];
  endif

  ## If A is a cell array of character vectors and INDICATOR is a character
  ## vector, convert it to a single-element cell array of character vectors
  if (iscellstr (A) && ischar (indicator) && ! iscellstr (indicator))
    indicator = {indicator};
  endif

  if ((! isempty (indicator)) &&
      ((isnumeric (A) && ! isnumeric (indicator)) ||
       (iscellstr (A) && ! iscellstr (indicator)) ||
       (ischar (A) && ! ischar (indicator)) ||
       (islogical (A) && ! (islogical (indicator) || isnumeric (indicator)))))
    error ("ismissing: 'indicator' and 'A' must have the same data type.");
  endif

  ## main logic
  if (isempty (indicator))

    if (isnumeric (A))
      ## Numeric matrix: just find the NaNs
      ## integer types have no missing value, but isnan will return false
      TF = isnan (A);

    elseif (iscellstr (A))
      ## Cell array of character vectors - find cells with empty char {''}
      TF = cellfun ('isempty', A);

    else
      ## No standard missing type defined, return false.
      ## Other datatypes with standard missing values are
      ## handled by their respective overloading methods.
      TF = false (size (A));
    endif

  else
    ## INDICATOR specified for missing data
    TF = false (size (A));
    if (isnumeric (A))
      for iter = 1:numel (indicator)
        if (isnan (indicator(iter)))
          TF(isnan (A)) = true;
        else
          TF(A == indicator(iter)) = true;
        endif
      endfor
    elseif (islogical (A) || ischar (A))
      for iter = 1:numel (indicator)
        TF(A == indicator(iter)) = true;
      endfor
    elseif (iscellstr (A))
      for iter = 1:numel (indicator)
        if (isempty (indicator{iter}))
          TF(cellfun ('isempty', A)) = true;
        else
          TF(strcmp (A, indicator(iter))) = true;
        endif
      endfor
    else
      error ("ismissing: indicators not supported for data type '%s'", ...
             class (A));
    endif
  endif
endfunction

%!assert (ismissing ([1, NaN, 3]), [false, true, false])
%!assert (ismissing ('abcd f'), [false, false, false, false, false, false])
%!assert (ismissing ({'xxx', '', 'xyz'}), [false, true, false])
%!assert (ismissing ({'x', '', 'y'}), [false, true, false])
%!assert (ismissing ({'x', '', 'y'; 'z', 'a', ''}), logical ([0, 1, 0; 0, 0, 1]))
%!assert (ismissing ([1, 2; NaN, 2]), [false, false; true, false])
%!assert (ismissing ([1, 2; NaN, 2], 2), [false, true; false, true])
%!assert (ismissing ([1, 2; NaN, 2], [1, 2]), [true, true; false, true])
%!assert (ismissing ([1, 2; NaN, 2], NaN), [false, false; true, false])

## test nD array data
%!assert (ismissing (cat (3, magic (2), magic (2))), logical (zeros (2, 2, 2)))
%!assert (ismissing (cat (3, magic (2), [1, 2; 3, NaN])), ...
%!        logical (cat (3, [0, 0; 0, 0], [0, 0; 0, 1])))
%!assert (ismissing ([1, 2; 3, 4], [5, 1; 2, 0]), logical ([1, 1; 0, 0]))
%!assert (ismissing (cat (3, 'f oo', 'ba r')), ...
%!        logical (cat (3, [0, 0, 0, 0], [0, 0, 0, 0])))
%!assert (ismissing (cat (3, {'foo'}, {''}, {'bar'})), logical (cat (3, 0, 1, 0)))

## test data type handling
%!assert (ismissing (double (NaN)), true)
%!assert (ismissing (single (NaN)), true)
%!assert (ismissing (' '), false)
%!assert (ismissing ({''}), true)
%!assert (ismissing ({' '}), false)
%!assert (ismissing (double (eye(3)), single (1)), logical (eye (3)))
%!assert (ismissing (double (eye(3)), int32 (1)), logical (eye (3)))
%!assert (ismissing (single (eye(3)), double (1)), logical (eye (3)))
%!assert (ismissing (single (eye(3)), int32 (1)), logical (eye (3)))

## test data types without missing values
%!assert (ismissing ({'123', '', 123}), [false, false, false])
%!assert (ismissing (logical ([1, 0, 1])), [false, false, false])
%!assert (ismissing (int32 ([1, 2, 3])), [false, false, false])
%!assert (ismissing (uint32 ([1, 2, 3])), [false, false, false])
%!assert (ismissing ({1, 2, 3}), [false, false, false])
%!assert (ismissing ([struct struct struct]), [false, false, false])
%!assert (ismissing (logical (eye(3)), true), logical (eye (3)))
%!assert (ismissing (logical (eye(3)), double (1)), logical (eye (3)))
%!assert (ismissing (logical (eye(3)), single (1)), logical (eye (3)))
%!assert (ismissing (logical (eye(3)), int32 (1)), logical (eye (3)))
%!assert (ismissing (int32 (eye(3)), int32 (1)), logical (eye (3)))
%!assert (ismissing (int32 (eye(3)), double (1)), logical (eye (3)))
%!assert (ismissing (int32 (eye(3)), single (1)), logical (eye (3)))

## test empty input handling
%!assert (ismissing ([]), logical([]))
%!assert (ismissing (''), logical([]))
%!assert (ismissing (ones (0,1)), logical(ones(0,1)))
%!assert (ismissing (ones (1,0)), logical(ones(1,0)))
%!assert (ismissing (ones (1,2,0)), logical(ones(1,2,0)))

## test indicators and standard missing values
%!assert (ismissing ([1, NaN, 0, 2]), [false, true, false, false])
%!assert (ismissing ([1, NaN, 0, 2], [0, 1]), [true, false, true, false])
%!assert (ismissing ([1, NaN, 0, 2], [0, NaN]), [false, true, true, false])
%!assert (ismissing ([true, false, true]), [false, false, false])
%!assert (ismissing ([true, false, true], 1), [true, false, true])
%!assert (ismissing ([true, false, true], 0), [false, true, false])
%!assert (ismissing ({'', 'a', 'f'}), [true, false, false])
%!assert (ismissing ({'', 'a', 'f'}, 'a'), [false, true, false])
%!assert (ismissing ({'', 'a', 'f'}, {'a', 'g'}), [false, true, false])
%!assert (ismissing ({'', 'a', 'f'}, {'a', 'f'}), [false, true, true])

## Test input validation
%!error ismissing ()
%!error ismissing (1, 2, 3)
%!error <ismissing: 'indicator' and 'A' must have the same data type.> ...
%!       ismissing ([1, 2; 3, 4], 'abc')
%!error <ismissing: 'indicator' and 'A' must have the same data type.> ...
%!       ismissing ({'', '', ''}, 1)
%!error <ismissing: 'indicator' and 'A' must have the same data type.> ...
%!       ismissing (1, struct)
%!error <ismissing: indicators not supported for data type 'struct'> ...
%!       ismissing (struct, 1)
%!error <ismissing: indicators not supported for data type 'cell'> ...
%!       ismissing ({1, 2, 3}, 2)
