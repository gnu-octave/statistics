## Copyright (C) 1995-2023 The Octave Project Developers
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
## @deftypefn  {statistics} {@var{B} =} standardizeMissing (@var{A}, @var{indicator})
##
## Replace data values specified by @var{indicator} in @var{A} by the
## standard 'missing' data value for that data type.
##
## @var{A} can be a numeric scalar or array, a character vector or array, or
## a cell array of character vectors (a.k.a. string cells).
##
## @var{indicator} can be a scalar or an array containing values to be
## replaced by the 'missing' value for the class of @var{A}, and should have
## a data type matching @var{A}.
##
## 'missing' values are defined as :
##
## @itemize
## @item
## @qcode{NaN}: @code{single}, @code{double}
##
## @item
## @qcode{" "} (white space): @code{char}
##
## @item
## @qcode{@{""@}} (empty string in cell): string cells.
## @end itemize
##
## Compatibility Notes:
## @itemize
## @item
## Octave's implementation of @code{standardizeMissing}
## does not restrict @var{indicator} of type @qcode{char} to row vectors.
##
## @item
## All numerical and logical inputs for @var{A} and @var{indicator} may
## be specified in any combination. The output will be the same class as
## @var{A}, with the @var{indicator} converted to that data type for
## comparison.  Only @code{single} and @code{double} have defined 'missing'
## values, so @var{A} of other data types will always output
## @var{B} = @var{A}.
## @end itemize
##
## @end deftypefn
##
## @seealso{fillmissing, ismissing, rmmissing}

function A = standardizeMissing (A, indicator)

  if (nargin != 2)
     print_usage ();
  endif

  input_class = class(A);
  do_nothing_flag = false;

  ## Set missing_val.
  ## Compatibility requirements:
  ## Numeric/logical: Only double and single have 'missing' values defined,
  ## so other numeric and logical pass through unchanged.
  ## Other: only char & cellstr have 'missing' values defined, others produce
  ## error.
  ##
  ## TODO: as implemented in Octave, add table, string, timetable,
  ## categorical, datetime, and duration to class checks and BISTs
  ##
  ## TODO: when 'missing' class is implemented, these switch blocks can all
  ## be removed and the final assignment updated to call missing instead of
  ## missing_val

  if (isnumeric(A) || islogical(A))
    switch (input_class)
      case "double"
        missing_val = NaN ("double");
      case "single"
        missing_val = NaN ("single");
      otherwise
        do_nothing_flag = true;
    endswitch
  else
    switch (input_class)
      case "char"
        missing_val = " ";
      case "cell"
        if iscellstr(A)
          missing_val = {""};
        else
          error ("standardizeMissing: only cells of strings are supported.");
        endif
    otherwise
      error ("standardizeMissing: unsupported data type %s.", input_class);
    endswitch
  endif

  if (! do_nothing_flag)

    ## if A is an array of cell strings and indicator just a string,
    ## convert indicator to a cell string with one element
    if (iscellstr (A) && ischar (indicator) && ! iscellstr (indicator))
      indicator = {indicator};
    endif

    if ((isnumeric (A) && ! (isnumeric (indicator) || islogical (indicator))) ||
        (ischar (A) && ! ischar (indicator)) ||
        (iscellstr (A) && ! (iscellstr (indicator))))
      error (strcat ("standardizeMissing: 'indicator' and 'A' must", ...
                     " have the same data type."));
    endif

    A(ismember (A, indicator)) = missing_val;

  endif
endfunction

## numeric tests
%!assert (standardizeMissing (1, 1), NaN)
%!assert (standardizeMissing (1, 0), 1)
%!assert (standardizeMissing (eye(2), 1), [NaN 0;0 NaN])
%!assert (standardizeMissing ([1:3;4:6], [2 3; 4 5]), [1, NaN, NaN; NaN, NaN, 6])
%!assert (standardizeMissing (cat (3,1,2,3,4), 3), cat (3,1,2,NaN,4))

## char and cellstr tests
%!assert (standardizeMissing ('foo', 'a'), 'foo')
%!assert (standardizeMissing ('foo', 'f'), ' oo')
%!assert (standardizeMissing ('foo', 'o'), 'f  ')
%!assert (standardizeMissing ('foo', 'oo'), 'f  ')

%!assert (standardizeMissing ({'foo'}, 'f'), {'foo'})
%!assert (standardizeMissing ({'foo'}, {'f'}), {'foo'})
%!assert (standardizeMissing ({'foo'}, 'test'), {'foo'})
%!assert (standardizeMissing ({'foo'}, {'test'}), {'foo'})
%!assert (standardizeMissing ({'foo'}, 'foo'), {''})
%!assert (standardizeMissing ({'foo'}, {'foo'}), {''})

## char and cellstr array tests
%!assert (standardizeMissing (['foo';'bar'], 'oar'), ['f  ';'b  '])
%!assert (standardizeMissing (['foo';'bar'], ['o';'a';'r']), ['f  ';'b  '])
%!assert (standardizeMissing (['foo';'bar'], ['o ';'ar']), ['f  ';'b  '])

%!assert (standardizeMissing ({'foo','bar'}, 'foo'), {'','bar'})
%!assert (standardizeMissing ({'foo','bar'}, 'f'), {'foo','bar'})
%!assert (standardizeMissing ({'foo','bar'}, {'foo', 'a'}), {'','bar'})
%!assert (standardizeMissing ({'foo'}, {'f', 'oo'}), {'foo'})
%!assert (standardizeMissing ({'foo','bar'}, {'foo'}), {'','bar'})
%!assert (standardizeMissing ({'foo','bar'}, {'foo', 'a'}), {'','bar'})

## numeric type preservation tests
%!assert (standardizeMissing (double (1), single (1)), double (NaN))
%!assert (standardizeMissing (single (1), single (1)), single (NaN))
%!assert (standardizeMissing (single (1), double (1)), single (NaN))
%!assert (standardizeMissing (single (1), true), single (NaN))
%!assert (standardizeMissing (double (1), int32(1)), double (NaN))

## Passtrough tests
%!assert (standardizeMissing (true, true), true)
%!assert (standardizeMissing (true, 1), true)
%!assert (standardizeMissing (int32 (1), int32 (1)), int32 (1))
%!assert (standardizeMissing (int32 (1), 1), int32 (1))
%!assert (standardizeMissing (uint32 (1), uint32 (1)), uint32 (1))
%!assert (standardizeMissing (uint32 (1), 1), uint32 (1))

## Test input validation
%!error standardizeMissing ();
%!error standardizeMissing (1);
%!error standardizeMissing (1,2,3);
%!error <only cells of strings> standardizeMissing ({'abc', 1}, 1);
%!error <unsupported data type> standardizeMissing (struct ('a','b'), 1);
%!error <'indicator' and 'A' must have > standardizeMissing ([1 2 3], {1});
%!error <'indicator' and 'A' must have > standardizeMissing ([1 2 3], 'a');
%!error <'indicator' and 'A' must have > standardizeMissing ([1 2 3], struct ('a', 1));
%!error <'indicator' and 'A' must have > standardizeMissing ('foo', 1);
%!error <'indicator' and 'A' must have > standardizeMissing ('foo', {1});
%!error <'indicator' and 'A' must have > standardizeMissing ('foo', {'f'});
%!error <'indicator' and 'A' must have > standardizeMissing ('foo', struct ('a', 1));
%!error <'indicator' and 'A' must have > standardizeMissing ({'foo'}, 1);
%!error <'indicator' and 'A' must have > standardizeMissing ({'foo'}, 1);

