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
## @deftypefn {statistics} {@var{B} =} standardizeMissing (@var{A}, @var{indicator})
##
## Replace selected values by standard missing values.
##
## @code{@var{Β} = standardizeMissing (@var{A}, @var{indicator})} returns a
## standardized array @var{B} of the same size and data type as the input array
## @var{A} and with all elements specified by @var{indicator} replaced by the
## standard missing value corresponding the data type of @var{A}.
## @var{indicator} can be either a scalar or a vector.
##
## Standard missing values and their corresponding data types are:
##
## @itemize
## @item @qcode{NaN} - for @qcode{double}, @qcode{single}, @qcode{duration}, and
## @qcode{calendarDuration} arrays.
## @item @qcode{NaT} - for @qcode{datetime} arrays.
## @item @qcode{<missing>} - for @qcode{string} arrays.
## @item @qcode{<undefined>} - for @qcode{categorical} arrays.
## @item @qcode{@{''@}} - for @qcode{cell} arrays of character vectors.
## @end itemize
##
## For any other data type input that does not support missing values,
## @code{standardizeMissing} returns @code{@var{B} = @var{A}} and any
## @var{indicator} value is ignored.
##
## The nonstandard missing value @var{indicator} must be of the same type as the
## data input @var{A} or have a compatible data types according to the following
## rules:
##
## @itemize
## @item all numeric indicators match both @qcode{double} and @qcode{single}
## data types in @var{A}.
## @item indicators specified as @qcode{string} arrays, @qcode{char} vectors,
## and @code{cell} arrays of character vectors match categorical data type in
## @var{A}.
## @item a @qcode{char} vector matches a @qcode{cell} array of character vectors
## in @var{A}.
## @end itemize
##
## Note: the generic @code{standardizeMissing} function from the statistics does
## not operate on table inputs, which is handled by the overloaded method of the
## table class.  Use @code{help table.standardizeMissing} to find more
## information about the functional specialization on tables.
##
## @seealso{fillmissing, ismissing, rmmissing}
## @end deftypefn

function A = standardizeMissing (A, indicator)

  if (nargin != 2)
    print_usage ();
  endif

  if (isnumeric (A))
    if (! isnumeric (indicator))
      error ("standardizeMissing: incompatible INDICATOR and input data A.");
    elseif (! isvector (indicator))
      error ("standardizeMissing: INDICATOR must be a scalar or a vector.");
    endif
    switch (class (A))
      case "double"
        A(ismember (A, indicator)) = NaN ("double");
      case "single"
        A(ismember (A, indicator)) = NaN ("single");
    endswitch

  elseif (iscellstr (A))
    if (ischar (indicator))
      if (! isrow (indicator))
        error ("standardizeMissing: character INDICATOR must be a row vector.");
      endif
      indicator = {indicator};
    elseif (! iscellstr (indicator))
      error ("standardizeMissing: incompatible INDICATOR and input data A.");
    endif
    A(ismember (A, indicator)) = {''};

  elseif (iscategorical (A))
    if (ischar (indicator))
      if (! isrow (indicator))
        error ("standardizeMissing: character INDICATOR must be a row vector.");
      endif
      indicator = {indicator};
    elseif (! (iscellstr (indicator) || isstring (indicator)
                                     || iscategorical (indicator)))
      error ("standardizeMissing: incompatible INDICATOR and input data A.");
    elseif (! isvector (indicator))
      error ("standardizeMissing: INDICATOR must be a scalar or a vector.");
    endif
    A(ismember (A, indicator)) = categorical (NaN);

  elseif (isdatetime (A))
    if (! isdatetime (indicator))
      error ("standardizeMissing: incompatible INDICATOR and input data A.");
    elseif (! isvector (indicator))
      error ("standardizeMissing: INDICATOR must be a scalar or a vector.");
    endif
    A(ismember (A, indicator)) = NaT;

  elseif (isduration (A))
    if (! isduration (indicator))
      error ("standardizeMissing: incompatible INDICATOR and input data A.");
    elseif (! isvector (indicator))
      error ("standardizeMissing: INDICATOR must be a scalar or a vector.");
    endif
    A(ismember (A, indicator)) = days (NaN);

  elseif (iscalendarduration (A))
    if (! iscalendarduration (indicator))
      error ("standardizeMissing: incompatible INDICATOR and input data A.");
    elseif (! isvector (indicator))
      error ("standardizeMissing: INDICATOR must be a scalar or a vector.");
    endif
    A(ismember (A, indicator)) = days (NaN);

  elseif (isstring (A))
    if (! isstring (indicator))
      error ("standardizeMissing: incompatible INDICATOR and input data A.");
    elseif (! isvector (indicator))
      error ("standardizeMissing: INDICATOR must be a scalar or a vector.");
    endif
    A(ismember (A, indicator)) = missing;
  endif

endfunction

## numeric tests
%!assert (standardizeMissing (1, 1), NaN)
%!assert (standardizeMissing (1, 0), 1)
%!assert (standardizeMissing (eye(2), 1), [NaN 0;0 NaN])
%!assert (standardizeMissing ([1:3;4:6], [2 3 4 5]), [1, NaN, NaN; NaN, NaN, 6])
%!assert (standardizeMissing (cat (3,1,2,3,4), 3), cat (3,1,2,NaN,4))

## char and cellstr tests
%!assert (standardizeMissing ('foo', 'a'), 'foo')
%!assert (standardizeMissing ('foo', 'f'), 'foo')
%!assert (standardizeMissing ('foo', 'o'), 'foo')
%!assert (standardizeMissing ('foo', 'oo'), 'foo')

%!assert (standardizeMissing ({'foo'}, 'f'), {'foo'})
%!assert (standardizeMissing ({'foo'}, {'f'}), {'foo'})
%!assert (standardizeMissing ({'foo'}, 'test'), {'foo'})
%!assert (standardizeMissing ({'foo'}, {'test'}), {'foo'})
%!assert (standardizeMissing ({'foo'}, 'foo'), {''})
%!assert (standardizeMissing ({'foo'}, {'foo'}), {''})

## char and cellstr array tests
%!assert (standardizeMissing (['foo';'bar'], 'oar'), ['foo';'bar'])
%!assert (standardizeMissing (['foo';'bar'], ['o';'a';'r']), ['foo';'bar'])
%!assert (standardizeMissing (['foo';'bar'], ['o ';'ar']), ['foo';'bar'])

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
%!assert (standardizeMissing (single (1), uint8 (1)), single (NaN))
%!assert (standardizeMissing (double (1), int32 (1)), double (NaN))

## Passtrough tests
%!assert (standardizeMissing (true, true), true)
%!assert (standardizeMissing (true, 1), true)
%!assert (standardizeMissing (int32 (1), int32 (1)), int32 (1))
%!assert (standardizeMissing (int32 (1), 1), int32 (1))
%!assert (standardizeMissing (uint32 (1), uint32 (1)), uint32 (1))
%!assert (standardizeMissing (uint32 (1), 1), uint32 (1))
%!assert (standardizeMissing ({'abc', 1}, 1), {'abc', 1})
%!assert (standardizeMissing (struct ('a','b'), 1), struct ('a','b'))

## categorical array tests
%!assert (double (standardizeMissing (categorical (1), categorical (1))), NaN)
%!assert (double (standardizeMissing (categorical (1), '1')), NaN)
%!assert (class (standardizeMissing (categorical (1), categorical (1))), 'categorical')
%!assert (double (standardizeMissing (categorical (1), categorical (2))), 1)
%!assert (double (standardizeMissing (categorical (1), '2')), 1)
%!assert (class (standardizeMissing (categorical (1), categorical (2))), 'categorical')

## datetime array tests
%!assert (isnat (standardizeMissing (datetime ('today'), datetime ('today'))), true)
%!assert (isnat (standardizeMissing (datetime ('today'), datetime ('yesterday'))), false)

## duration array tests
%!assert (days (standardizeMissing (days (1), days (1))), NaN)
%!assert (days (standardizeMissing (days (1), days (2))), 1)

## string array tests
%!assert (cellstr (standardizeMissing (string (1), string (1))), {''})
%!assert (cellstr (standardizeMissing (string (1), string (2))), {'1'})

## Test input validation
%!error <Invalid call> standardizeMissing ();
%!error <Invalid call> standardizeMissing (1);
%!error <standardizeMissing: function called with too many inputs> standardizeMissing (1, 2, 3);
%!error <standardizeMissing: incompatible INDICATOR and input data A.> ...
%!       standardizeMissing ([1, 2, 3], {1});
%!error <standardizeMissing: incompatible INDICATOR and input data A.> ...
%!       standardizeMissing ([1, 2, 3], 'a');
%!error <standardizeMissing: incompatible INDICATOR and input data A.> ...
%!       standardizeMissing ([1, 2, 3], struct ('a', 1));
%!error <standardizeMissing: incompatible INDICATOR and input data A.> ...
%!       standardizeMissing (categorical (1), 1);
%!error <standardizeMissing: incompatible INDICATOR and input data A.> ...
%!       standardizeMissing ({'foo'}, string ('foo'));
%!error <standardizeMissing: character INDICATOR must be a row vector.> ...
%!       standardizeMissing ({'foo'}, ['a';'b']);

