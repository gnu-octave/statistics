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
## Find missing data in a matrix or a string array.
##
## Given an input vector, matrix or array of cell strings @var{A},
## @code{ismissing} returns a logical vector or matrix @var{TF} with the same
## dimensions as @var{A}, where @code{true} values match missing values in the
## input data.
##
## The optional input @var{indicator} is an array of values, which represent
## missing values in the input data.  The values which represent missing data by
## default depend on the data type of @var{A}:
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
## @end deftypefn
##
## @seealso{all, any, isempty, isnan, rmmissing}

function TF = ismissing (A, indicator)

  ## check "indicator"
  if (nargin != 2)
     indicator = [];
  else
    if (! isvector (indicator))
      error ("ismissing: invalid format for 'indicator'");
    endif

    if ((isnumeric (A) && ! isnumeric (indicator)) ||
        (ischar (A) && ! ischar (indicator)) ||
        (iscellstr (A) && ! (iscellstr (indicator) || ischar (indicator))))
      error ("ismissing: 'indicator' and 'A' must have the same data type");
    endif

    ## if A is an array of cell strings and indicator just a string,
    ## convert indicator to a cell string with one element
    if (iscellstr (A) && ischar (indicator) && ! iscellstr (indicator))
      tmpstr = indicator;
      indicator = {};

      indicator = {tmpstr};
    endif
  endif

  ## main logic
  if (iscellstr (A) && isempty (indicator))
    TF = false (1, length (A));
    ## remove all empty strings
    for iter = 1 : length (A)
      if (isempty (A{iter}))
        TF(iter) = true;
      endif
    endfor
  elseif (ismatrix (A) && isempty (indicator))
    if (isnumeric (A))
      ## numeric matrix: just remove the NaNs
      TF = isnan (A);
    elseif (ischar (A))
      ## char matrix: remove the white spaces
      TF = isspace (A);
    else
      error ("ismissing: unsupported data type");
    endif

    TF = logical (TF);
  elseif (! isempty (indicator))
    ## special cases with custom values for missing data
    [r, c] = size (A);
    TF = false (r, c);

    if (iscellstr (A))
      for iter = 1 : length (indicator)
        TF(find (strcmp (A, indicator(iter)))) = true;
      endfor
    elseif (ismatrix (A))
      for iter = 1 : length (indicator)
        TF(find (A == indicator(iter))) = true;
      endfor
    else
      error ("ismissing: unsupported data format");
    endif
  else
    error ("ismissing: unsupported data format");
  endif
endfunction


%!assert (ismissing ([1,NaN,3]), [false,true,false])
%!assert (ismissing ('abcd f'), [false,false,false,false,true,false])
%!assert (ismissing ({'xxx','','xyz'}), [false,true,false])
%!assert (ismissing ([1,2;NaN,2]), [false,false;true,false])

## Test input validation
%!error ismissing ();
%!error ismissing ({1, 2, 3});
%!error ismissing ([1 2; 3 4], [5 1; 2 0]);
%!error ismissing ([1 2; 3 4], "abc");
%!error ismissing ({"", "", ""}, 1);
%!error ismissing ({1, 2, 3});
