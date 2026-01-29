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
## @deftypefn  {statistics} {@var{R} =} rmmissing (@var{A})
## @deftypefnx {statistics} {@var{R} =} rmmissing (@var{A}, @var{dim})
## @deftypefnx {statistics} {@var{R} =} rmmissing (@dots{}, @var{Name}, @var{Value})
## @deftypefnx {statistics} {[@var{R}, @var{TF}] =} rmmissing (@dots{})
##
## Remove missing data from arrays.
##
## Given an input vector or matrix (2-D array) @var{A}, @code{@var{R} =
## rmmissing (@var{A})} returns an output vector or matrix @var{R} of the same
## type as input @var{A} and any missing elements removed.  If @var{A} is a
## vector, missing elements are removed individually, if @var{A} is a matrix,
## then rows containing missing elements are removed.
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
## For any data types that do not support missing values, @code{rmmissing}
## returns @code{@var{R} == @var{A}} and if a second output argument is
## requested it also returns @code{@var{TF} = false (size (@var{A}))}.
##
## Given an input matrix (2-D array) @var{A}, @code{@var{R} = rmmissing
## (@var{A}, @var{dim})} further specifies whether rows or columns containing
## missing data are removed from the output @var{R} based on the value of
## @var{dim}, which must be either 1 or 0.
##
## @itemize
## @item
## @qcode{1}: remove rows.
##
## @item
## @qcode{2}: remove columns.
## @end itemize
##
## @code{@var{R} = rmmissing (@dots{}, @var{Name}, @var{Value})} also accepts
## the following paired arguments.
##
## @multitable @columnfractions 0.2 0.05 0.75
## @headitem Name @tab @tab Value
## @item @qcode{'MinNumMissing'} @tab @tab A positive integer scalar value
## specifying the required minimum number of missing values for removing any
## particular row or column from a matrix input.  Note that this argument is
## ignored if input @var{A} is a vector.
##
## @item @qcode{'MissingLocations'} @tab @tab A logical array of the same size
## as input @var{A} indexing the locations of missing values in input array
## @var{A}.  Note that specifying @qcode{'MissingLocations'} overrides any
## standard missing values in @var{A}.
## @end multitable
##
## Optional return value @var{TF} is a logical array where @code{true} values
## represent removed entries, rows or columns from the original data @var{A}.
##
## @seealso{fillmissing, ismissing, standardizeMissing}
## @end deftypefn

function [R, TF] = rmmissing (A, varargin)

  ## Validate data in A
  if (nargin < 1)
    print_usage ();
  endif
  if (ndims(A) > 2)
    error ("rmmissing: A must be a matrix; no more than 2 dimensions allowed.");
  endif
  if (isempty (A))
    R = A;
    TF = false (size (A));
    return;
  endif

  ## Parse optional Name-Value paired arguments
  optNames = {'MinNumMissing', 'MissingLocations'};
  dfValues = {1, []};
  [MinNumMissing, MissingLocations, args] = pairedArgs (optNames, dfValues, ...
                                                        varargin(:));

  ## Validate optional Name-Value paired arguments
  if (! (isscalar (MinNumMissing) && isnumeric (MinNumMissing) &&
         MinNumMissing > 0 && fix (MinNumMissing) == MinNumMissing))
    error ("rmmissing: 'MinNumMissing' must be a positive integer value.");
  endif
  if (! isempty (MissingLocations))
    if (! (islogical (MissingLocations) &&
           isequal (size (MissingLocations), size (A))))
      error (strcat ("rmmissing: 'MissingLocations' must be a", ...
                     " logical matrix of the same size as input A."));
    endif
  endif

  ## Check for DIM
  if (isempty (args))
    dim = 2;
  elseif (isscalar (args))
    dim = args{1};
    if (! (isscalar (dim) && (dim == 1 || dim == 2)))
      error ("rmmissing: specified DIM must be either 1 or 2.");
    endif
    ## Swap DIM to operate orthogonal to specified DIM
    if (dim == 1)
      dim = 2;
    else
      dim = 1;
    endif
  else
    error ("rmmissing: too many input arguments.");
  endif

  ## Get missing values
  if (isempty (MissingLocations))
    TF = ismissing (A);
  else
    TF = MissingLocations;
  endif

  ## Remove missing values
  if (isvector (A))
    R = A(TF == 0); # MinNumMissing does not matter here
  else
    ## matrix: ismissing returns an array, so it must be converted
    ## to a row or column vector according to the "dim" of choice
    if (MinNumMissing > 1)
      TF = sum (TF, dim);
      TF(TF < MinNumMissing) = 0; # true only if at least MinNumMissing
      TF = logical (TF);
    else
      TF = any (TF, dim);
    endif

    if (dim == 2)
      ## remove the rows
      R = A((TF == 0), :);
    else
      ## remove the columns
      R = A(:, (TF == 0));
    endif
  endif

endfunction

%!assert (rmmissing ([1, NaN, 3]), [1, 3])
%!assert (rmmissing ('abcd f'), 'abcd f')
%!assert (rmmissing ({'xxx', '', 'xyz'}), {'xxx', 'xyz'})
%!assert (rmmissing ({'xxx', ''; 'xyz', 'yyy'}), {'xyz', 'yyy'})
%!assert (rmmissing ({'xxx', ''; 'xyz', 'yyy'}, 2), {'xxx'; 'xyz'})
%!assert (rmmissing ([1, 2; NaN, 2]), [1, 2])
%!assert (rmmissing ([1, 2; NaN, 2], 2), [2, 2]')
%!assert (rmmissing ([1, 2; NaN, 4; NaN, NaN],"MinNumMissing", 2), [1, 2; NaN, 4])

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
%!assert (rmmissing (single ([1, 2, NaN; 3, 4, 5])), single ([3, 4, 5]))
%!assert (rmmissing (logical (ones (3))), logical (ones (3)))
%!assert (rmmissing (int32 (ones (3))), int32 (ones (3)))
%!assert (rmmissing (uint32 (ones (3))), uint32 (ones (3)))
%!assert (rmmissing ({1, 2, 3}), {1, 2, 3})
%!assert (rmmissing ([struct, struct, struct]), [struct, struct, struct])

## Test empty input handling
%!assert (rmmissing ([]), [])
%!assert (rmmissing (ones (1, 0)), ones (1, 0))
%!assert (rmmissing (ones (1, 0), 1), ones (1, 0))
%!assert (rmmissing (ones (1, 0), 2), ones (1, 0))
%!assert (rmmissing (ones (0, 1)), ones (0, 1))
%!assert (rmmissing (ones (0, 1), 1), ones (0, 1))
%!assert (rmmissing (ones (0, 1), 2), ones (0, 1))
%!error <rmmissing: A must be a matrix; no more than 2 dimensions allowed.> ...
%!       rmmissing (ones (0, 1, 2))

## Test input validation
%!error rmmissing ()
%!error <rmmissing: A must be a matrix; no more than 2 dimensions allowed.> ...
%!       rmmissing (ones(2, 2, 2))
%!error <rmmissing: 'MinNumMissing' must be a positive integer value.> ...
%!       rmmissing (ones(2, 2), 'MinNumMissing', 0)
%!error <rmmissing: 'MinNumMissing' must be a positive integer value.> ...
%!       rmmissing ([1, 2; 3, 4], 2, "MinNumMissing", -2)
%!error <rmmissing: 'MinNumMissing' must be a positive integer value.> ...
%!       rmmissing ([1, 2; 3, 4], "MinNumMissing", 3.8)
%!error <rmmissing: 'MinNumMissing' must be a positive integer value.> ...
%!       rmmissing ([1, 2; 3, 4], "MinNumMissing", [1, 2, 3])
%!error <rmmissing: 'MinNumMissing' must be a positive integer value.> ...
%!       rmmissing ([1, 2; 3, 4], "MinNumMissing", 'xxx')
%!error <rmmissing: 'MissingLocations' must be a logical matrix of the same size as input A.> ...
%!       rmmissing ([1, 2; 3, 4], 'MissingLocations', false ([1, 1, 1]))
%!error <rmmissing: specified DIM must be either 1 or 2.> rmmissing ([1, 2; 3, 4], 5)
%!error <rmmissing: too many input arguments.> rmmissing ([1, 2; 3, 4], 'XXX', 1)
