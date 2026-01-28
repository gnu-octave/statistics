## Copyright (C) 2015 Carnë Draug <carandraug@octave.org>
## Copyright (C) 2022-2026 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {statistics} {@var{g} =} grp2idx (@var{s})
## @deftypefnx {statistics} {[@var{g}, @var{gn}] =} grp2idx (@var{s})
## @deftypefnx {statistics} {[@var{g}, @var{gn}, @var{gl}] =} grp2idx (@var{s})
##
## Get index for grouping variable.
##
## @code{@var{g} = grp2idx (@var{s})} returns a numeric column vector of integer
## values @var{g} indexing the distinct groups in the grouping variable @var{s}.
## @var{s} can specified as any of the following data types:
##
## @itemize
## @item categorical vector
## @item cell array of character vectors
## @item character array
## @item duration vector
## @item logical vector
## @item numeric vector
## @end itemize
##
## @var{s} must be a vector, unless it is a 2-D character array.  In the case of
## numerical and logical data types, the group indices are ordered in sorted
## order of @var{s}.  In the case of categorical arrays, the group indices are
## allocated by the order of the categories in @var{s}.  For the rest of the
## data types, the group indices are allocated by order of first appearance in
## @var{s}.  Note that in case of a categorical grouping variable, the indexing
## integer values might not be continuous, since @var{s} may contain unassigned
## categories.  For every other data type, @var{g} will contain integer values
## in the the range @math{[1:K]}, where @math{K} is the number of distinct
## groups in @var{s}.
##
## @code{[@var{g}, @var{gn}] = grp2idx (@var{s})} also retuns a cell array of
## character vectors @var{gn} representing the list of group names.  The order
## of the group names in @var{gn} follow the same pattern as the group indices
## in @var{g} according to the data type of @var{s}, as described above.
##
## @code{[@var{g}, @var{gn}, @var{gl}] = grp2idx (@var{s})} further returns a
## column vector @var{gl} representing the list of the group levels with the
## same data type as @var{s}.
##
## Note that standard missing values in @var{s} appear as NaN in @var{g} and are
## not present on either @var{gn} and @var{gl}.
##
## @seealso{grpstats}
## @end deftypefn

function [g, gn, gl] = grp2idx (s)

  if (nargin != 1)
    print_usage ();
  endif
  if (ndims (s) != 2)
    error ("grp2idx: S must be either a vector or a matrix.");
  endif

  is_categorical = false;
  is_char_array = false;
  #is_datetime = false;
  is_duration = false;
  is_string = false;

  if (iscategorical (s))
    if (! isvector (s))
      error ("grp2idx: 'categorical' grouping variable must be a vector.");
    endif
    is_categorical = true;
    undef = isundefined (s);
    cats = categories (s);
    s = cellstr (s);
    s(undef) = {""};
  elseif (ischar (s))
    is_char_array = true;
    s = cellstr (s);
  elseif (isdatetime (s))
    error ("grp2idx: 'datetime' grouping variable is not supported yet.");
  elseif (isduration (s))
    if (! isvector (s))
      error ("grp2idx: 'duration' grouping variable must be a vector.");
    endif
    is_duration = true;
  elseif (isstring (s))
    if (! isvector (s))
      error ("grp2idx: 'string' grouping variable must be a vector.");
    endif
    is_string = true;
    s = cellstr (s);
  elseif (! (isnumeric (s) || islogical (s) || iscellstr (s)))
    error ("grp2idx: unsupported type for input S.");
  endif

  [gl, I, g] = unique (s(:));
  ## Fix order in here, since unique does not support this yet
  if (iscellstr (s) && ! is_categorical)
    I = sort (I);
    for i = 1:length (gl)
      gl_s(i) = gl(g(I(i)));
      idx(i,:) = (g == g(I(i)));
    endfor
    for i = 1:length (gl)
      g(idx(i,:)) = i;
    endfor
    gl = gl_s;
    gl = gl';
  endif

  ## handle NaNs and empty strings
  if (iscellstr (s))
    empties = cellfun (@isempty, s);
    if (any (empties))
      g(empties) = NaN;
      rm = find (cellfun (@isempty, gl));
      to_decrement = ! isnan (g) & g > rm;
      g(to_decrement) -= 1;
    endif
    empties = cellfun (@isempty, gl);
    if (any (empties))
      gl(empties) = [];
    endif
  else
    ## This works fine because NaN come at the end after sorting, we don't
    ## have to worry about change on the indices.
    g(isnan (s)) = NaN;
    gl(isnan (gl)) = [];
  endif

  if (nargout > 1)
    if (is_categorical)
      gn = cellstr (cats);
    elseif (is_duration)
      gn = cellstr (gl);
    elseif (iscellstr (gl))
      gn = gl;
    elseif (iscell (gl))
      gn = cellfun (@num2str, gl, "UniformOutput", false);
    else
      gn = arrayfun (@num2str, gl, "UniformOutput", false);
    endif
    if (isempty (gn))
      gn = cell (0,1);
    endif
  endif

  if (nargout > 2)
    if (is_categorical)
      gl = categorical (cats);
    elseif (is_char_array)
      if (isempty (gl))
        gl = char (cell (0,1));
      else
        gl = char (gn);
      endif
    elseif (is_duration)
      if (isempty (gl))
        gl = duration (NaN (0,3));
      endif
    elseif (is_string)
      gl = string (gn);
    elseif (iscell (gl))
      if (isempty (gl))
        gl = cell (0,1);
      endif
    endif
  endif

endfunction

# test for one output argument
%!test
%! g = grp2idx ([3 2 1 2 3 1]);
%! assert (g, [3; 2; 1; 2; 3; 1]);

# test for two output arguments
%!test
%! [g, gn] = grp2idx (['b'; 'a'; 'c'; 'a']);
%! assert (g, [1; 2; 3; 2]);
%! assert (gn, {'b'; 'a'; 'c'});

## test boolean input and note that row or column vector makes no difference
%!test
%! in = [true, false, false, true];
%! out = {[2; 1; 1; 2] {'0'; '1'} [false; true]};
%! assert (nthargout (1:3, @grp2idx, in), out)
%! assert (nthargout (1:3, @grp2idx, in), nthargout (1:3, @grp2idx, in'))

## test that boolean groups are ordered in order of appearance
%!test
%! assert (nthargout (1:3, @grp2idx, [false, true]),
%!         {[1; 2] {'0'; '1'} [false; true]});
%! assert (nthargout (1:3, @grp2idx, [true, false]),
%!         {[2; 1] {'0'; '1'} [false; true]});

## test char matrix and cell array of strings
%!assert (nthargout (1:3, @grp2idx, ['oct'; 'sci'; 'oct'; 'oct'; 'sci']),
%!        {[1; 2; 1; 1; 2] {'oct'; 'sci'} ['oct'; 'sci']})
## and cell array of strings
%!assert (nthargout (1:3, @grp2idx, {'oct'; 'sci'; 'oct'; 'oct'; 'sci'}),
%!        {[1; 2; 1; 1; 2] {'oct'; 'sci'} {'oct'; 'sci'}})

## test numeric arrays
%!assert (nthargout (1:3, @grp2idx, [1, -3, -2, -3, -3,  2,  1, -1,  3, -3]),
%!        {[4; 1; 2; 1; 1; 5; 4; 3; 6; 1], {"-3"; "-2"; "-1"; "1"; "2"; "3"}, ...
%!         [-3; -2; -1; 1; 2; 3]})

%!test
%! s = [1e6, 2e6, 1e6, 3e6];
%! [g, gn, gl] = grp2idx (s);
%! assert (g, [1; 2; 1; 3]);
%! assert (gn, {'1000000'; '2000000'; '3000000'});
%! assert (gl, [1000000; 2000000; 3000000]);

%!test
%! s = [0.1, 0.2, 0.3, 0.1, 0.2];
%! [g, gn, gl] = grp2idx (s);
%! assert (g, [1; 2; 3; 1; 2]);
%! assert (gn, {'0.1'; '0.2'; '0.3'});
%! assert (gl, [0.1; 0.2; 0.3]);

%!test
%! s = [-5 -10 0 5 10 -5];
%! [g, gn, gl] = grp2idx (s);
%! assert (g, [2; 1; 3; 4; 5; 2]);
%! assert (gn, {'-10'; '-5'; '0'; '5'; '10'});
%! assert (gl, [-10; -5; 0; 5; 10]);


## test for NaN and empty strings
%!assert (nthargout (1:3, @grp2idx, [2, 2, 3, NaN, 2, 3]), ...
%!        {[1; 1; 2; NaN; 1; 2] {'2'; '3'} [2; 3]})
%!assert (nthargout (1:3, @grp2idx, {'et', 'sa', 'sa', '', 'et'}), ...
%!        {[1; 2; 2; NaN; 1] {'et'; 'sa'} {'et'; 'sa'}})

%!assert (nthargout (1:3, @grp2idx, [2, 2, 3, NaN, 2, 4]), ...
%!        {[1; 1; 2; NaN; 1; 3] {'2'; '3'; '4'} [2; 3; 4]})

%!test
%! s = [NaN, NaN, NaN];
%! [g, gn, gl] = grp2idx (s);
%! assert (g, [NaN; NaN; NaN]);
%! assert (gn, cell (0,1));
%! assert (gl, zeros (0,1));

%!test
%! s = single ([NaN, NaN, NaN]);
%! [g, gn, gl] = grp2idx (s);
%! assert (g, [NaN; NaN; NaN]);
%! assert (gn, cell (0,1));
%! assert (gl, single (zeros(0,1)));

%!test
%! s = {''; ''; ''; ''};
%! [g, gn, gl] = grp2idx (s);
%! assert (g, [NaN; NaN; NaN; NaN]);
%! assert (gn, cell(0,1));
%! assert (gl, cell(0,1));

%!test
%! s = {'a'; ''; 'b'; ''; 'c'};
%! [g, gn, gl] = grp2idx (s);
%! assert (g, [1; NaN; 2; NaN; 3]);
%! assert (gn, {'a'; 'b'; 'c'});
%! assert (gl, {'a'; 'b'; 'c'});

%!test
%! s = categorical ({''; ''; ''; ''});
%! [g, gn, gl] = grp2idx (s);
%! assert (g, [NaN; NaN; NaN; NaN]);
%! assert (gn, cell (0,1));
%! assert (isequaln (gl, categorical (cell (0,1))));

%!test
%! s = string ({missing, missing, missing});
%! [g, gn, gl] = grp2idx (s);
%! assert (g, [NaN; NaN; NaN]);
%! assert (gn, cell (0,1));
%! assert (isequal (gl, string (cell (0,1))));

%!test
%! s = [duration(NaN, 0, 0), duration(NaN, 0, 0), duration(NaN, 0, 0)];
%! [g, gn, gl] = grp2idx (s);
%! assert (g, [NaN; NaN; NaN]);
%! assert (gn, cell (0,1));
%! assert (isequal (gl, duration (NaN (0,3))));

## Test that order when handling strings is by order of appearance
%!test assert (nthargout (1:3, @grp2idx, ['sci'; 'oct'; 'sci'; 'oct'; 'oct']),
%!        {[1; 2; 1; 2; 2] {'sci'; 'oct'} ['sci'; 'oct']});
%!test assert (nthargout (1:3, @grp2idx, {'sci'; 'oct'; 'sci'; 'oct'; 'oct'}),
%!        {[1; 2; 1; 2; 2] {'sci'; 'oct'} {'sci'; 'oct'}});
%!test assert (nthargout (1:3, @grp2idx, {'sa' 'et' 'et' '' 'sa'}),
%!        {[1; 2; 2; NaN; 1] {'sa'; 'et'} {'sa'; 'et'}})

## test for categorical arrays
%!test
%! [g, gn, gl] = grp2idx (categorical ({'low', 'med', 'high', 'low'}));
%! assert (g, [2; 3; 1; 2]);
%! assert (gn, {'high'; 'low'; 'med'});
%! assert (isequal (gl, categorical ({'high'; 'low'; 'med'})));

%!test
%! [g, gn, gl] = grp2idx (categorical ([10, 20, 10, 30, 20]));
%! assert (g, [1; 2; 1; 3; 2]);
%! assert (gn, {'10'; '20'; '30'});
%! assert (isequal (gl, categorical ([10; 20; 30])));

%!test
%! cats = categorical ({'high', '<undefined>', 'low', '<undefined>'});
%! [g, gn, gl] = grp2idx (cats);
%! assert (g, [2; 1; 3; 1]);
%! assert (gn, {'<undefined>'; 'high'; 'low'});
%! assert (isequal (gl, categorical ({'<undefined>'; 'high'; 'low'})));

%!test
%! s = categorical ({''; ''; ''; ''}, {'1', '2', '3'}, {'1', '2', '3'});
%! [g, gn, gl] = grp2idx (s);
%! assert (g, nan (4, 1));
%! assert (gn, {'1'; '2'; '3'});
%! assert (iscategorical (gl), true);
%! assert (cellstr (gl), gn);

%!test
%! s = categorical ({''; '1'; ''; '2'}, {'1', '2', '3'}, {'1', '2', '3'});
%! [g, gn, gl] = grp2idx (s);
%! assert (g, [NaN; 1; NaN; 2]);
%! assert (gn, {'1'; '2'; '3'});
%! assert (iscategorical (gl), true);
%! assert (cellstr (gl), gn);

%!test
%! s = categorical ({''; '2'; ''; '1'}, {'1', '2', '3'}, {'1', '2', '3'});
%! [g, gn, gl] = grp2idx (s);
%! assert (g, [NaN; 2; NaN; 1]);
%! assert (gn, {'1'; '2'; '3'});
%! assert (iscategorical (gl), true);
%! assert (cellstr (gl), gn);

## test for duration arrays
%!test
%! g = gn = gl = [];
%! [g, gn, gl] = grp2idx (seconds ([1.234, 1.234, 2.5, 3.000]));
%! assert (g, [1; 1; 2; 3]);
%! assert (gn, {'1.234 sec'; '2.5 sec'; '3 sec'});
%! assert (isequal (gl, seconds ([1.234; 2.5; 3.000])));

%!test
%! [g, gn, gl] = grp2idx ([hours(1); hours(2); hours(1); hours(3)]);
%! assert (g, [1; 2; 1; 3]);
%! assert (gn, {'1 hr'; '2 hr'; '3 hr'});
%! assert (isequal (gl, [hours(1); hours(2); hours(3)]));

%!test
%! in = [duration(1, 30, 0); duration(0, 45, 30); duration(1, 30, 0); duration(2, 15, 15)];
%! [g, gn, gl] = grp2idx(in);
%! assert (g, [2; 1; 2; 3]);
%! assert (gn, {'00:45:30'; '01:30:00'; '02:15:15'});
%! assert (isequal (gl, [duration(0, 45, 30); duration(1, 30, 0); duration(2, 15, 15)]));

## Inconsistency Note: following test is inconsistent with MATLAB due to a
## probable bug in their implementation, where they include multiple NaNs
## in the output group labels for duration array inputs.
%!test
%! in = [hours(1); NaN; minutes(30); hours(1); NaN; seconds(90)];
%! [g, gn, gl] = grp2idx(in);
%! assert (g, [3; NaN; 2; 3; NaN; 1]);
%! assert (gn, {'0.025 hr'; '0.5 hr'; '1 hr'});
%! assert (isequal(gl, [seconds(90); minutes(30);  hours(1)]));

## test for string arrays
%!test
%! [g, gn, gl] = grp2idx (string ({'123', 'erw', missing, '', '234'}));
%! assert (g, [1; 2; NaN; NaN; 3]);
%! assert (gn, {'123'; 'erw'; '234'});
%! assert (isequal (gl, string ({'123'; 'erw'; '234'})));

%!test
%! [g, gn, gl] = grp2idx (string ({'medium', 'low', 'high', 'medium', 'medium'}));
%! assert (g, [1; 2; 3; 1; 1]);
%! assert (gn, {'medium'; 'low'; 'high'});
%! assert (isequal (gl, string ({'medium'; 'low'; 'high'})));

%!test
%! [g, gn, gl] = grp2idx (string ({'', 'high', 'low', ''}));
%! assert (g, [NaN; 1; 2; NaN]);
%! assert (gn, {'high'; 'low'});
%! assert (isequal (gl, string ({'high'; 'low'})));

%!test
%! [g, gn, gl] = grp2idx (string ({'a', 'a', 'b', 'c'}));
%! assert (g, [1; 1; 2; 3]);
%! assert (gn, {'a'; 'b'; 'c'});
%! assert (isstring (gl), true);
%! assert (cellstr (gl), gn);

## Test input validation
%!error <grp2idx: S must be either a vector or a matrix.> grp2idx (ones (3, 3, 3))
%!error <grp2idx: 'categorical' grouping variable must be a vector.> ...
%! grp2idx (categorical ([1, 2; 1, 3]))
%!error <grp2idx: 'datetime' grouping variable is not supported yet.> ...
%! grp2idx (datetime ('now'))
%!error <grp2idx: 'duration' grouping variable must be a vector.> ...
%! grp2idx (days ([1, 2; 1, 3]))
%!error <grp2idx: 'string' grouping variable must be a vector.> ...
%! grp2idx (string ({'a', 'a'; 'b', 'c'}))
%!error <grp2idx: unsupported type for input S.> grp2idx ({1})
