## Copyright (C) 2015 Carnë Draug <carandraug@octave.org>
## Copyright (C) 2022 Andreas Bertsatos <abertsatos@biol.uoa.gr>
##
## This file is part of the statistics package for GNU Octave.
##
## This program is free software; you can redistribute it and/or
## modify it under the terms of the GNU General Public License as
## published by the Free Software Foundation; either version 3 of the
## License, or (at your option) any later version.
##
## This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
## General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, see
## <http:##www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn  {statistics} {[@var{g}, @var{gn}, @var{gl}] =} grp2idx (@var{s})
##
## Get index for group variables.
##
## For variable @var{s}, returns the indices @var{g}, into the variable
## groups @var{gn} and @var{gl}.  The first has a string representation of
## the groups while the later has its actual values.  In the case of numerical
## and logical data types, the group indices are ordered in sorted order of
## @var{s}.  In the case of categorical arrays, the group indices are allocated
## by the order of the categories in @var{s}.  For the rest of the data types,
## the group indices are allocated by order of first appearance in @var{s}.
##
## NaNs and empty strings in @var{s} appear as NaN in @var{g} and are
## not present on either @var{gn} and @var{gl}.
##
## @seealso{grpstats}
## @end deftypefn

function [g, gn, gl] = grp2idx (s)
  if (nargin != 1)
    print_usage ();
  endif

  s_was_char = false;
  s_was_categorical = false;
  s_was_duration = false;
  s_was_string = false;
  if (ischar (s))
    s_was_char = true;
    s = cellstr (s);
  elseif (isstring (s) && isvector (s))
    s_was_string = true;
    s = cellstr (s);
  elseif (iscategorical (s) && isvector (s))
    s_was_categorical = true;
    undef = isundefined (s);
    cats = categories (s);
    s = cellstr (string (s));
    s(undef) = {''};
  elseif (isduration (s) && isvector (s))
    s_was_duration = true;
  elseif (! isvector (s))
    error ("grp2idx: S must be a vector, cell array of strings, or char matrix");
  endif

  [gl, I, g] = unique (s(:));
  ## Fix order in here, since unique does not support this yet
  if (iscellstr (s) && ! s_was_categorical)
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
    if (s_was_categorical)
      gn = categories (categorical (s));
    elseif (s_was_duration)
      gn = cellstr(gl);
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

  if (nargout > 2 && s_was_char)
    if (isempty (gl))
      gl = char (cell (0,1));
    else
      gl = char (gn);
    endif
  elseif (nargout > 2 && s_was_categorical)
    gl = categorical (gn);
  elseif (nargout > 2 && s_was_string)
    gl = string (gn);
  elseif (nargout > 2 && s_was_duration)
    if (isempty (gl))
      gl = duration (NaN (0,3));
    endif
  elseif (nargout > 2 && iscell (gl))
    if (isempty (gl))
      gl = cell (0,1);
    endif
  endif

endfunction

# test for one output argument
%!test
%! g = grp2idx ([3 2 1 2 3 1]);
%! assert (isequal (g, [3; 2; 1; 2; 3; 1]));
# test for two output arguments
%!test
%! [g, gn] = grp2idx (['b'; 'a'; 'c'; 'a']);
%! assert (isequal (g, [1; 2; 3; 2]));
%! assert (isequal (gn, {'b'; 'a'; 'c'}));
## test boolean input and note that row or column vector makes no difference
%!test
%! in = [true false false true];
%! out = {[2; 1; 1; 2] {'0'; '1'} [false; true]};
%! assert (nthargout (1:3, @grp2idx, in), out)
%! assert (nthargout (1:3, @grp2idx, in), nthargout (1:3, @grp2idx, in'))
## test that boolean groups are ordered in order of appearance
%!test
%! assert (nthargout (1:3, @grp2idx, [false, true]), {[1; 2] {'0'; '1'} ...
%!                                                   [false; true]});
%! assert (nthargout (1:3, @grp2idx, [true, false]), {[2; 1] {'0'; '1'} ...
%!                                                   [false; true]});
## test char matrix and cell array of strings
%!assert (nthargout (1:3, @grp2idx, ['oct'; 'sci'; 'oct'; 'oct'; 'sci']),
%!        {[1; 2; 1; 1; 2] {'oct'; 'sci'} ['oct'; 'sci']});
## and cell array of strings
%!assert (nthargout (1:3, @grp2idx, {'oct'; 'sci'; 'oct'; 'oct'; 'sci'}),
%!        {[1; 2; 1; 1; 2] {'oct'; 'sci'} {'oct'; 'sci'}});
## test numeric arrays
%!assert (nthargout (1:3, @grp2idx, [ 1 -3 -2 -3 -3  2  1 -1  3 -3]),
%!        {[4; 1; 2; 1; 1; 5; 4; 3; 6; 1], {'-3'; '-2'; '-1'; '1'; '2'; '3'}, ...
%!         [-3; -2; -1; 1; 2; 3]});
%!test
%! s = [1e6 2e6 1e6 3e6];
%! [g, gn, gl] = grp2idx (s);
%! assert (g, [1; 2; 1; 3]);
%! assert (gn, {'1000000'; '2000000'; '3000000'});
%! assert (gl, [1000000; 2000000; 3000000]);
%!test
%! s = [0.1 0.2 0.3 0.1 0.2];
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
%!assert (nthargout (1:3, @grp2idx, [2 2 3 NaN 2 3]),
%!        {[1; 1; 2; NaN; 1; 2] {'2'; '3'} [2; 3]})
%!assert (nthargout (1:3, @grp2idx, {'et' 'sa' 'sa' '' 'et'}),
%!        {[1; 2; 2; NaN; 1] {'et'; 'sa'} {'et'; 'sa'}})
%!assert (nthargout (1:3, @grp2idx, [2 2 3 NaN 2 4]),
%!        {[1; 1; 2; NaN; 1; 3] {'2'; '3'; '4'} [2; 3; 4]})
%!test
%! s = [NaN, NaN, NaN];
%! [g, gn, gl] = grp2idx (s);
%! assert (isequaln (g, [NaN; NaN; NaN]));
%! assert (isequaln (gn, cell (0,1)));
%! assert (isequaln (gl, zeros (0,1)));
%!test
%! s = single ([NaN, NaN, NaN]);
%! [g, gn, gl] = grp2idx (s);
%! assert (isequaln (g, [NaN; NaN; NaN]));
%! assert (isequaln (gn, cell (0,1)));
%! assert (isequaln (gl, single (zeros (0,1))));
%!test
%! s = {''; ''; ''; ''};
%! [g, gn, gl] = grp2idx (s);
%! assert (isequaln (g, [NaN; NaN; NaN; NaN]), true);
%! assert (isequaln (gn, cell (0,1)));
%! assert (isequaln (gl, cell (0,1)));
%!test
%! s = {'a'; ''; 'b'; ''; 'c'};
%! [g, gn, gl] = grp2idx (s);
%! assert (isequaln (g, [1; NaN; 2; NaN; 3]), true);
%! assert (gn, {'a'; 'b'; 'c'});
%! assert (gl, {'a'; 'b'; 'c'});
%!test
%! s = categorical ({''; ''; ''; ''});
%! [g, gn, gl] = grp2idx (s);
%! assert (isequaln (g, [NaN; NaN; NaN; NaN]), true);
%! assert (isequaln (gn, cell (0,1)));
%! assert (isequaln (gl, categorical (cell (0,1))));
%!test
%! s = string ({missing, missing, missing});
%! out = string (cell (0,1));
%! [g, gn, gl] = grp2idx (s);
%! assert (isequaln (g, [NaN; NaN; NaN]), true);
%! assert (isequal (gn, cell (0,1)));
%!test
%! s = [duration(NaN,0,0), duration(NaN,0,0), duration(NaN,0,0)];
%! [g, gn, gl] = grp2idx (s);
%! assert (isequaln (g, [NaN; NaN; NaN]));
%! assert (isequal (gn, cell (0,1)));
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
%! assert (isequal (g, [2; 3; 1; 2]));
%! assert (isequal (gn, {'high'; 'low'; 'med'}));
%! assert (isequal (gl, categorical ({'high'; 'low'; 'med'})));
%!test
%! [g, gn, gl] = grp2idx (categorical ([10, 20, 10, 30, 20]));
%! assert (isequal (g, [1; 2; 1; 3; 2]));
%! assert (isequal (gn, {'10'; '20'; '30'}));
%! assert (isequal (gl, categorical ([10; 20; 30])));
%!test
%! cats = categorical ({'high', '<undefined>', 'low', '<undefined>'});
%! [g, gn, gl] = grp2idx (cats);
%! assert (isequal (g, [2; 1; 3; 1]));
%! assert (isequal (gn, {'<undefined>'; 'high'; 'low'}));
%! assert (isequal (gl, categorical ({'<undefined>'; 'high'; 'low'})));
## test for duration arrays
%!test
%! g = gn = gl = [];
%! [g, gn, gl] = grp2idx (seconds ([1.234, 1.234, 2.5, 3.000]));
%! assert (isequal (g, [1; 1; 2; 3]));
%! assert (isequal (gn, {'1.234 sec'; '2.5 sec'; '3 sec'}));
%! assert (isequal (gl, seconds ([1.234; 2.5; 3.000])));
%!test
%! [g, gn, gl] = grp2idx ([hours(1); hours(2); hours(1); hours(3)]);
%! assert(isequal (g, [1; 2; 1; 3]));
%! assert(isequal (gn, {'1 hr'; '2 hr'; '3 hr'}));
%! assert(isequal (gl, [hours(1); hours(2); hours(3)]));
%!test
%! in = [duration(1, 30, 0); duration(0, 45, 30); duration(1, 30, 0); duration(2, 15, 15)];
%! [g, gn, gl] = grp2idx (in);
%! assert(isequal (g, [2; 1; 2; 3]));
%! assert(isequal (gn, {'00:45:30'; '01:30:00'; '02:15:15'}));
%! assert(isequal (gl, [duration(0,45,30); duration(1,30,0); duration(2,15,15)]));
## Inconsistency Note: following test is inconsistent with MATLAB due to a
## probable bug in their implementation, where they include multiple NaNs
## in the output group labels for duration array inputs.
%!test
%! in = [hours(1); NaN; minutes(30); hours(1); NaN; seconds(90)];
%! [g, gn, gl] = grp2idx (in);
%! assert(isequaln (g, [3; NaN; 2; 3; NaN; 1]));
%! assert(isequal (gn, {'0.025 hr'; '0.5 hr'; '1 hr'}));
%! assert(isequal (gl, [seconds(90); minutes(30);  hours(1)]));
## test for string arrays
%!test
%! [g, gn, gl] = grp2idx (string ({'123', 'erw', missing, '', '234'}));
%! assert (isequaln (g, [1; 2; NaN; NaN; 3]));
%! assert (isequal (gn, {'123'; 'erw'; '234'}));
%! assert (isequal (gl, string ({'123'; 'erw'; '234'})));
%!test
%! [g, gn, gl] = grp2idx (string ({'medium', 'low', 'high', 'medium', 'medium'}));
%! assert (isequaln (g, [1; 2; 3; 1; 1]));
%! assert (isequal (gn, {'medium'; 'low'; 'high'}));
%! assert (isequal (gl, string ({'medium'; 'low'; 'high'})));
%!test
%! [g, gn, gl] = grp2idx (string ({'', 'high', 'low', ''}));
%! assert (isequaln (g, [NaN; 1; 2; NaN]));
%! assert (isequal (gn, {'high'; 'low'}));
%! assert (isequal (gl, string ({'high'; 'low'})));
%!test
%! s = [true false false true];
%! [g, gn, gl] = grp2idx (s);
%! assert (g, [2; 1; 1; 2]);
%! assert (gn, {'0'; '1'});
%! assert (gl, [false; true]);
%!test
%! s = [true; false; false; true];
%! [g, gn, gl] = grp2idx (s);
%! assert (g, [2; 1; 1; 2]);
%! assert (gn, {'0'; '1'});
%! assert (gl, [false; true]);
%!test
%! s = [false true];
%! [g, gn, gl] = grp2idx (s);
%! assert (g, [1; 2]);
%! assert (gn, {'0'; '1'});
%! assert (gl, [false; true]);
%!test
%! s = [true false];
%! [g, gn, gl] = grp2idx (s);
%! assert (g, [2; 1]);
%! assert (gn, {'0'; '1'});
%! assert (gl, [false; true]);
%!test
%! s = [true true true];
%! [g, gn, gl] = grp2idx (s);
%! assert (g, [1; 1; 1]);
%! assert (gn, {'1'});
%! assert (gl, true);
%!test
%! s = [false false false];
%! [g, gn, gl] = grp2idx (s);
%! assert (g, [1; 1; 1]);
%! assert (gn, {'0'});
%! assert (gl, false);
%!test
%! s = ['oct'; 'sci'; 'oct'; 'oct'; 'sci'];
%! [g, gn, gl] = grp2idx (s);
%! assert (g, [1; 2; 1; 1; 2]);
%! assert (gn, {'oct'; 'sci'});
%! assert (gl, ['oct'; 'sci']);
%!test
%! s = ['sci'; 'oct'; 'sci'; 'oct'; 'oct'];
%! [g, gn, gl] = grp2idx (s);
%! assert (g, [1; 2; 1; 2; 2]);
%! assert (gn, {'sci'; 'oct'});
%! assert (gl, ['sci'; 'oct']);
%!test
%! s = {'oct'; 'sci'; 'oct'; 'oct'; 'sci'};
%! [g, gn, gl] = grp2idx (s);
%! assert (g, [1; 2; 1; 1; 2]);
%! assert (gn, {'oct'; 'sci'});
%! assert (gl, {'oct'; 'sci'});
%!test
%! s = {'sci'; 'oct'; 'sci'; 'oct'; 'oct'};
%! [g, gn, gl] = grp2idx (s);
%! assert (g, [1; 2; 1; 2; 2]);
%! assert (gn, {'sci'; 'oct'});
%! assert (gl, {'sci'; 'oct'});
%!test
%! s = {'et' 'sa' 'sa' '' 'et'};
%! [g, gn, gl] = grp2idx (s);
%! assert (isequaln (g, [1; 2; 2; NaN; 1]), true);
%! assert (gn, {'et'; 'sa'});
%! assert (gl, {'et'; 'sa'});
%!test
%! s = {'sa' 'et' 'et' '' 'sa'};
%! [g, gn, gl] = grp2idx (s);
%! assert (isequaln (g, [1; 2; 2; NaN; 1]), true);
%! assert (gn, {'sa'; 'et'});
%! assert (gl, {'sa'; 'et'});
%!test
%! s = {'a'; 'b'; 'c'; 'a'; 'b'; 'c'};
%! [g, gn, gl] = grp2idx (s);
%! assert (g, [1; 2; 3; 1; 2; 3]);
%! assert (gn, {'a'; 'b'; 'c'});
%! assert (gl, {'a'; 'b'; 'c'});
%!test
%! s = {'c'; 'b'; 'a'; 'c'; 'b'; 'a'};
%! [g, gn, gl] = grp2idx (s);
%! assert (g, [1; 2; 3; 1; 2; 3]);
%! assert (gn, {'c'; 'b'; 'a'});
%! assert (gl, {'c'; 'b'; 'a'});
%!test
%! s = {''; ''; ''; ''};
%! [g, gn, gl] = grp2idx (s);
%! assert (isequaln (g, [NaN; NaN; NaN; NaN]), true);
%! assert (isequal (size (gn), [0 1]))
%! assert (isequal (size (gl), [0 1]))
%!test
%! s = {'a'; ''; 'b'; ''; 'c'};
%! [g, gn, gl] = grp2idx (s);
%! assert (isequaln (g, [1; NaN; 2; NaN; 3]), true);
%! assert (gn, {'a'; 'b'; 'c'});
%! assert (gl, {'a'; 'b'; 'c'});
%!test
%! s = [1 -3 -2 -3 -3 2 1 -1 3 -3];
%! [g, gn, gl] = grp2idx (s);
%! assert (g, [4; 1; 2; 1; 1; 5; 4; 3; 6; 1]);
%! assert (gn, {'-3'; '-2'; '-1'; '1'; '2'; '3'});
%! assert (gl, [-3; -2; -1; 1; 2; 3]);
%!test
%! s = [1 2 3 4 5];
%! [g, gn, gl] = grp2idx (s);
%! assert (g, [1; 2; 3; 4; 5]);
%! assert (gn, {'1'; '2'; '3'; '4'; '5'});
%! assert (gl, [1; 2; 3; 4; 5]);
%!test
%! s = [5 4 3 2 1];
%! [g, gn, gl] = grp2idx (s);
%! assert (g, [5; 4; 3; 2; 1]);
%! assert (gn, {'1'; '2'; '3'; '4'; '5'});
%! assert (gl, [1; 2; 3; 4; 5]);
%!test
%! s = [1 1 1 1 1];
%! [g, gn, gl] = grp2idx (s);
%! assert (g, [1; 1; 1; 1; 1]);
%! assert (gn, {'1'});
%! assert (gl, 1);
%!test
%! s = [10; 20; 30; 10; 20];
%! [g, gn, gl] = grp2idx (s);
%! assert (g, [1; 2; 3; 1; 2]);
%! assert (gn, {'10'; '20'; '30'});
%! assert (gl, [10; 20; 30]);
%!test
%! s = [2 2 3 NaN 2 3];
%! [g, gn, gl] = grp2idx (s);
%! assert (isequaln (g, [1; 1; 2; NaN; 1; 2]), true);
%! assert (gn, {'2'; '3'});
%! assert (gl, [2; 3]);
%!test
%! s = [1 NaN 2 NaN 3];
%! [g, gn, gl] = grp2idx (s);
%! assert (isequaln (g, [1; NaN; 2; NaN; 3]), true);
%! assert (gn, {'1'; '2'; '3'});
%! assert (gl, [1; 2; 3]);
%!test
%! s = [NaN NaN NaN];
%! [g, gn, gl] = grp2idx (s);
%! assert (isequaln (g, [NaN; NaN; NaN]), true);
%! assert (isequal (size (gn), [0 1]))
%! assert (gl, zeros (0, 1));
%!test
%! s = [1.5 2.5 1.5 3.5 2.5];
%! [g, gn, gl] = grp2idx (s);
%! assert (g, [1; 2; 1; 3; 2]);
%! assert (gn, {'1.5'; '2.5'; '3.5'});
%! assert (gl, [1.5; 2.5; 3.5]);
%!test
%! s = [-5 -10 0 5 10 -5];
%! [g, gn, gl] = grp2idx (s);
%! assert (g, [2; 1; 3; 4; 5; 2]);
%! assert (gn, {'-10'; '-5'; '0'; '5'; '10'});
%! assert (gl, [-10; -5; 0; 5; 10]);
%!test
%! s = [0 0 0 0];
%! [g, gn, gl] = grp2idx (s);
%! assert (g, [1; 1; 1; 1]);
%! assert (gn, {'0'});
%! assert (gl, 0);
%!test
%! s = [1e6 2e6 1e6 3e6];
%! [g, gn, gl] = grp2idx (s);
%! assert (g, [1; 2; 1; 3]);
%! assert (gn, {'1000000'; '2000000'; '3000000'});
%! assert (gl, [1000000; 2000000; 3000000]);
%!test
%! s = categorical ({'low', 'med', 'high', 'low'});
%! [g, gn, gl] = grp2idx (s);
%! assert (isequal (g, [2; 3; 1; 2]));
%! assert (isequal (gn, {'high'; 'low'; 'med'}));
%! assert (isequal (gl, categorical ({'high'; 'low'; 'med'})));
%!test
%! s = categorical ([10, 20, 10, 30, 20]);
%! [g, gn, gl] = grp2idx (s);
%! assert (isequal (g, [1; 2; 1; 3; 2]));
%! assert (isequal (gn, {'10'; '20'; '30'}));
%! assert (isequal (gl, categorical ({'10'; '20'; '30'})));
%!test
%! s = categorical ({'high', '', 'low', ''});
%! [g, gn, gl] = grp2idx (s);
%! assert (isequaln (g, [1; NaN; 2; NaN]), true);
%! assert (isequal (gn, {'high'; 'low'}));
%! assert (isequal (gl, categorical ({'high'; 'low'})));
%!test
%! s = categorical ({'a', 'b', 'c', 'a', 'b', 'c'});
%! [g, gn, gl] = grp2idx (s);
%! assert (isequal (g, [1; 2; 3; 1; 2; 3]));
%! assert (isequal (gn, {'a'; 'b'; 'c'}));
%! assert (isequal (gl, categorical ({'a'; 'b'; 'c'})));
%!test
%! s = seconds ([1.234, 1.234, 2.5, 3.000]);
%! [g, gn, gl] = grp2idx (s);
%! assert (isequal (g, [1; 1; 2; 3]));
%! assert (isequal (gn, {'1.234 sec'; '2.5 sec'; '3 sec'}));
%! assert (isequal (gl, seconds ([1.234; 2.5; 3])));
%!test
%! s = [hours(1); hours(2); hours(1); hours(3)];
%! [g, gn, gl] = grp2idx (s);
%! assert (isequal (g, [1; 2; 1; 3]));
%! assert (isequal (gn, {'1 hr'; '2 hr'; '3 hr'}));
%! assert (isequal (gl, [hours(1); hours(2); hours(3)]));
%!test
%! s = [duration(1, 30, 0); duration(0, 45, 30); duration(1, 30, 0); duration(2, 15, 15)];
%! [g, gn, gl] = grp2idx (s);
%! assert (isequal (g, [2; 1; 2; 3]));
%! assert (isequal (gn, {'00:45:30'; '01:30:00'; '02:15:15'}));
%! assert (isequal (gl, [duration(0,45,30); duration(1,30,0); duration(2,15,15)]));
%!test
%! s = [minutes(10); minutes(20); minutes(10); minutes(30)];
%! [g, gn, gl] = grp2idx (s);
%! assert (isequal (g, [1; 2; 1; 3]));
%! assert (isequal (gn, {'10 min'; '20 min'; '30 min'}));
%! assert (isequal (gl, [minutes(10); minutes(20); minutes(30)]));
%!test
%! s = string ({'123', 'erw', missing, '', '234'});
%! [g, gn, gl] = grp2idx (s);
%! assert (isequaln (g, [1; 2; NaN; NaN; 3]), true);
%! assert (isequal (gn, {'123'; 'erw'; '234'}));
%! assert (isequal (gl, string ({'123'; 'erw'; '234'})));
%!test
%! s = string ({'medium', 'low', 'high', 'medium', 'medium'});
%! [g, gn, gl] = grp2idx (s);
%! assert (isequaln (g, [1; 2; 3; 1; 1]), true);
%! assert (isequal (gn, {'medium'; 'low'; 'high'}));
%! assert (isequal (gl, string ({'medium'; 'low'; 'high'})));
%!test
%! s = string ({'', 'high', 'low', ''});
%! [g, gn, gl] = grp2idx (s);
%! assert (isequaln (g, [NaN; 1; 2; NaN]), true);
%! assert (isequal (gn, {'high'; 'low'}));
%! assert (isequal (gl, string ({'high'; 'low'})));
%!test
%! s = string ({'a', 'b', 'c', 'a', 'b'});
%! [g, gn, gl] = grp2idx (s);
%! assert (isequaln (g, [1; 2; 3; 1; 2]), true);
%! assert (isequal (gn, {'a'; 'b'; 'c'}));
%! assert (isequal (gl, string ({'a'; 'b'; 'c'})));
%!test
%! s = string ({'test', 'test', 'test'});
%! [g, gn, gl] = grp2idx (s);
%! assert (isequal (g, [1; 1; 1]));
%! assert (isequal (gn, {'test'}));
%! assert (isequal (gl, string ({'test'})));
%!test
%! s = string ({missing, missing, missing});
%! [g, gn, gl] = grp2idx (s);
%! assert (isequaln (g, [NaN; NaN; NaN]), true);
%! assert (isequal (size (gn), [0 1]))
%! assert (isequal (size (gl), [0 1]));
%!test
%! s = [1 2 3];
%! g = grp2idx (s);
%! assert (g, [1; 2; 3]);
%!test
%! s = {'a'; 'b'; 'c'};
%! g = grp2idx (s);
%! assert (g, [1; 2; 3]);
%!test
%! s = [true false true];
%! g = grp2idx (s);
%! assert (g, [2; 1; 2]);
%!test
%! s = [1 2 3 4 5];
%! [g, gn] = grp2idx (s);
%! assert (g, [1; 2; 3; 4; 5]);
%! assert (gn, {'1'; '2'; '3'; '4'; '5'});
%!test
%! s = {'x'; 'y'; 'z'};
%! [g, gn] = grp2idx (s);
%! assert (g, [1; 2; 3]);
%! assert (gn, {'x'; 'y'; 'z'});
%!test
%! s = categorical ({'cat1', 'cat2', 'cat3'});
%! [g, gn] = grp2idx (s);
%! assert (g, [1; 2; 3]);
%! assert (gn, {'cat1'; 'cat2'; 'cat3'});
%!test
%! s = [1; 1; 2; 2; 3; 3; 4; 4; 5; 5];
%! [g, gn, gl] = grp2idx (s);
%! assert (g, [1; 1; 2; 2; 3; 3; 4; 4; 5; 5]);
%! assert (gn, {'1'; '2'; '3'; '4'; '5'});
%! assert (gl, [1; 2; 3; 4; 5]);
%!test
%! s = {'group1'; 'group2'; 'group1'; 'group3'; 'group2'; 'group3'};
%! [g, gn, gl] = grp2idx (s);
%! assert (g, [1; 2; 1; 3; 2; 3]);
%! assert (gn, {'group1'; 'group2'; 'group3'});
%! assert (gl, {'group1'; 'group2'; 'group3'});
%!test
%! s = [10 20 30 40 50];
%! [g, gn, gl] = grp2idx (s);
%! assert (g, [1; 2; 3; 4; 5]);
%! assert (gn, {'10'; '20'; '30'; '40'; '50'});
%! assert (gl, [10; 20; 30; 40; 50]);
%!test
%! s = logical ([1 0 1 0 1 0]);
%! [g, gn, gl] = grp2idx (s);
%! assert (g, [2; 1; 2; 1; 2; 1]);
%! assert (gn, {'0'; '1'});
%! assert (gl, [false; true]);
%!test
%! s = [0.1 0.2 0.3 0.1 0.2];
%! [g, gn, gl] = grp2idx (s);
%! assert (g, [1; 2; 3; 1; 2]);
%! assert (gn, {'0.1'; '0.2'; '0.3'});
%! assert (gl, [0.1; 0.2; 0.3]);