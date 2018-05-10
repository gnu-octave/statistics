## Copyright (C) 2015 Carnë Draug <carandraug@octave.org>
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
## @deftypefn {Function File} {[@var{g}, @var{gn}, @var{gl}] =} grp2idx (@var{s})
## Get index for group variables.
##
## For variable @var{s}, returns the indices @var{g}, into the variable
## groups @var{gn} and @var{gl}.  The first has a string representation of
## the groups while the later has its actual values.
##
## NaNs and empty strings in @var{s} appear as NaN in @var{g} and are
## not present on either @var{gn} and @var{gl}.
##
## @seealso{cellstr, num2str, unique}
## @end deftypefn

function [g, gn, gl] = grp2idx (s)
  if (nargin != 1)
    print_usage ();
  endif

  s_was_char = false;
  if (ischar (s))
    s_was_char = true;
    s = cellstr (s);
  elseif (! isvector (s))
    error ("grp2idx: S must be a vector, cell array of strings, or char matrix");
  endif

  ## FIXME once Octave core implements "sorted" and "stable" argument to
  ##        unique(), we can use the following snippet so that we are fully
  ##        Matlab compatible.
#  set_order = "sorted";
#  if (iscellstr (s))
#    set_order = "stable";
#  endif
#  [gl, ~, g] = unique (s(:), set_order);

  [gl, ~, g] = unique (s(:));

  ## handle NaNs and empty strings
  if (iscellstr (s))
    ## FIXME empty strings appear at the front because unique is sorting
    ##        them, so we only need to subtract one.  However, when fix the
    ##        order for strings (when core's unique has the stable option),
    ##        then we'll have to come up with something clever.
    empties = cellfun (@isempty, s);
    if (any (empties))
      g(empties) = NaN;
      g--;
      gl(1) = [];
    endif
  else
    ## This works fine because NaN come at the end after sorting, we don't
    ## have to worry about change on the indices.
    g(isnan (s)) = NaN;
    gl(isnan (gl)) = [];
  endif

  if (isargout (2))
    if (iscellstr (gl))
      gn = gl;
    elseif (iscell (gl))
      gn = cellfun (@num2str, gl, "UniformOutput", false);
    else
      gn = arrayfun (@num2str, gl, "UniformOutput", false);
    endif
  endif

  if (isargout (3) && s_was_char)
    gl = char (gl);
  endif

endfunction

## test boolean input and note that row or column vector makes no difference
%!test
%! in = [true false false true];
%! out = {[2; 1; 1; 2] {"0"; "1"} [false; true]};
%! assert (nthargout (1:3, @grp2idx, in), out)
%! assert (nthargout (1:3, @grp2idx, in), nthargout (1:3, @grp2idx, in'))

## test that groups are ordered in boolean
%!test
%! assert (nthargout (1:3, @grp2idx, [false true]),
%!         {[1; 2] {"0"; "1"} [false; true]});
%! assert (nthargout (1:3, @grp2idx, [true false]),
%!         {[2; 1] {"0"; "1"} [false; true]});

## test char matrix and cell array of strings
%!assert (nthargout (1:3, @grp2idx, ["oct"; "sci"; "oct"; "oct"; "sci"]),
%!        {[1; 2; 1; 1; 2] {"oct"; "sci"} ["oct"; "sci"]});
## and cell array of strings
%!assert (nthargout (1:3, @grp2idx, {"oct"; "sci"; "oct"; "oct"; "sci"}),
%!        {[1; 2; 1; 1; 2] {"oct"; "sci"} {"oct"; "sci"}});

## test numeric arrays
%!assert (nthargout (1:3, @grp2idx, [ 1 -3 -2 -3 -3  2  1 -1  3 -3]),
%!        {[4; 1; 2; 1; 1; 5; 4; 3; 6; 1] {"-3"; "-2"; "-1"; "1"; "2"; "3"} ...
%!         [-3; -2; -1; 1; 2; 3]});

## test for NaN and empty strings
%!assert (nthargout (1:3, @grp2idx, [2 2 3 NaN 2 3]),
%!        {[1; 1; 2; NaN; 1; 2] {"2"; "3"} [2; 3]})
%!assert (nthargout (1:3, @grp2idx, {"et" "sa" "sa" "" "et"}),
%!        {[1; 2; 2; NaN; 1] {"et"; "sa"} {"et"; "sa"}})

## FIXME this fails because unique() in core does not yet have set_order
##        option implemented.  See code for code to uncomment once it is
##        implemented in core.
## Test that order when handling strings is by order of appearance
%!xtest <51928> assert (nthargout (1:3, @grp2idx, ["sci"; "oct"; "sci"; "oct"; "oct"]),
%!        {[1; 2; 1; 2; 2] {"sci"; "oct"} ["sci"; "oct"]});
%!xtest <51928> assert (nthargout (1:3, @grp2idx, {"sci"; "oct"; "sci"; "oct"; "oct"}),
%!        {[1; 2; 1; 2; 2] {"sci"; "oct"} {"sci"; "oct"}});
%!xtest <51928> assert (nthargout (1:3, @grp2idx, {"sa" "et" "et" "" "sa"}),
%!        {[1; 2; 2; NaN; 1] {"sa"; "et"} {"sa"; "et"}})
