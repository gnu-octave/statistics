## Copyright (C) 2014 Nir Krakauer
##
## This file is part of the statistics package for GNU Octave.
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; If not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn  {statistics} {@var{Cnew} =} set (@var{C}, @var{field}, @var{value})
##
## Set @var{field}, in a @samp{cvpartition} object, @var{C}.
##
## @seealso{@cvpartition/cvpartition}
## @end deftypefn

function s = set (c, varargin)
  s = struct(c);
  if (length (varargin) < 2 || rem (length (varargin), 2) != 0)
    error ("set: expecting field/value pairs.");
  endif
  while (length (varargin) > 1)
    prop = varargin{1};
    val = varargin{2};
    varargin(1:2) = [];
    if (ischar (prop))
      switch (prop)
        case {"classes", "inds", "n_classes", "NumObservations", ...
              "NumTestSets", "TestSize", "TrainSize", "Type"}
          s = setfield (s, prop, val);
        otherwise
          error ("set: invalid field %s.", prop);
      endswitch
    else
      error ("set: expecting the field to be a string.");
    endif
  endwhile
  s = class (s, "cvpartition");
endfunction

%!shared C
%! C = cvpartition (ones (10, 1), "KFold", 5);
%!test
%! Cnew = set (C, "inds", [1 2 2 2 3 4 3 4 5 5]');
%! assert (get (Cnew, "inds"), [1 2 2 2 3 4 3 4 5 5]');
%!error<set: expecting field/value pairs.> set (C)
%!error<set: expecting field/value pairs.> set (C, "NumObservations")
%!error<set: invalid field some.> set (C, "some", 15)
%!error<set: expecting the field to be a string.> set (C, 15, 15)
