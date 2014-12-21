## Copyright (C) 2014 Nir Krakauer
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
## @deftypefn{Function File}@var{s} = set (@var{c}, @var{varargin})
## Set field(s) in a @samp{cvpartition} object.
##
## @seealso{cvpartition}
## @end deftypefn

function s = set (c, varargin)
  s = struct(c);
  if (length (varargin) < 2 || rem (length (varargin), 2) != 0)
    error ("set: expecting property/value pairs");
  endif
  while (length (varargin) > 1)
    prop = varargin{1};
    val = varargin{2};
    varargin(1:2) = [];
    if (ischar (prop))
      switch (prop)
        case {"classes", "inds", "n_classes", "NumObservations", "NumTestSets", "TestSize", "TrainSize", "Type"}
          s = setfield (s, prop, val);
        otherwise
          error ("set: invalid property %s", f);
      endswitch
    else
      error ("set: expecting the property to be a string");
    endif
  endwhile
  s = class (s, "cvpartition");
endfunction

