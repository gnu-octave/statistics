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
## @deftypefn {statistics} {[@var{f}, @var{dw}] =} parseDistanceWeight (@var{DistanceWeight}, @var{classname})
##
## Parse a nearest-neighbour distance weight.
##
## @var{dw} is the weight as it is reported, a character vector naming a
## built-in weight or the @code{func2str} form of a supplied handle, and
## @var{f} is the callable that applies it.  The property holds the text, as
## MATLAB does, and the callable is kept out of sight beside it.
##
## @end deftypefn

function [f, dw] = parseDistanceWeight (DistanceWeight, classname)

  if (! (ischar (DistanceWeight) || is_function_handle (DistanceWeight)))
    error (strcat ("%s: 'DistanceWeight' must be a character vector or a", ...
                   " function handle."), classname);
  endif

  if (is_function_handle (DistanceWeight))
    m = eye (5);
    if (! isequal (size (m), size (DistanceWeight (m))))
      error (strcat ("%s: function handle for 'DistanceWeight' must return", ...
                     " the same size as its input."), classname);
    endif
    f = DistanceWeight;
    dw = func2str (DistanceWeight);
    return;
  endif

  switch (lower (DistanceWeight))
    case 'equal'
      f = @(d) ones (size (d));
    case 'inverse'
      f = @(d) d .^ (-1);
    case 'squaredinverse'
      f = @(d) d .^ (-2);
    otherwise
      error (strcat ("%s: 'DistanceWeight' must be 'equal', 'inverse',", ...
                     " 'squaredinverse', or a function handle."), classname);
  endswitch
  dw = lower (DistanceWeight);

endfunction
