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
## @deftypefn {Private Function} {[@var{S}, @var{E}, @var{errmsg}] =} bagMds (@var{P}, @var{g}, @var{colors}, @var{coords})
##
## Classical multidimensional scaling of a proximity matrix, with its plot.
##
## @var{P} is an @math{NxN} proximity matrix, scaled as the distances
## @code{1 - @var{P}}.  @var{S} holds the scaled coordinates, one column per
## positive eigenvalue, and @var{E} the eigenvalues, as @code{cmdscale}
## returns them.
##
## @var{coords} holds two or three column indices of @var{S}, which must exist
## whether or not anything is drawn.  @var{colors}, when not empty, is a
## character vector with one color letter per class, and those coordinates are
## drawn as overlaid scatter plots, one per class.  @var{g} holds
## the class index of each observation, or is empty to draw every observation
## in the first color; a class beyond the number of letters is not drawn.
## @var{errmsg} is the body of the message the caller should raise, or empty.
##
## @end deftypefn

function [S, E, errmsg] = bagMds (P, g, colors, coords)

  S = [];
  E = [];
  errmsg = "";
  if (isa (colors, 'string') && isscalar (colors))
    colors = char (colors);
  endif
  if (! (isempty (colors) || (ischar (colors) && isrow (colors))))
    errmsg = "'Colors' must be a character vector or a string scalar.";
    return;
  endif
  if (! (isnumeric (coords) && isvector (coords) && isreal (coords)
         && any (numel (coords) == [2, 3]) && all (coords >= 1)
         && all (coords == fix (coords))))
    errmsg = strcat ("'MDSCoordinates' must be a vector of two or three", ...
                     " positive integers.");
    return;
  endif

  [S, E] = cmdscale (1 - P);
  if (any (coords > columns (S)))
    S = [];
    E = [];
    errmsg = strcat ("'MDSCoordinates' must not exceed the number of", ...
                     " scaled coordinates.");
    return;
  endif
  if (isempty (colors))
    return;
  endif

  if (isempty (g))
    g = ones (rows (S), 1);
  endif
  held = ishold ();
  for k = 1:min (numel (colors), max ([0; g(! isnan (g))]))
    idx = g == k;
    if (numel (coords) == 2)
      plot (S(idx,coords(1)), S(idx,coords(2)), [colors(k), '.']);
    else
      plot3 (S(idx,coords(1)), S(idx,coords(2)), S(idx,coords(3)), ...
             [colors(k), '.']);
    endif
    hold on;
  endfor
  if (! held)
    hold off;
  endif

endfunction
