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
## @deftypefn {Private Function} {@var{label} =} missingLabels (@var{Y}, @var{n})
##
## A column of @var{n} missing labels in the type of the response @var{Y}.
##
## A cell array of character vectors gives empty names, a character matrix
## rows of blanks as wide as @var{Y}, a numeric response @code{NaN}, a
## categorical response @qcode{<undefined>} keeping the categories of @var{Y},
## and a string response @qcode{<missing>}.  A logical response has no missing
## value and gives @code{false}.
##
## @end deftypefn

function label = missingLabels (Y, n)

  if (iscellstr (Y))
    label = repmat ({''}, n, 1);
  elseif (islogical (Y))
    label = false (n, 1);
  elseif (isnumeric (Y))
    label = nan (n, 1);
  elseif (ischar (Y))
    label = repmat (' ', n, columns (Y));
  elseif (isa (Y, 'categorical'))
    label = categorical (repmat ({''}, n, 1), categories (Y), ...
                         'Ordinal', isordinal (Y), ...
                         'Protected', isprotected (Y));
  elseif (isa (Y, 'string'))
    label = string (nan (n, 1));
  endif

endfunction
