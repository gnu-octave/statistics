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
## FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
## more details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn {Private Function} {@var{Yx} =} namesAsLabels (@var{Y}, @var{names})
##
## The response @var{Y} with the class names @var{names} appended as labels in
## its own type.
##
## A class removed for its prior of 0 may be absent from the response while
## the model still names it; appending it lets the classes be ordered and
## counted as if it were present.  Text names become numbers for a numeric
## response, as @code{namedClasses} matches them, and @qcode{'true'} or
## @qcode{'false'} for a logical one.
## @end deftypefn

function Yx = namesAsLabels (Y, names)

  textN = ! (isnumeric (names) || islogical (names));
  if (textN)
    names = cellstr (names);
  endif
  names = names(:);

  if (iscellstr (Y))
    Yx = [Y(:); names];
  elseif (ischar (Y))
    Yx = char ([cellstr(Y); names]);
  elseif (isa (Y, 'string'))
    Yx = [Y(:); string(names)];
  elseif (isa (Y, 'categorical'))
    Yx = [Y(:); categorical(names)];
  elseif (islogical (Y))
    if (textN)
      names = strcmp (names, 'true');
    endif
    Yx = [Y(:); logical(names)];
  else
    if (textN)
      names = str2double (names);
    endif
    Yx = [Y(:); cast(names, class (Y))];
  endif

endfunction
