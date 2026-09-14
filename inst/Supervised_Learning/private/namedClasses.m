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
## @deftypefn {Private Function} {[@var{pos}, @var{errmsg}] =} namedClasses (@var{C}, @var{ClassNames})
## Where each class a @qcode{'ClassNames'} option names sits among the classes
## of a response.
##
## @var{C} holds the distinct classes of the response and @var{pos} the index
## into it of each name in @var{ClassNames}, zero for a name it does not hold.
## Text matches text, and numbers and logicals match by value.  A text name
## matches a numeric class by the text @code{num2str} writes and a logical
## class by @qcode{'true'} or @qcode{'false'}; a numeric or logical name never
## matches a text class.  @var{errmsg} is empty when every name matched, and
## otherwise the body of the message the caller raises under its own name.
## @end deftypefn

function [pos, errmsg] = namedClasses (C, ClassNames)

  errmsg = "";
  textC = ! (isnumeric (C) || islogical (C));
  textN = ! (isnumeric (ClassNames) || islogical (ClassNames));
  if (textC && ! textN)
    pos = zeros (numel (ClassNames), 1);
  elseif (! textC && textN)
    if (islogical (C))
      words = {'false'; 'true'};
      T = words(double (C(:)) + 1);
    else
      T = arrayfun (@num2str, C(:), 'UniformOutput', false);
    endif
    pos = labelIndices (T, ClassNames);
  else
    pos = labelIndices (C, ClassNames);
  endif
  if (isempty (pos) || ! all (pos > 0))
    errmsg = "not all 'ClassNames' are present in Y.";
  endif

endfunction
