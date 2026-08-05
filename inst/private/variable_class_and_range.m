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
## @deftypefn {Private Function} {[@var{cname}, @var{range}] =} variable_class_and_range (@var{col})
##
## Report the class of a model variable and the range of values it takes.
##
## @var{col} is one variable of the data as it was given.
##
## @var{cname} is its class.  @var{range} is the two-element vector
## @code{[min, max]} for a numeric or logical variable, and otherwise the list
## of levels it takes, as a row and in the variable's own type: a
## @code{categorical} keeps its category order and comes back a
## @code{categorical}, while a @code{string} or a cell array of character
## vectors lists its levels in the order the data presents them and comes back
## as a @code{string} or a cell array respectively.
##
## The level order matters: it is the order the design matrix codes them in, so
## the reference level is the first entry.
##
## This helper is shared by the @code{LinearModel} and
## @code{GeneralizedLinearModel} classes.
##
## @end deftypefn

function [cname, range] = variable_class_and_range (col)

  cname = class (col);
  if (isnumeric (col) || islogical (col))
    fv = double (col(:));
    fv = fv(isfinite (fv));
    if (isempty (fv))
      range = [NaN, NaN];
    else
      range = [min(fv), max(fv)];
    endif
    ## The range keeps the variable's own type, as the level list does.
    if (islogical (col))
      range = logical (range);
    endif
  elseif (isa (col, 'categorical'))
    lv    = categories (col);
    range = categorical (lv(:)', lv);
  elseif (isa (col, 'string'))
    range = unique (col(:)', 'stable');
  elseif (iscell (col))
    lv    = unique (col, 'stable');
    range = lv(:)';
  else
    range = {};
  endif

endfunction
