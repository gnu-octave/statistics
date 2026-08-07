## Copyright (C) 2026 Andreas Bertsatos <abertsatos@biol.uoa.gr>
##
## This file is part of the statistics package for GNU Octave.
##
## This program is free software: you can redistribute it and/or
## modify it under the terms of the GNU General Public License as
## published by the Free Software Foundation, either version 3 of the
## License, or (at your option) any later version.
##
## This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; see the file COPYING.  If not, see
## <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn {Private Function} {[@var{inds}, @var{ngroups}, @var{gsize}, @var{missidx}, @var{grpvars}] =} __resolve_groups__ (@var{grpvars})
##
## Reduce grouping variables to a group index per observation.
##
## @var{inds} numbers each observation by the group it belongs to, @var{ngroups}
## counts the groups and @var{gsize} their sizes.  @var{missidx} flags the
## observations dropped for holding a missing value, and @var{grpvars} is
## returned with those observations removed.
##
## One grouping variable is given as a vector, several as a matrix or cell
## array whose columns are the variables; a row then belongs to a group only if
## it matches on every one of them.
##
## @end deftypefn

function [inds, ngroups, gsize, missidx, grpvars] = __resolve_groups__ (grpvars)

  if (isvector (grpvars))
    missidx = ismissing (grpvars);
    if (any (missidx))
      grpvars(missidx) = [];
    endif
    if (isa (grpvars, 'categorical'))
      [~, idx, inds] = unique (grpvars, 'stable');
    else
      [~, idx, inds] = __unique__ (grpvars, 'stable');
    endif
  elseif (ismatrix (grpvars))
    missidx = any (ismissing (grpvars), 2);
    if (any (missidx))
      grpvars(missidx, :) = [];
    endif
    if (isa (grpvars, 'categorical'))
      [~, idx, inds] = unique (grpvars, 'rows', 'stable');
    else
      [~, idx, inds] = __unique__ (grpvars, 'rows', 'stable');
    endif
  else
    error (strcat ("cvpartition: invalid value for optional", ...
                   " paired argument 'GroupingVariables'."));
  endif

  ngroups = numel (idx);
  gsize = zeros (1, ngroups);
  for i = 1:ngroups
    gsize(i) = sum (inds == i);
  endfor

endfunction
