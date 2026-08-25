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
## @deftypefn {Private Function} {@var{out} =} restore_tables (@var{in})
##
## Rebuild any @code{table} that @code{load} handed back as a plain structure.
##
## Saving to a binary file keeps a table's contents but not its class, so it
## comes back as a structure carrying @qcode{VariableNames} and
## @qcode{VariableValues} among the table's other properties.  This helper
## turns such a structure back into the table it was, and recurses through the
## fields of any other structure and the elements of any cell, so that a table
## nested inside one is restored with it.
##
## Anything that is not a saved table is returned unchanged.
##
## This helper is shared by the @code{load_model} methods of the learner
## classes, whose saved properties may hold tables.
##
## @end deftypefn

function out = restore_tables (in)

  if (isstruct (in) && isscalar (in) && is_saved_table (in))
    out = table (in.VariableValues{:}, 'VariableNames', in.VariableNames);
    if (! isempty (in.RowNames))
      out.Properties.RowNames = in.RowNames;
    endif
    return;
  endif

  if (isstruct (in) && isscalar (in))
    out = in;
    for [val, key] = in
      out.(key) = restore_tables (val);
    endfor
    return;
  endif

  if (iscell (in))
    out = cellfun (@restore_tables, in, 'UniformOutput', false);
    return;
  endif

  out = in;

endfunction

## A structure is a saved table when it carries the fields the table class
## writes out; testing for the pair that holds the data is enough, since no
## structure the models save has both.
function tf = is_saved_table (s)
  tf = isfield (s, 'VariableNames') && isfield (s, 'VariableValues') ...
       && isfield (s, 'RowNames') && iscell (s.VariableValues);
endfunction
