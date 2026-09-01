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
## @deftypefn {Private Function} {@var{out} =} objective_history (@var{action}, @var{val})
##
## Collect the objective values @code{fminunc} walks through.
##
## @var{action} is @qcode{'reset'} to begin a fit, @qcode{'add'} to append
## @var{val}, or @qcode{'get'} to return everything gathered since the reset.
## Only @qcode{'get'} returns anything, the others leaving @var{out} empty.
##
## @code{fminunc} offers no history of its own and its @code{OutputFcn} cannot
## return one, so the values are gathered here and reported as
## @code{FitInfo.Objective}, which would otherwise hold the final value alone.
## One fit runs at a time, which is what lets a single buffer serve both
## classes.
##
## @end deftypefn

function out = objective_history (action, val)

  persistent buf

  out = [];
  switch (action)
    case 'reset'
      buf = zeros (0, 1);
    case 'add'
      buf(end+1, 1) = val;
    case 'get'
      out = buf;
    otherwise
      error ("objective_history: unknown action '%s'.", action);
  endswitch

endfunction
