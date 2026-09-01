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
## @deftypefn {Private Function} {@var{stop} =} objective_history_fcn (@var{x}, @var{optimValues}, @var{state})
##
## Append the current objective to the history, as @code{fminunc}'s
## @code{OutputFcn}.
##
## @var{stop} is always false, so a fit is never halted from here, and @var{x}
## is unused.  The @qcode{'init'} call carries the objective at the starting
## weights, which is the first entry of the history MATLAB reports.
##
## @end deftypefn

function stop = objective_history_fcn (x, optimValues, state)

  stop = false;
  if (any (strcmp (state, {"init", "iter"})))
    objective_history ("add", optimValues.fval);
  endif

endfunction
