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
## @deftypefn {Private Function} {@var{txt} =} treeLevelText (@var{name}, @var{levels})
##
## The condition a categorical cut tests, as a tree's @code{view} prints it.
##
## @var{txt} reads @qcode{"x1=4"} for a single level and
## @qcode{"x1 in @{1 4 6@}"} for several, as MATLAB R2024a prints them.
##
## @end deftypefn

function txt = treeLevelText (name, levels)

  if (isscalar (levels))
    txt = sprintf ("%s=%g", name, levels);
  else
    txt = sprintf ("%s in {%s}", name, strtrim (sprintf ("%g ", levels)));
  endif

endfunction
