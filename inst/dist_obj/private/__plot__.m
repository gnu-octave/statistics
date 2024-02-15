## Copyright (C) 2024 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn {Private Function} [@var{h}] = __plot__ (@var{pd}, @var{Discrete}, @var{varargin})
##
## Plot a probability distribution object.
##
## @end deftypefn

function h = __plot__ (pd, Discrete, varargin)

  ## Add defaults (Discrete is passed by the calling method)
  h = [];
  PlotType = "pdf";

  ## Parse optional agruments
  if (mod (numel (varargin), 2) != 0)
    error ("plot: optional arguments must be in NAME-VALUE pairs.");
  endif
  while (numel (varargin) > 0)
    switch (tolower (varargin{1}))

      case "plottype"
        ValidTypes = {"pdf", "cdf", "probability"};
        try
          selected_T = strcmpi (varargin{2}, ValidTypes);
        catch
          error ("plot: invalid VALUE size for 'Parameter' argument.");
        end_try_catch
        if (! any (selected_T) || sum (selected_T) > 1)
          error ("plot: invalid VALUE for 'PlotType' argument.");
        endif
        PlotType = ValidTypes{selected_T};

      case "discrete"
        if (! (islogical (varargin{2}) && isscalar (varargin{2})))
          error ("plot: invalid VALUE for 'Discrete' argument.");
        endif
        ## Only for discrete distributions this can be changed by the user
        if (Discrete)
          Discrete = varargin{2};
        endif

      case "parent"
        if (! isaxes (varargin{2}))
          error ("plot: invalid VALUE for 'Parent' argument.");
        endif
        h = varargin{2};

      otherwise
        error ("plot: invalid NAME for optional argument.");

    endswitch
    varargin([1:2]) = [];
  endwhile

  ## Get current axes or create new ones
  if (isempty (h))
    h = cga;
  endif



endfunction
