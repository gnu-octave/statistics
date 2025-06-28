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
## @deftypefn {Private Function} __disp__ (@var{pd}, @var{distname})
##
## Display summary of a probability distribution object.
##
## @end deftypefn

function ci = __disp__ (pd, distname)

  if (isscalar (pd))

    ## Handle special case of PiecewiseLinearDistribution
    if (strcmpi (pd.DistributionName, "PiecewiseLinearDistribution"))
      ## Print distribution header
      fprintf ("  %s\n\n", pd.DistributionName);
      ## Print parameter values
      for i = 1:numel (pd.x)
        fprintf ("F(%g) = %g\n", pd.x(i), pd.Fx(i));
      endfor
      ## Print truncation interval if applicable
      if (pd.IsTruncated)
        fprintf ("  Truncated to the interval [%g, %g]\n\n", pd.Truncation);
      else
        fprintf ("\n");
      endif

    ## Handle special case of MultinomialDistribution
    elseif (strcmpi (pd.DistributionName, "MultinomialDistribution"))
      ## Print distribution header
      fprintf ("  %s\n\n", pd.DistributionName);
      ## Print parameter values
      fprintf ("  Probabilities:\n");
      fprintf ("    %0.4f", pd.Probabilities);
      fprintf ("\n\n");
      ## Print truncation interval if applicable
      if (pd.IsTruncated)
        fprintf ("  Truncated to the interval [%g, %g]\n\n", pd.Truncation);
      endif

    ## Handle all other cases
    else
      ## Get required length for parameter values
      PVlen = max (arrayfun (@(x) numel (sprintf ("%g", x)), ...
                             pd.ParameterValues));
      PVstr = sprintf ("%%%dg", PVlen);

      ## Prepare template for fitted and not fitted distributions
      pat1 = ["  %+7s = ", PVstr, "   [%g, %g]\n"];
      pat2 = ["  %+7s = ", PVstr, "\n"];

      ## Grad distributions that are non fittable
      if (any (strcmpi (pd.DistributionCode, {"unif", "tri", "logu"})))
        fitted = false;
        ParameterIsFixed = true;
      elseif (all (pd.ParameterIsFixed))
        fitted = false;
        ParameterIsFixed = pd.ParameterIsFixed;
      else
        fitted = true;
        ParameterIsFixed = pd.ParameterIsFixed;
      endif

      ## Print distribution header
      fprintf ("  %s\n\n", pd.DistributionName);
      fprintf ("  %s\n", distname);
      ## Print parameter values
      for i = 1:pd.NumParameters
        if (fitted && ! ParameterIsFixed(i))
          fprintf (pat1, pd.ParameterNames{i}, pd.ParameterValues(i), ...
                         pd.ParameterCI(1,i), pd.ParameterCI(2,i));
        else
          fprintf (pat2, pd.ParameterNames{i}, pd.ParameterValues(i));
        endif
      endfor
      ## Print truncation interval if applicable
      if (pd.IsTruncated)
        fprintf ("  Truncated to the interval [%g, %g]\n\n", pd.Truncation);
      else
        fprintf ("\n");
      endif
    endif
  else
    fprintf ("%dx%d %s array\n", size (pd), class (pd));
  endif

endfunction
