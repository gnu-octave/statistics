## Copyright (C) 2022 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn {Function File} @var{hCDF} = cdfplot (@var{x})
## @deftypefnx {Function File} [@var{hCDF}, @var{stats}] = cdfplot (@var{x})
##
## Display an empirical cumulative distribution function.
##
## @code{@var{hCDF} = cdfplot (@var{x})} plots an empirical cumulative
## distribution function (CDF) of the observations in the data sample vector
## @var{x}.  @var{x} may be a row or column vector, and represents a random
## sample of observations from some underlying distribution.
##
## @code{cdfplot} plots F(x), the empirical (or sample) CDF versus the
## observations in @var{x}. The empirical CDF, F(x), is defined as follows:
##
## F(x) = (Number of observations <= x) / (Total number of observations)
##
## for all values in the sample vector @var{x}.  NaNs are ignored.  @var{hCDF}
## is the handle of the empirical CDF curve (a handle hraphics 'line' object).
##
## @code{[@var{hCDF}, @var{stats}] = cdfplot (@var{x})} also returns a structure
## with the following fields as a statistical summary.
##
## @multitable @columnfractions 0.05 0.3 0.65
## @item @tab STATS.min @tab minimum value of @var{x}
## @item @tab STATS.max @tab maximum value of @var{x}
## @item @tab STATS.mean @tab sample mean of @var{x}
## @item @tab STATS.median @tab sample median (50th percentile) of @var{x}
## @item @tab STATS.std @tab sample standard deviation of @var{x}
## @end multitable
##
## @seealso{qqplot, cdfcalc}
## @end deftypefn

function [hCDF, stats] = cdfplot (x)
  ## Check number of input arguments
  narginchk (1,1);
  ## Calculate sample cdf
  [yy, xx, ~, ~, eid] = cdfcalc (x);
  ## Check for errors returned from cdfcalc
  if (strcmpi (eid, "VectorRequired"))
    error ("cdfplot: vector required as input.");
  elseif (strcmpi (eid, "NotEnoughData"))
    error("cdfplot: not enough data.");
  endif
  ## Create vectors for plotting
  k = length (xx);
  n = reshape (repmat (1:k, 2, 1), 2*k, 1);
  xCDF = [-Inf; xx(n); Inf];
  yCDF = [0; 0; yy(1+n)];
  ## Plot cdf
  h = plot (xCDF, yCDF);
  grid  ('on')
  xlabel ("x")
  ylabel ("F(x)")
  title ("CDF plot of x");

  ## Return requested output arguments
  if (nargout > 0)
    hCDF = h;
  endif
  if (nargout > 1)
    stats.min = nanmin(x);
    stats.max = nanmax(x);
    stats.mean = mean(x, "omitnan");
    stats.median = nanmedian(x);
    stats.std = nanstd(x);
  endif
endfunction

%!demo
%! x = randn(100,1);
%! cdfplot (x);

## Get current figure visibility so it can be restored after tests
%!shared visibility_setting
%! visibility_setting = get (0, "DefaultFigureVisible");
%! set (0, "DefaultFigureVisible", "off");
%!test
%! x = [2, 4, 3, 2, 4, 3, 2, 5, 6, 4];
%! [hCDF, stats] = cdfplot (x);
%! assert (stats.min, 2);
%! assert (stats.max, 6);
%! assert (stats.median, 3.5);
%! assert (stats.std, 1.35400640077266, 1e-14);
%!error cdfplot ();
%!error cdfplot ([x',x']);
%!error cdfplot ([NaN, NaN, NaN, NaN]);
%!test
%! x = randn(100,1);
%! cdfplot (x);
%! set (0, "DefaultFigureVisible", visibility_setting);
