## Copyright (C) 2003 Alberto Terruzzi <t-albert@libero.it>
## Copyright (C) 2022-2023 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {statistics} {} histfit (@var{x}, @var{nbins})
## @deftypefnx {statistics} @var{h} = histfit (@var{x}, @var{nbins})
##
## Plot histogram with superimposed fitted normal density.
##
## @code{histfit (@var{x}, @var{nbins})} plots a histogram of the values in
## the vector @var{x} using @var{nbins} bars in the histogram.  With one input
## argument, @var{nbins} is set to the square root of the number of elements in
## @var{x}.
##
## @code{@var{h} = histfit (@var{x}, @var{nbins})} returns the bins and fitted
## line handles of the plot in @var{h}.
##
## Example
##
## @example
## histfit (randn (100, 1))
## @end example
##
## @seealso{bar, hist, pareto}
## @end deftypefn

function [varargout] = histfit (x, nbins)

  if (nargin < 1 || nargin > 2)
    print_usage;
  endif

  if (! isnumeric (x) || ! isreal (x) || ! isvector (x) || isscalar (x))
    error ("histfit: X must be a numeric vector of real numbers.");
  endif

  row = sum (! isnan (x));

  if (nargin < 2)
    nbins = ceil (sqrt (row));
  endif

  [n, xbin] = hist (x, nbins);
  if (any (abs (diff (xbin, 2)) > 10 * max (abs (xbin)) * eps))
    error ("histfit: bins must have uniform width.");
  endif

  ## Compute mu and sigma parameters
  mr = mean (x, "omitnan");
  sr = std (x);
  ## Evenly spaced samples of the expected range in X
  x = (-3*sr+mr:0.1*sr:3*sr+mr)';
  [xb, yb] = bar (xbin, n);
  y = normpdf (x, mr, sr);
  binwidth = xbin(2) - xbin(1);
  ## Necessary normalization to overplot the histogram
  y = row * y * binwidth;
  ## Plot density line over histogram.
  h = plot (xb, yb, ";;b", x, y, ";;r-");

  ## Return the plot's handle if requested
  if (nargout == 1)
    varargout{1} = h;
  endif
endfunction

%!demo
%! histfit (randn (100, 1))

## testing
%!shared visibility_setting
%! visibility_setting = get (0, "DefaultFigureVisible");
%!test
%! set (0, "DefaultFigureVisible", "off");
%! x = [2, 4, 3, 2, 4, 3, 2, 5, 6, 4, 7, 5, 9, 8, 10, 4, 11];
%! histfit (x);
%! x = [2, 4, 3, 2, NaN, 3, 2, 5, 6, 4, 7, 5, 9, 8, 10, 4, 11];
%! histfit (x);
%! histfit (x, 3);
%! set (0, "DefaultFigureVisible", visibility_setting);
%!error histfit ();
%!error histfit ([x',x']);
