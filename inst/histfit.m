## Copyright (C) 2003 Alberto Terruzzi <t-albert@libero.it>
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
## @deftypefn {Function File} histfit (@var{data}, @var{nbins})
##
## Plot histogram with superimposed fitted normal density.
##
## @code{histfit (@var{data}, @var{nbins})} plots a histogram of the values in
## the vector @var{data} using @var{nbins} bars in the histogram. With one input
## argument, @var{nbins} is set to the square root of the number of elements in
## data.
##
## Example
##
## @example
## histfit (randn (100, 1))
## @end example
##
## @seealso{bar, hist, pareto}
## @end deftypefn

function histfit (data, nbins)

  if (nargin < 1 || nargin > 2)
    print_usage;
  endif

  if (! isnumeric (x) || ! isreal (x) || ! isvector (x) || isscalar (x))
    error ("histfit: data must be a numeric vector of real numbers.");
  endif

  row = sum (! isnan (data));

  if (nargin < 2)
    nbins = ceil (sqrt (row));
  endif

  [n, xbin] = hist (data, nbins);
  if (any (abs (diff (xbin, 2)) > 10 * max (abs (xbin)) * eps))
    error ("histfit: bins must have uniform width.");
  endif

  ## Compute mu and sigma parameters
  mr = mean (data, "omitnan");
  sr = std (data);
  ## Evenly spaced samples of the expected data range
  x = (-3*sr+mr:0.1*sr:3*sr+mr)';
  [xb, yb] = bar (xbin, n);
  y = normpdf (x, mr, sr);
  binwidth = xbin(2) - xbin(1);
  ## Necessary normalization to overplot the histogram
  y = row * y * binwidth;
  ## Plot density line over histogram.
  plot (xb, yb, ";;b", x, y, ";;r-");

endfunction

%!demo
%! histfit (randn (100, 1))

## testing
%!shared visibility_setting
%! visibility_setting = get (0, "DefaultFigureVisible");
%! set (0, "DefaultFigureVisible", "off");
%!test
%! x = [2, 4, 3, 2, 4, 3, 2, 5, 6, 4, 7, 5, 9, 8, 10, 4, 11];
%! histfit (x);
%! x = [2, 4, 3, 2, NaN, 3, 2, 5, 6, 4, 7, 5, 9, 8, 10, 4, 11];
%! histfit (x);
%! histfit (x, 2);
%!error histfit ();
%!error histfit ([x',x']);
%! set (0, "DefaultFigureVisible", visibility_setting);
