## Copyright (C) 1995-2017 Kurt Hornik
## Copyright (C) 2022-2023 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {statistics} {} ppplot (@var{x}, @var{dist})
## @deftypefnx {statistics} {} ppplot (@var{x}, @var{dist}, @var{params})
## @deftypefnx {statistics} {[@var{p}, @var{y}] =} ppplot (@var{x}, @var{dist}, @var{params})
##
## Perform a PP-plot (probability plot).
##
## If F is the CDF of the distribution @var{dist} with parameters
## @var{params} and @var{x} a sample vector of length @var{n}, the PP-plot
## graphs ordinate @var{y}(@var{i}) = F (@var{i}-th largest element of
## @var{x}) versus abscissa @var{p}(@var{i}) = (@var{i} - 0.5)/@var{n}.  If
## the sample comes from F, the pairs will approximately follow a straight
## line.
##
## The default for @var{dist} is the standard normal distribution.
##
## The optional argument @var{params} contains a list of parameters of
## @var{dist}.
##
## For example, for a probability plot of the uniform distribution on [2,4]
## and @var{x}, use
##
## @example
## ppplot (x, "unif", 2, 4)
## @end example
##
## @noindent
## @var{dist} can be any string for which a function @var{distcdf} that
## calculates the CDF of distribution @var{dist} exists.
##
## If no output is requested then the data are plotted immediately.
## @seealso{qqplot}
## @end deftypefn

function [p, y] = ppplot (x, dist, varargin)

  if (nargin < 1)
    print_usage ();
  endif

  if (! isnumeric (x) || ! isreal (x) || ! isvector (x) || isscalar (x))
    error ("ppplot: X must be a numeric vector of real numbers");
  endif

  s = sort (x);
  n = length (x);
  p = ((1 : n)' - 0.5) / n;
  if (nargin == 1)
    F = @stdnormal_cdf;
  elseif (! ischar (dist))
    error ("ppplot: DIST must be a string");
  else
    F = str2func ([dist "cdf"]);
  endif

  if (nargin <= 2)
    y = feval (F, s);
  else
    y = feval (F, s, varargin{:});
  endif

  if (nargout == 0)
    plot (p, y);
    axis ([0, 1, 0, 1]);
  endif

endfunction

function p = stdnormal_cdf (x)
  p = 0.5 * erfc (x ./ sqrt (2));
endfunction


## Test plotting
%!test
%! hf = figure ("visible", "off");
%! unwind_protect
%!   ppplot ([2 3 3 4 4 5 6 5 6 7 8 9 8 7 8 9 0 8 7 6 5 4 6 13 8 15 9 9]);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

## Test input validation
%!error ppplot ()
%!error <ppplot: X must be a numeric vector> ppplot (ones (2,2))
%!error <ppplot: X must be a numeric vector> ppplot (1, 2)
%!error <ppplot: DIST must be a strin> ppplot ([1 2 3 4], 2)
