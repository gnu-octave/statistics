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
## @deftypefn  {statistics} {[@var{q}, @var{s}] =} qqplot (@var{x})
## @deftypefnx {statistics} {[@var{q}, @var{s}] =} qqplot (@var{x}, @var{y})
## @deftypefnx {statistics} {[@var{q}, @var{s}] =} qqplot (@var{x}, @var{dist})
## @deftypefnx {statistics} {[@var{q}, @var{s}] =} qqplot (@var{x}, @var{y}, @var{params})
## @deftypefnx {statistics} {} qqplot (@dots{})
##
## Perform a QQ-plot (quantile plot).
##
## If F is the CDF of the distribution @var{dist} with parameters
## @var{params} and G its inverse, and @var{x} a sample vector of length
## @var{n}, the QQ-plot graphs ordinate @var{s}(@var{i}) = @var{i}-th
## largest element of x versus abscissa @var{q}(@var{i}f) = G((@var{i} -
## 0.5)/@var{n}).
##
## If the sample comes from F, except for a transformation of location
## and scale, the pairs will approximately follow a straight line.
##
## If the second argument is a vector @var{y} the empirical CDF of @var{y}
## is used as @var{dist}.
##
## The default for @var{dist} is the standard normal distribution.  The
## optional argument @var{params} contains a list of parameters of
## @var{dist}.  For example, for a quantile plot of the uniform
## distribution on [2,4] and @var{x}, use
##
## @example
## qqplot (x, "unif", 2, 4)
## @end example
##
## @noindent
## @var{dist} can be any string for which a function @var{distinv} or
## @var{dist_inv} exists that calculates the inverse CDF of distribution
## @var{dist}.
##
## If no output arguments are given, the data are plotted directly.
## @seealso{ppplot}
## @end deftypefn

function [qout, sout] = qqplot (x, dist, varargin)

  if (nargin < 1)
    print_usage ();
  endif

  if (! isnumeric (x) || ! isreal (x) || ! isvector (x) || isscalar (x))
    error ("qqplot: X must be a numeric vector of real numbers");
  endif

  if (nargin == 1)
    f = @probit;
  else
    if (isnumeric (dist))
      f = @(y) empirical_inv (y, dist);
    elseif (ischar (dist) && (exist (invname = [dist "inv"])
                              || exist (invname = [dist "_inv"])))
      f = str2func (invname);
    else
      error ("qqplot: no inverse CDF found for distribution DIST");
    endif
  endif;

  s = sort (x);
  n = length (x);
  t = ((1 : n)' - .5) / n;
  if (nargin <= 2)
    q = f (t);
    q_label = func2str (f);
  else
    q = f (t, varargin{:});
    if (nargin == 3)
      q_label = sprintf ("%s with parameter %g", func2str (f), varargin{1});
    else
      q_label = sprintf ("%s with parameters %g", func2str (f), varargin{1});
      param_str = sprintf (", %g", varargin{2:end});
      q_label = [q_label param_str];
    endif
  endif

  if (nargout == 0)
    plot (q, s, "-x");
    q_label = strrep (q_label, '_inv', '\_inv');
    if (q_label(1) == '@')
      q_label = q_label(6:end);  # Strip "@(y) " from anon. function
    endif
    xlabel (q_label);
    ylabel ("sample points");
  else
    qout = q;
    sout = s;
  endif

endfunction


## Test plotting
%!test
%! hf = figure ("visible", "off");
%! unwind_protect
%!   qqplot ([2 3 3 4 4 5 6 5 6 7 8 9 8 7 8 9 0 8 7 6 5 4 6 13 8 15 9 9]);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

## Test input validation
%!error qqplot ()
%!error <qqplot: X must be a numeric vector> qqplot ({1})
%!error <qqplot: X must be a numeric vector> qqplot (ones (2,2))
%!error <qqplot: X must be a numeric vector> qqplot (1, "foobar")
%!error <qqplot: no inverse CDF found> qqplot ([1 2 3], "foobar")
