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
## @deftypefn  {statistics} {[@var{n}, @var{c}] =} ecdfhist (@var{f}, @var{x})
## @deftypefnx {statistics} {[@var{n}, @var{c}] =} ecdfhist (@var{f}, @var{x}, @var{m})
## @deftypefnx {statistics} {[@var{n}, @var{c}] =} ecdfhist (@var{f}, @var{x}, @var{centers})
## @deftypefnx {statistics} {} ecdfhist (@dots{})
## @deftypefnx {statistics} {} ecdfhist (@var{ax}, @dots{})
##
## Create a histogram from the output of @code{ecdf}.
##
## @code{[@var{n}, @var{c}] = ecdfhist (@var{f}, @var{x})} takes the empirical
## cumulative distribution function @var{f} evaluated at the points @var{x}, as
## computed by @code{ecdf}, and returns the heights @var{n} of histogram bars
## for 10 equally spaced bins together with their centers @var{c}.  Unlike a
## count histogram, the bar heights are normalized so that the area of the
## histogram is equal to 1, giving an estimate of the probability density
## function.
##
## @code{[@var{n}, @var{c}] = ecdfhist (@var{f}, @var{x}, @var{m})} uses @var{m}
## equally spaced bins.
##
## @code{[@var{n}, @var{c}] = ecdfhist (@var{f}, @var{x}, @var{centers})} uses
## bins with the specified centers, given as a vector of monotonically
## increasing values.
##
## @code{ecdfhist (@dots{})} without output arguments plots the histogram.
##
## @code{ecdfhist (@var{ax}, @dots{})} plots into the axes @var{ax} instead of
## the current axes.
##
## The probability mass assigned to each bin is the sum of the increments of the
## empirical cdf, @code{diff (@var{f})}, over the points @var{x} that fall
## closest to the corresponding bin center; ties are assigned to the lower
## center.  Each bar height is that mass divided by the bin width.
##
## @seealso{ecdf, cdfplot, hist, histogram}
## @end deftypefn

function [nout, cout] = ecdfhist (varargin)

  ## Detect a leading axes handle
  ax = [];
  if (numel (varargin) > 0 && isaxes (varargin{1}))
    ax = varargin{1};
    varargin(1) = [];
  endif

  if (numel (varargin) < 2 || numel (varargin) > 3)
    print_usage ();
  endif

  f = varargin{1};
  x = varargin{2};
  if (! isnumeric (f) || ! isreal (f) || ! isnumeric (x) || ! isreal (x))
    error ("ecdfhist: F and X must be real numeric vectors.");
  endif
  f = f(:);
  x = x(:);
  if (numel (f) != numel (x))
    error ("ecdfhist: F and X must have the same length.");
  endif
  if (numel (f) < 2)
    error ("ecdfhist: F and X must have at least two elements.");
  endif

  ## Probability masses are the cdf increments, located at x(2:end)
  xm = x(2:end);
  masses = diff (f);

  ## Determine the bin centers and widths
  if (numel (varargin) < 3 || isempty (varargin{3}))
    m = 10;
    centers = [];
  elseif (isscalar (varargin{3}))
    m = varargin{3};
    if (! isnumeric (m) || ! isreal (m) || m < 1 || fix (m) != m)
      error ("ecdfhist: M must be a positive integer.");
    endif
    centers = [];
  else
    centers = varargin{3};
    if (! isnumeric (centers) || ! isreal (centers) || ! isvector (centers))
      error ("ecdfhist: CENTERS must be a real numeric vector.");
    endif
    centers = centers(:).';
  endif

  if (isempty (centers))
    xlo = min (xm);
    xhi = max (xm);
    if (xlo == xhi)
      xlo -= 0.5;
      xhi += 0.5;
    endif
    edges = linspace (xlo, xhi, m + 1);
    centers = edges(1:end-1) + diff (edges) / 2;
    binwidth = diff (edges);
  elseif (numel (centers) == 1)
    binwidth = 1;
  else
    ec = (centers(1:end-1) + centers(2:end)) / 2;
    binwidth = diff ([2 * centers(1) - ec(1), ec, 2 * centers(end) - ec(end)]);
  endif
  ncnt = numel (centers);

  ## Assign each mass to the nearest center (ties go to the lower center)
  [~, bin] = min (abs (xm - centers), [], 2);
  binmass = accumarray (bin, masses, [ncnt, 1]);
  heights = binmass.' ./ binwidth;

  if (nargout == 0)
    if (isempty (ax))
      ax = newplot ();
    endif
    bar (ax, centers, heights, 'hist');
  else
    nout = heights;
    cout = centers;
  endif

endfunction

%!demo
%! ## Histogram (density estimate) from the empirical cdf of a random sample.
%!
%! x = randn (100, 1);
%! [f, xx] = ecdf (x);
%! ecdfhist (f, xx);
%! title ("ecdfhist of a standard normal sample");

%!demo
%! ## Compare the empirical density with a finer set of bins.
%!
%! x = exprnd (2, 200, 1);
%! [f, xx] = ecdf (x);
%! ecdfhist (f, xx, 20);
%! title ("ecdfhist of an exponential sample (20 bins)");

## Test output
%!test
%! f = [0 0.1 0.3 0.6 0.8 0.9 1.0]';
%! x = [1 1 2 3 4 5 8]';
%! [n, c] = ecdfhist (f, x);
%! assert (c, 1.35:0.7:7.65, 1e-12);
%! assert (n, [0.1 0.2 0.3 0 0.2 0.1 0 0 0 0.1] / 0.7, 1e-12);
%! assert (sum (n .* 0.7), 1, 1e-12);
%!test
%! f = [0 0.1 0.3 0.6 0.8 0.9 1.0]';
%! x = [1 1 2 3 4 5 8]';
%! [n, c] = ecdfhist (f, x, 5);
%! assert (c, [1.7 3.1 4.5 5.9 7.3], 1e-12);
%! assert (n, [0.3 0.3 0.3 0 0.1] / 1.4, 1e-12);
%!test
%! f = [0 0.1 0.3 0.6 0.8 0.9 1.0]';
%! x = [1 1 2 3 4 5 8]';
%! [n, c] = ecdfhist (f, x, [1.5 2.5 3.5 4.5 5.5 6.5 7.5]);
%! assert (c, [1.5 2.5 3.5 4.5 5.5 6.5 7.5], 1e-12);
%! assert (n, [0.3 0.3 0.2 0.1 0 0 0.1], 1e-12);
%!test  # returns only the heights when a single output is requested
%! f = [0 0.1 0.3 0.6 0.8 0.9 1.0]';
%! x = [1 1 2 3 4 5 8]';
%! n = ecdfhist (f, x, 5);
%! assert (n, [0.3 0.3 0.3 0 0.1] / 1.4, 1e-12);

## Test plotting
%!test
%! hf = figure ("visible", "off");
%! unwind_protect
%!   f = [0 0.1 0.3 0.6 0.8 0.9 1.0]';
%!   x = [1 1 2 3 4 5 8]';
%!   ecdfhist (f, x);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect
%!test
%! hf = figure ("visible", "off");
%! unwind_protect
%!   ax = axes ("parent", hf);
%!   f = [0 0.1 0.3 0.6 0.8 0.9 1.0]';
%!   x = [1 1 2 3 4 5 8]';
%!   ecdfhist (ax, f, x);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

## Test input validation
%!error <Invalid call to ecdfhist> ecdfhist ()
%!error <Invalid call to ecdfhist> ecdfhist ([0 1])
%!error <ecdfhist: F and X must be real numeric vectors.> ecdfhist ([0 1], {1 2})
%!error <ecdfhist: F and X must have the same length.> ecdfhist ([0 1], [1 2 3])
%!error <ecdfhist: F and X must have at least two elements.> ecdfhist (1, 1)
%!error <ecdfhist: M must be a positive integer.> ecdfhist ([0 1], [1 2], 0)
%!error <ecdfhist: M must be a positive integer.> ecdfhist ([0 1], [1 2], 2.5)
