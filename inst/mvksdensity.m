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
## @deftypefn  {statistics} {@var{f} =} mvksdensity (@var{x}, @var{pts}, @var{Name}, @var{Value})
##
## Multivariate kernel smoothing density estimate.
##
## @code{@var{f} = mvksdensity (@var{x}, @var{pts})} computes a probability
## density estimate of the sample in the @math{NxD} matrix @var{x}, evaluated at
## the points in the @math{MxD} matrix @var{pts}.  Each row of @var{x} is a
## single @math{D}-dimensional observation, and each row of @var{pts} is a point
## at which to evaluate the estimate.  The result @var{f} is an @math{Mx1}
## vector, with one density value per row of @var{pts}.
##
## The density estimate uses a product kernel: the multivariate kernel is the
## product of the univariate kernels applied to each dimension, each with its
## own bandwidth.
##
## The following @qcode{Name-Value} pairs are supported:
##
## @multitable @columnfractions 0.18 0.82
## @headitem @var{Name} @tab @var{Value}
##
## @item @qcode{'Bandwidth'} @tab The kernel bandwidth, either a positive scalar
## applied to every dimension or a @math{1xD} vector of positive values, one per
## dimension.  The default is a diagonal normal-reference (Silverman) rule
## computed from @var{x}.
##
## @item @qcode{'Kernel'} @tab The smoothing kernel applied in each dimension:
## @qcode{'normal'} (default), @qcode{'box'}, @qcode{'triangle'}, or
## @qcode{'epanechnikov'}.
##
## @item @qcode{'Function'} @tab The function to estimate: @qcode{'pdf'}
## (default) or @qcode{'cdf'}.
##
## @item @qcode{'Weights'} @tab A vector of non-negative weights, one for each
## row of @var{x}.  The default weights are all equal.
## @end multitable
##
## @seealso{ksdensity}
## @end deftypefn

function f = mvksdensity (x, pts, varargin)

  if (nargin < 2)
    print_usage ();
  endif
  if (! (isnumeric (x) && isreal (x) && ismatrix (x) && ndims (x) == 2))
    error ("mvksdensity: X must be a matrix of real values.");
  endif
  n = rows (x);
  d = columns (x);
  if (n < 2)
    error ("mvksdensity: X must contain at least two observations.");
  endif
  if (! (isnumeric (pts) && isreal (pts) && ismatrix (pts) && ndims (pts) == 2))
    error ("mvksdensity: PTS must be a matrix of real values.");
  endif
  if (columns (pts) != d)
    error ("mvksdensity: PTS must have the same number of columns as X.");
  endif

  ## Defaults and Name-Value parsing.
  kernel = 'normal';
  bw = [];
  func = 'pdf';
  weights = [];
  if (mod (numel (varargin), 2) != 0)
    error ("mvksdensity: optional arguments must be Name-Value pairs.");
  endif
  for k = 1:2:numel (varargin)
    if (! ischar (varargin{k}))
      error ("mvksdensity: parameter names must be character vectors.");
    endif
    switch (lower (varargin{k}))
      case 'bandwidth'
        bw = varargin{k+1};
      case 'kernel'
        kernel = varargin{k+1};
      case 'function'
        func = lower (varargin{k+1});
      case 'weights'
        weights = varargin{k+1};
      otherwise
        error ("mvksdensity: unknown parameter name '%s'.", varargin{k});
    endswitch
  endfor

  ## Validate the kernel.
  if (! (ischar (kernel) && any (strcmpi (kernel, ...
         {'normal', 'box', 'triangle', 'epanechnikov'}))))
    error ("mvksdensity: unrecognised 'Kernel' value.");
  endif
  kernel = lower (kernel);
  if (! any (strcmp (func, {'pdf', 'cdf'})))
    error ("mvksdensity: unrecognised 'Function' value.");
  endif

  ## Bandwidth: normal-reference rule per dimension unless supplied.
  if (isempty (bw))
    bw = default_bw (x, n, d);
  elseif (isnumeric (bw) && isreal (bw) && isscalar (bw) && bw > 0)
    bw = repmat (bw, 1, d);
  elseif (! (isnumeric (bw) && isreal (bw) && isvector (bw)
             && numel (bw) == d && all (bw > 0)))
    error (strcat ("mvksdensity: 'Bandwidth' must be a positive scalar", ...
                   " or a vector with one element per column of X."));
  else
    bw = bw(:)';
  endif

  ## Weights (normalised to sum to one).
  if (isempty (weights))
    w = ones (n, 1) / n;
  elseif (isnumeric (weights) && isreal (weights) && isvector (weights)
          && numel (weights) == n && all (weights >= 0) && any (weights > 0))
    w = weights(:) / sum (weights);
  else
    error (strcat ("mvksdensity: 'Weights' must be a non-negative vector", ...
                   " with one element per observation in X."));
  endif

  ## Evaluate the product-kernel estimate at each query point.
  m = rows (pts);
  f = zeros (m, 1);
  ispdf = strcmp (func, 'pdf');
  for j = 1:m
    ## Standardized distance from every observation to this query point.
    z = (pts(j,:) - x) ./ bw;
    if (ispdf)
      ## Product of the per-dimension kernels, normalised once by prod (bw).
      k = prod_dim (mvks_kpdf (z, kernel), d) ./ prod (bw);
      f(j) = sum (w .* k);
    else
      c = prod_dim (mvks_kcdf (z, kernel), d);
      f(j) = sum (w .* c);
    endif
  endfor

endfunction

## Diagonal normal-reference (Silverman) bandwidth, one element per dimension.
## The per-column spread is a robust estimate of the standard deviation (the
## same rule ksdensity uses), falling back to the ordinary standard deviation
## for any degenerate column, matching MATLAB.
function bw = default_bw (x, n, d)
  sigma = median (abs (x - median (x, 1)), 1) / 0.6745;
  bad = ! (sigma > 0);
  if (any (bad))
    s = std (x, 0, 1);
    sigma(bad) = s(bad);
  endif
  sigma(! (sigma > 0)) = 1;
  bw = sigma * (4 / ((d + 2) * n)) ^ (1 / (d + 4));
endfunction

## Row-wise product across the D columns of a per-dimension kernel matrix.
function p = prod_dim (k, d)
  if (d == 1)
    p = k;
  else
    p = prod (k, 2);
  endif
endfunction

## Product-kernel per-dimension unit-variance kernel density at Z (NxD).
function k = mvks_kpdf (z, kernel)
  s = mvks_kscale (kernel);
  z = s * z;
  switch (kernel)
    case 'normal'
      k = exp (-0.5 * z .^ 2) / sqrt (2 * pi);
    case 'box'
      k = 0.5 * (abs (z) <= 1);
    case 'triangle'
      k = max (1 - abs (z), 0);
    case 'epanechnikov'
      k = 0.75 * max (1 - z .^ 2, 0);
  endswitch
  k = s * k;
endfunction

## Product-kernel per-dimension unit-variance kernel cdf at Z (NxD).
function c = mvks_kcdf (z, kernel)
  z = mvks_kscale (kernel) * z;
  switch (kernel)
    case 'normal'
      c = 0.5 * erfc (-z / sqrt (2));
    case 'box'
      c = min (max ((z + 1) / 2, 0), 1);
    case 'triangle'
      zc = max (min (z, 1), -1);
      c = (zc < 0) .* (0.5 * (zc + 1) .^ 2) ...
          + (zc >= 0) .* (0.5 + zc - 0.5 * zc .^ 2);
    case 'epanechnikov'
      zc = max (min (z, 1), -1);
      c = 0.75 * zc - 0.25 * zc .^ 3 + 0.5;
  endswitch
endfunction

## Standardized-distance scale of a named kernel (its canonical std).
function s = mvks_kscale (kernel)
  switch (kernel)
    case 'box'
      s = 1 / sqrt (3);
    case 'triangle'
      s = 1 / sqrt (6);
    case 'epanechnikov'
      s = 1 / sqrt (5);
    otherwise
      s = 1;
  endswitch
endfunction

%!demo
%! ## Bivariate kernel density estimate over a grid, drawn as a contour plot.
%! x = [randn(60, 2); randn(40, 2) + 3];
%! [gx, gy] = meshgrid (linspace (-4, 7, 60));
%! f = mvksdensity (x, [gx(:), gy(:)]);
%! contourf (gx, gy, reshape (f, size (gx)));
%! hold on;  plot (x(:,1), x(:,2), 'k.');  hold off;
%! title ('Bivariate kernel density estimate');

%!shared X, pts
%! X = [1 1; 2 1; 1 2; 3 2; 2 3; 4 3; 3 4; 5 4];
%! pts = [2 2; 3 3; 1 1; 4 4];

%!test  ## MATLAB parity: default (robust normal-reference) bandwidth
%! assert_equal (mvksdensity (X, pts), ...
%!    [0.0569998; 0.0520004; 0.0448930; 0.0382812], 1e-6);

%!test  ## MATLAB parity: fixed vector and scalar bandwidth (normal kernel)
%! assert_equal (mvksdensity (X, pts, "Bandwidth", [1 1]), ...
%!    [0.0589; 0.0535; 0.0474; 0.0395], 1e-3);
%! assert_equal (mvksdensity (X, pts, "Bandwidth", 1), ...
%!    mvksdensity (X, pts, "Bandwidth", [1 1]));

%!test  ## MATLAB parity: box and epanechnikov product kernels
%! assert_equal (mvksdensity (X, pts, "Bandwidth", [1 1], "Kernel", "box"), ...
%!    [0.0521; 0.0417; 0.0313; 0.0313], 1e-3);
%! assert_equal (mvksdensity (X, pts, "Bandwidth", [1 1], ...
%!    "Kernel", "epanechnikov"), [0.0585; 0.0523; 0.0411; 0.0383], 1e-3);

%!test  ## MATLAB parity: cumulative distribution and weighted estimate
%! assert_equal (mvksdensity (X, pts, "Bandwidth", [1 1], "Function", "cdf"), ...
%!    [0.2144; 0.4504; 0.0520; 0.6893], 1e-3);
%! assert_equal (mvksdensity (X, pts, "Bandwidth", [1 1], ...
%!    "Weights", [2 1 1 1 1 1 1 1]), [0.0588; 0.0479; 0.0598; 0.0351], 1e-3);

%!test  ## the density integrates to ~1 over a wide grid
%! [gx, gy] = meshgrid (linspace (-6, 11, 220));
%! f = mvksdensity (X, [gx(:), gy(:)]);
%! dx = gx(1,2) - gx(1,1);
%! assert_equal (sum (f) * dx ^ 2, 1, 1e-2);

## Test input validation
%!error <Invalid call to mvksdensity> mvksdensity (ones (3, 2))
%!error <mvksdensity: X must be a matrix of real values.> ...
%! mvksdensity (ones (2, 2, 2), [1 1])
%!error <mvksdensity: X must contain at least two observations.> ...
%! mvksdensity ([1 2], [1 1])
%!error <mvksdensity: PTS must have the same number of columns as X.> ...
%! mvksdensity ([1 1; 2 2], [1 1 1])
%!error <mvksdensity: optional arguments must be Name-Value pairs.> ...
%! mvksdensity ([1 1; 2 2], [1 1], "Bandwidth")
%!error <mvksdensity: unrecognised 'Kernel' value.> ...
%! mvksdensity ([1 1; 2 2], [1 1], "Kernel", "cosine")
%!error <mvksdensity: unrecognised 'Function' value.> ...
%! mvksdensity ([1 1; 2 2], [1 1], "Function", "icdf")
%!error <mvksdensity: 'Bandwidth' must be a positive scalar or a vector with one element per column of X.> ...
%! mvksdensity ([1 1; 2 2], [1 1], "Bandwidth", [1 2 3])
%!error <mvksdensity: 'Bandwidth' must be a positive scalar or a vector with one element per column of X.> ...
%! mvksdensity ([1 1; 2 2], [1 1], "Bandwidth", -1)
%!error <mvksdensity: 'Weights' must be a non-negative vector with one element per observation in X.> ...
%! mvksdensity ([1 1; 2 2], [1 1], "Weights", [1 2 3])
