## Copyright (C) 2026 Andreas Bertsatos <abertsatos@biol.uoa.gr>
## Copyright (C) 2026 Avanish Salunke <avanishsalunke16@gmail.com>
##
## This file is part of the statistics package for GNU Octave.
##
## Octave is free software; you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or
## (at your option) any later version.
##
## Octave is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with Octave; see the file COPYING. If not,
## see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn  {statistics} {@var{yi} =} makima (@var{x}, @var{y}, @var{xq})
## @deftypefnx {statistics} {@var{yi} =} makima (@var{y}, @var{xq})
## @deftypefnx {statistics} {@var{yi} =} makima (@dots{}, @qcode{"extrap"})
##
## Compute the 1-D Modified Akima piecewise cubic Hermite interpolant of
## sample data @var{x} and @var{y}.
##
## The Modified Akima (MAKIMA) algorithm generates a shape-preserving 
## piecewise cubic interpolant. It differs from standard splines by avoiding
## excessive local undulations and overshoots, and it connects collinear points
## (flat regions) with straight lines. It is particularly well-suited for
## oscillatory data where @code{pchip} might aggressively flatten local
## extrema.
##
## The sample points @var{x} must be a vector of unique values. If @var{x}
## is not sorted, the function will automatically sort it and rearrange
## @var{y} accordingly.
##
## The sample values @var{y} can be a scalar, vector, or an N-dimensional
## array. If @var{y} is an N-dimensional array, the interpolation is
## performed along its last dimension, which must have the same length as
## @var{x}. Complex values for @var{y} are supported.
##
## If query points @var{xq} are provided, the function evaluates the 
## interpolant and returns the interpolated values @var{yi}. By default,
## @code{makima} uses the boundary polynomials to extrapolate for points 
## outside the range of @var{x}. The optional string argument @qcode{"extrap"} 
## is accepted for compatibility with other interpolation functions.
##
## If only @var{x} and @var{y} are provided, the function returns a 
## piecewise polynomial structure @var{pp} that represents the interpolant.
## This structure can be evaluated later at specific query points using
## @code{ppval}.
##
## Evaluating the interpolant at query points outside the domain of @var{x}
## automatically extrapolates using the boundary polynomials.
##
## @seealso{interp1, pchip, spline}
## @end deftypefn

function yi = makima (x, y, xq, varargin)
  if (nargin < 2 || nargin > 4)
    error ("makima: invalid number of inputs");
  endif

  if (nargin == 4 && ! strcmpi (varargin{1}, "extrap"))
    error ("makima: unknown option '%s'", varargin{1});
  endif

  return_pp = (nargin == 2);

  if (! return_pp)
    size_xq = size (xq);
  endif

  x = x(:);
  n = numel (x);

  is_y_vector = isvector (y);
  size_y = size (y);

  if (is_y_vector)
    if (numel (y) != n)
      error ("makima: The number of sample points X, %d, is incompatible with the number of values Y, %d.", n, numel (y));
    endif
    y = y(:);
    dim_y = 1;
    nc = 1;
  else
    dim_y = ndims (y);
    if (size_y(dim_y) != n)
      error ("makima: The number of sample points X, %d, is incompatible with the number of values Y, %d.", n, size_y(dim_y));
    endif

    ## permute the interpolation dimension to be the first
    perm_order = 1:numel (size_y);
    perm_order(dim_y) = 1;
    perm_order(1) = dim_y;

    y = permute (y, perm_order);
    nc = numel (y) / n;
    y = reshape (y, n, nc);
  endif

  if (iscomplex (y))
    if (return_pp)
      ## Build complex pp struct
      pp_real = makima (x, real (y));
      pp_imag = makima (x, imag (y));
      yi = pp_real;
      yi.coefs = pp_real.coefs + 1i * pp_imag.coefs;
    else
      yi = makima (x, real (y), xq, varargin{:}) + 1i * makima (x, imag (y), xq, varargin{:});
    endif
    return;
  endif

  if (! return_pp)
    xqv = xq(:);
    nq = numel (xqv);
  endif

  if (n < 2)
    error ("makima: The first two inputs must have at least two elements.");
  endif

  if (! issorted (x))
    [x, sort_idx] = sort (x);
    y = y(sort_idx, :);
  endif

  math_done = false;
  if (n == 2)
    if (return_pp)
      ## Linear coefficients for 2-point pp struct
      coefs = zeros (nc, 4, class (y));
      coefs(:, 3) = (y(2, :).' - y(1, :).') ./ (x(2) - x(1));
      coefs(:, 4) = y(1, :).';
      
      if (is_y_vector)
        dim_out = 1;
      else
        dim_out = size_y;
        dim_out(dim_y) = [];
      endif
      yi = mkpp (x.', coefs, dim_out);
      return;
    else
      yi = interp1 (x, y, xqv, "linear", "extrap");
      yi = reshape (yi, [nq, nc]); 
      math_done = true;
    endif
  endif

  if (! math_done)
    dx = diff (x);
  
    if (any (dx <= 0))
      error ("makima: The sample points x must be unique.");
    endif
  dy = diff (y);
  m = dy ./ dx;

  m_0  = 2 * m(1, :) - m(2, :);
  m_m1 = 2 * m_0     - m(1, :);
  m_n  = 2 * m(end, :) - m(end-1, :);
  m_n1 = 2 * m_n       - m(end, :);

  m_ext = [m_m1; m_0; m; m_n; m_n1];

  d = zeros (n, nc, class (y));
  k_idx = (1 : n)';

  s_im2 = m_ext(k_idx    , :);
  s_im1 = m_ext(k_idx + 1, :);
  s_i   = m_ext(k_idx + 2, :);
  s_ip1 = m_ext(k_idx + 3, :);

  w1 = abs (s_ip1 - s_i)   + abs (s_ip1 + s_i)   / 2;
  w2 = abs (s_im1 - s_im2) + abs (s_im1 + s_im2) / 2;

  W = w1 + w2;

  numer = (w1 .* s_im1 + w2 .* s_i);
  denom = max (W, eps);
  d = numer ./ denom;

  zero_mask = (W == 0);
  if (any (zero_mask(:)))
    fallback = (s_im1 + s_i) / 2;
    d(zero_mask) = fallback(zero_mask);
  endif

  if (return_pp)
    hseg = dx;
    delta = m;
    
    d0 = d(1:end-1, :);
    d1 = d(2:end, :);
    y0 = y(1:end-1, :);
    
    c3 = (d0 + d1 - 2*delta) ./ (hseg .* hseg);
    c2 = (3*delta - 2*d0 - d1) ./ hseg;
    c1 = d0;
    c0 = y0;
    
    c3_t = c3.';
    c2_t = c2.';
    c1_t = c1.';
    c0_t = c0.';
    
    coefs = [c3_t(:), c2_t(:), c1_t(:), c0_t(:)];
    
    if (is_y_vector)
      dim_out = 1;
    else
      dim_out = size_y;
      dim_out(dim_y) = [];
      if (isempty (dim_out))
        dim_out = 1;
      endif
    endif

    yi = mkpp (x.', coefs, dim_out);
    return;
  endif

  yi = NaN (nq, nc, class (y));

  if (nq > 0)
    idx = lookup (x, xqv);

    idx (idx >= n) = n - 1;
    idx (idx == 0) = 1;

    x_left = x(idx);
    hseg = dx(idx);

    s = xqv - x_left;

    for c = 1:nc
      y0 = y(idx, c);
      y1 = y(idx + 1, c);
      d0 = d(idx, c);
      d1 = d(idx + 1, c);

      delta = (y1 - y0) ./ hseg;
      c2 = (3*delta - 2*d0 - d1) ./ hseg;
      c3 = (d0 + d1 - 2*delta) ./ (hseg.^2);

      yi(:, c) = y0 + s .* (d0 + s .* (c2 + s .* c3));
    endfor

  endif

  endif

  if (! return_pp)
    if (is_y_vector)
      yi = reshape (yi, size_xq);
    else
      out_shape = size_y;
      out_shape(dim_y) = nq;
      
      yi = reshape (yi, out_shape(perm_order));
      yi = ipermute (yi, perm_order);
      
      ## only append size_xq if multi-dimensional array is given
      if (! isvector (xq))
        final_shape = size_y;
        final_shape(dim_y) = [];
        final_shape = [final_shape, size_xq];
        
        yi = reshape (yi, final_shape);
      endif
    endif
  endif

endfunction


%!test
%! ## Basic linear-like data
%! x = [1; 2; 3; 4];
%! y = [2; 4; 6; 8];
%! xi = [1.5; 2.5; 3.5];
%! yi = makima (x, y, xi);
%! assert (yi, [3; 5; 7], 1e-12);
%!test
%! ## Nonlinear dataset (finite check)
%! x = [0; 1; 2; 3; 4];
%! y = [0; 1; 0; 1; 0];
%! xi = linspace (0,4,20)';
%! yi = makima (x, y, xi);
%! assert (all (isfinite (yi)));
%!test
%! ## pp structure output
%! x = [1; 2; 3; 4];
%! y = [2; 4; 6; 8];
%! pp = makima (x, y);
%! assert (isstruct (pp));
%! assert (strcmp (pp.form, "pp"));
%! assert (pp.pieces, 3);
%! assert (pp.order, 4);
%!test
%! ## Matrix y input.
%! x = [1; 3; 5];
%! y = [1 3 2; 2 4 6]; 
%! xi = 2;
%! yi = makima (x, y, xi);
%! assert (size(yi), [2, 1]);
%! assert (all (isfinite (yi)));
%! assert (yi(1), 2.304086538461538, 1e-12);
%! assert (yi(2), 3.000000000000000, 1e-12);
%!test
%! ## Extrapolation through default method.
%! x = [1; 2; 3];
%! y = [5; 10; 15];
%! xi = [0; 4];
%! yi = makima (x, y, xi);
%! assert (all (isfinite (yi)));
%! assert (yi, [0; 20], 1e-12);
%!test
%! ## Complex interpolation.
%! x = [1; 2; 4];
%! y = [1+2i; 2+3i; 4+8i];
%! xi = 3;
%! yi = makima (x, y, xi);
%! assert (yi, 3 + 5.09767206477733i, 1e-12);
%! assert (iscomplex (yi));
%!test
%! ## Two-point interpolation.
%! x = [1; 5];
%! y = [10; 30];
%! xi = 3;
%! yi = makima (x, y, xi);
%! assert (yi, 20, 1e-12);
%!test
%! ## Single Precision Input.
%! x = single ([1; 2; 3]);
%! y = single ([10; 20; 30]);
%! xi = single (1.5);
%! yi = makima (x, y, xi);
%! assert (isa (yi, "single"));
%! assert (yi, single (15), 1e-6);
%!test
%! ## Row vector inputs.
%! x = [1 2 3];
%! y = [4 5 6];
%! xi = [1.5 2.5];
%! yi = makima (x, y, xi);
%! assert (yi, [4.5 5.5], 1e-12);
%!test
%! ## Step function.
%! x = [1 2 3 4 5 6];
%! y = [0 0 1 1 0 0];
%! xi = [2.5 3.5 4.5];
%! yi = makima (x, y, xi);
%! expected_11 = [0.5000, 1.1250, 0.5000];
%! assert (yi, expected_11, 1e-12);
%!test
%! ## Runge function (Oscillation Check)
%! x = linspace (-1, 1, 7)';
%! y = 1 ./ (1 + 25 * x.^2);
%! xi = [-0.5; 0.1; 0.5];
%! yi = makima (x, y, xi);
%! expected_12 = [0.148690385982729; 0.857734549516009; 0.148690385982729];
%! assert (yi, expected_12, 1e-12);
%!test
%! ## Constant Slopes / Zero Weights
%! x = [1; 2; 3; 4; 5];
%! y = [1; 1; 1; 1; 1];
%! xi = 3.5;
%! yi = makima (x, y, xi);
%! expected_13 = [1];
%! assert (yi, expected_13, 1e-12);
%!test
%! ## Empty xq input
%! x = [1; 2; 3];
%! y = [4; 5; 6];
%! xi = [];
%! yi = makima (x, y, xi);
%! assert (isempty (yi));
%! assert (! (iscolumn (yi)));
%!test
%! ## Wide range of y-values
%! x = [1e-10; 2e-10; 3e-10; 4e-10];
%! y = [1e10; 2e10; 3e10; 4e10];
%! xi = 2.5e-10;
%! yi = makima (x, y, xi);
%! assert (yi, 2.5e10, 1e-12); 
%!test
%! ## Single column matrix input.
%! x = [1; 2; 3];
%! y = [10; 20; 30];
%! xi = [1.5 2.5]; % Row input
%! yi = makima (x, y, xi);
%! assert (yi, [15 25], 1e-12);
%! assert (isrow (yi));
%!test
%! ## Evaluate pp structure with ppval
%! x = [1; 2; 3; 4];
%! y = [2; 4; 6; 8];
%! xi = [1.5; 2.5; 3.5];
%! pp = makima (x, y);
%! yi_ppval = ppval (pp, xi);
%! yi_direct = makima (x, y, xi);
%! assert (yi_ppval, yi_direct, 1e-12);
%!test
%! ## xq is a 2x2 matrix
%! x = [1; 2; 3; 4; 5];
%! y = [10; 20; 15; 5; 25];
%! xq = [1.5, 2.5; 3.5, 4.5];
%! yi = makima (x, y, xq);
%! assert (size (yi), [2, 2]);
%! expected = [16.85897435897436, 18.22916666666667; 9.81182795698925, 10.84522332506203];
%! assert (yi, expected, 1e-12);
%!test
%! ## xq is a 3D array
%! x = [1; 2; 3; 4; 5];
%! y = [10; 20; 15; 5; 25];
%! xq = ones (2, 2, 2) * 2.5;
%! yi = makima (x, y, xq);
%! assert (size (yi), [2, 2, 2]);
%!test
%! ## pp structure with matrix y input
%! x = [1; 3; 5];
%! y = [1 3 2; 2 4 6]; 
%! pp = makima (x, y);
%! assert (isstruct (pp));
%! assert (pp.pieces, 2);
%! assert (pp.dim, 2);
%! yi_ppval = ppval (pp, 2);
%! yi_direct = makima (x, y, 2);
%! assert (yi_ppval, yi_direct, 1e-12);
%!test
%! ## y is a 3D array [2x3x4] and x is length 4
%! x = [1, 2, 3, 4];
%! xq = [1.5, 2.5, 3.5];
%! y3 = reshape (1:24, [2, 3, 4]);
%! yi = makima (x, y3, xq);
%! assert (size (yi), [2, 3, 3]);
%!test
%! ## Unsorted 'x' inputs 
%! x_unsorted = [3; 1; 2; 4];
%! y_unsorted = [9; 1; 4; 16];
%! xq = [1.5; 2.5];
%! x_sorted = [1; 2; 3; 4];
%! y_sorted = [1; 4; 9; 16];
%! assert (makima (x_unsorted, y_unsorted, xq), makima (x_sorted, y_sorted, xq), 1e-12);
%!test
%! ## Complex piecewise polynomial (pp) structure
%! x = [1 2 3];
%! y = [1 4 9] + 1i * [2 8 18];
%! pp = makima (x, y);
%! assert (isstruct (pp));
%! assert (iscomplex (pp.coefs));
%! assert (ppval (pp, 1.5), makima (x, y, 1.5), 1e-12);
%!test
%! ## N-dimensional y (3D) with N-dimensional xq (2x2 matrix)
%! x = [1 2 3 4];
%! y3 = reshape (1:24, [2 3 4]);
%! xq = [1.5 2.5; 3.5 1.5];
%! yi = makima (x, y3, xq);
%! assert (size (yi), [2 3 2 2]);
%!test
%! ## 2-point pp struct 
%! x = [1; 5];
%! y = [10; 30];
%! pp = makima (x, y);
%! assert (pp.pieces, 1);
%! assert (pp.order, 4);
%! assert (ppval (pp, 3), 20, 1e-12);
%!test
%! ## Exact Collinearity
%! x = [1 2 3 4];
%! y = [2 4 6 8];
%! xi = 2.5;
%! yi = makima (x, y, xi);
%! assert (yi, 5, 1e-12);
%!error <makima: The sample points x must be unique.> makima ([1 1 2], [3 4 5], 1.5)
%!error <makima: invalid number of inputs> makima (1)
%!error <makima: The first two inputs must have at least two elements.> makima (1, 2, 1.5)
%!error <makima: The number of sample points X, 4, is incompatible with the number of values Y, 5.> makima ([1 2 3 4], [1 2 3 4 5], 2)
%!error <makima: unknown option 'linear'> makima ([1 2 3], [1 2 3], 2, "linear")
