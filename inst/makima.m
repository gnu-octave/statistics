## Copyright (C) 2025 Avanish Salunke <avanishsalunke16@gmail.com>
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
## @deftypefn {Function File} {@var{yi} =} makima (@var{x}, @var{y}, @var{xq})
## @deftypefnx {Function File} {@var{yi} =} makima (@var{y}, @var{xq})
## @deftypefnx {Function File} {@var{yi} =} makima (@dots{}, @qcode{"extrap"})
##
## Compute the 1-D Modified Akima (MAKIMA) interpolant of sample data
## @var{x}, @var{y} evaluated at query points @var{xq}.
##
## The MAKIMA method generates a **shape-preserving piecewise cubic Hermite**
## interpolant using a modified Akima derivative estimation method.
##
## For a strictly increasing vector of sample points @var{x} and corresponding
## values @var{y}, the function returns the interpolated values @var{yi} at the
## query points @var{xq}.
##
## If only two inputs are given, @var{x} is assumed to be $1:\text{length}(\var{y})$.
##
## The optional string $\qcode{"extrap"}$ enables extrapolation for query points
## outside the range of @var{x}. By default, extrapolation returns $\text{NaN}$.
##
## @seealso{interp1, pchip, spline}
## @end deftypefn


function yi = makima (x, y, xq, varargin)

  if (nargin == 2)
    xq = y;
    y = x;
    x = (1:numel (y)).';
  elseif (nargin != 3 && nargin != 4)
    error ("makima: invalid number of inputs");
  endif

  is_row_input = isrow (xq);

  extrapolate = false;
  if (nargin == 4)
    if (ischar (varargin{1}) && strcmpi (varargin{1}, "extrap"))
      extrapolate = true;
    else
      error ("makima: unknown option '%s'", varargin{1});
    endif
  endif

  x = x(:);
  n = numel (x);

  if (isvector (y))
    y = y(:);
  endif
  [ry, nc] = size (y);
  if (ry != n)
    error ("makima: x and y must have the same number of rows");
  endif

  if (iscomplex (y))
    yi = makima (x, real (y), xq, varargin{:}) + 1i * makima (x, imag (y), xq, varargin{:});
    if (is_row_input && nc == 1)
      yi = yi.';
    endif
    return;
  endif

  if (n == 1)
    yi = repmat (y(1, :), numel (xq(:)), 1);
    if (is_row_input && nc == 1)
      yi = yi.';
    endif
    return;
  elseif (n == 2)
    yi = interp1 (x, y, xq, "linear", extrapolate);
    if (is_row_input && nc == 1 && iscolumn (yi))
       yi = yi.';
    endif
    return;
  endif

  dx = diff (x);
  if (any (dx <= 0))
    error ("makima: x must be strictly increasing");
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

  xqv = xq(:);
  nq = numel (xqv);
  yi = NaN (nq, nc, class (y));

  if (nq > 0)
    idx = lookup (x, xqv);
    left_mask  = (xqv < x(1));
    right_mask = (xqv > x(end));

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

    if (! extrapolate)
      yi(left_mask, :)  = NaN;
      yi(right_mask, :) = NaN;
    endif
  endif

  if (is_row_input && nc == 1)
    yi = yi.';
  endif

endfunction


%!test
%! % 1. Basic linear-like data
%! x = [1; 2; 3; 4];
%! y = [2; 4; 6; 8];
%! xi = [1.5; 2.5; 3.5];
%! yi = makima (x, y, xi);
%! assert (yi, [3; 5; 7], 1e-12);

%!test
%! % 2. Nonlinear dataset (finite check)
%! x = [0; 1; 2; 3; 4];
%! y = [0; 1; 0; 1; 0];
%! xi = linspace (0,4,20)';
%! yi = makima (x, y, xi);
%! assert (all (isfinite (yi)));

%!test
%! % 3. makima(y, xq) syntax
%! y = [10; 20; 30];
%! xq = [1; 2; 3];
%! yi = makima (y, xq);
%! assert (yi, y, 1e-12);

%!test
%! % 4. Matrix y input (multiple columns)
%! x = [1; 3; 5];
%! y = [1 2; 3 4; 2 6];
%! xi = 2;
%! yi = makima (x, y, xi);
%! assert (columns(yi), 2);
%! assert (all (isfinite (yi)));
%! % Column 1 exact value (1917/832)
%! assert (yi(1), 2.3040865384615385, 1e-12);
%! % Column 2 exact value (Linear)
%! assert (yi(2), 3.000, 1e-12);

%!test
%! % 5. Extrapolation enabled
%! x = [1; 2; 3];
%! y = [5; 10; 15];
%! xi = [0; 4];
%! yi = makima (x, y, xi, "extrap");
%! assert (all (isfinite (yi)));
%! assert (yi, [0; 20], 1e-12);

%!test
%! % 6. Default NaN extrapolation
%! x = [1; 2; 3];
%! y = [5; 10; 15];
%! xi = [0; 4];
%! yi = makima (x, y, xi);
%! assert (isnan (yi(1)));
%! assert (isnan (yi(2)));

%!test
%! % 7. Complex interpolation
%! x = [1; 2; 4];
%! y = [1+2i; 2+3i; 4+8i];
%! xi = 3;
%! yi = makima (x, y, xi);
%! assert (yi, 3 + 5.09767206477733i, 1e-12);
%! assert (iscomplex (yi));

%!test
%! % 8. Single Point Input (1x1)
%! x = 5;
%! y = 12;
%! xi = 5;
%! yi = makima (x, y, xi);
%! assert (yi, 12, 1e-12);

%!test
%! % 9. Two-point interpolation (Linear fallback)
%! x = [1; 5];
%! y = [10; 30];
%! xi = 3;
%! yi = makima (x, y, xi);
%! assert (yi, 20, 1e-12);

%!test
%! % 10. Single Precision Input
%! x = single ([1; 2; 3]);
%! y = single ([10; 20; 30]);
%! xi = single (1.5);
%! yi = makima (x, y, xi);
%! assert (isa (yi, "single"));
%! assert (yi, single (15), 1e-6);

%!test
%! % 11. Row vector inputs (Orientation check)
%! x = [1 2 3];
%! y = [4 5 6];
%! xi = [1.5 2.5];
%! yi = makima (x, y, xi);
%! assert (yi, [4.5 5.5], 1e-12); 

%!test
%! % 12. Step function (Overshoot Check)
%! x = [1 2 3 4 5 6];
%! y = [0 0 1 1 0 0];
%! xi = [2.5 3.5 4.5];
%! yi = makima (x, y, xi);
%! expected_12 = [0.5000, 1.1250, 0.5000];
%! assert (yi, expected_12, 1e-12);

%!test
%! % 13. Runge function (Oscillation Check)
%! x = linspace (-1, 1, 7)';
%! y = 1 ./ (1 + 25 * x.^2);
%! xi = [-0.5; 0.1; 0.5];
%! yi = makima (x, y, xi);
%! expected_13 = [0.148690385982729; 0.857734549516009; 0.148690385982729];
%! assert (yi, expected_13, 1e-12);

%!test
%! % 14. Constant Slopes / Zero Weights
%! x = [1; 2; 3; 4; 5];
%! y = [1; 1; 1; 1; 1];
%! xi = 3.5;
%! yi = makima (x, y, xi);
%! expected_14 = [1];
%! assert (yi, expected_14, 1e-12);

%!test
%! % 15. Empty xq input
%! x = [1; 2; 3];
%! y = [4; 5; 6];
%! xi = [];
%! yi = makima (x, y, xi);
%! assert (isempty (yi));
%! assert (iscolumn (yi)); % Ensure it defaults to column output

%!test
%! % 16. Wide range of y-values
%! x = [1e-10; 2e-10; 3e-10; 4e-10];
%! y = [1e10; 2e10; 3e10; 4e10];
%! xi = 2.5e-10;
%! yi = makima (x, y, xi);
%! assert (yi, 2.5e10, 1e-6); % Using a slightly relaxed tolerance

%!test
%! % 17. Single column matrix input (ensures nc=1 logic works)
%! x = [1; 2; 3];
%! y = [10; 20; 30];
%! xi = [1.5 2.5]; % Row input
%! yi = makima (x, y, xi);
%! assert (yi, [15 25], 1e-12);
%! assert (isrow (yi));

%!error<makima: x must be strictly increasing>
%! makima ([1 1 2], [3 4 5], 1.5)

%!error<makima: x and y must have the same number of rows>
%! makima ([1;2;3], [1;2], 2)

%!error<makima: invalid number of inputs>
%! makima (1)