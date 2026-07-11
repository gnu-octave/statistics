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
## @deftypefn  {statistics} {@var{p} =} dwtest (@var{r}, @var{x})
## @deftypefnx {statistics} {@var{p} =} dwtest (@var{r}, @var{x}, @var{name}, @var{value})
## @deftypefnx {statistics} {[@var{p}, @var{d}] =} dwtest (@dots{})
##
## Durbin-Watson test for autocorrelation in linear regression residuals.
##
## @code{@var{p} = dwtest (@var{r}, @var{x})} performs the Durbin-Watson test on
## the residuals @var{r} of a linear regression with design matrix @var{x}
## (which should include a column of ones if the model has a constant term).
## The null hypothesis is that the residuals are uncorrelated, against the
## alternative that they are autocorrelated.  @var{r} is an @math{N*1} vector and
## @var{x} is an @math{N*P} matrix.  @var{p} is the p-value of the test.
##
## The Durbin-Watson statistic is
## @tex
## $ d = \sum_{i=1}^{n-1} (r_{i+1} - r_i)^2 / \sum_{i=1}^{n} r_i^2 $.
## @end tex
## @ifnottex
## @code{@var{d} = sum ((diff (@var{r})) .^ 2) / sum (@var{r} .^ 2)}.
## @end ifnottex
## Values near 2 indicate no autocorrelation, values towards 0 positive
## autocorrelation, and values towards 4 negative autocorrelation.
##
## @code{@var{p} = dwtest (@var{r}, @var{x}, @var{name}, @var{value})} specifies
## additional options using @qcode{Name-Value} pair arguments:
##
## @multitable @columnfractions 0.18 0.8
## @headitem @var{Name} @tab @var{Value}
##
## @item @qcode{'Method'} @tab @qcode{'exact'} to compute the exact p-value
## from the null distribution of the statistic (a ratio of quadratic forms,
## evaluated with Imhof's method), or @qcode{'approximate'} to use a normal
## approximation based on the mean and variance of the statistic.  The default is
## @qcode{'exact'} for @math{n < 400} and @qcode{'approximate'} otherwise.
##
## @item @qcode{'Tail'} @tab The alternative hypothesis: @qcode{'both'}
## (default) for a nonzero autocorrelation, @qcode{'right'} for a positive
## autocorrelation, or @qcode{'left'} for a negative autocorrelation.
## @end multitable
##
## @code{[@var{p}, @var{d}] = dwtest (@dots{})} also returns the Durbin-Watson
## statistic @var{d}.
##
## @seealso{regress, fitlm, runstest}
## @end deftypefn

function [pval, d] = dwtest (r, x, varargin)

  if (nargin < 2)
    print_usage ();
  endif
  if (! (isnumeric (r) && isreal (r) && isvector (r)))
    error ("dwtest: R must be a real vector of residuals.");
  endif
  r = r(:);
  n = numel (r);
  if (! (isnumeric (x) && isreal (x) && ismatrix (x) && rows (x) == n))
    error ("dwtest: X must be a real matrix with one row per residual.");
  endif

  ## Parse Name-Value pairs
  method = "";
  tail = "both";
  if (mod (numel (varargin), 2) != 0)
    error ("dwtest: optional arguments must be given as Name-Value pairs.");
  endif
  for k = 1:2:numel (varargin)
    if (! ischar (varargin{k}))
      error ("dwtest: parameter names must be character vectors.");
    endif
    switch (lower (varargin{k}))
      case "method"
        method = lower (varargin{k+1});
      case "tail"
        tail = lower (varargin{k+1});
      otherwise
        error ("dwtest: unknown parameter name '%s'.", varargin{k});
    endswitch
  endfor
  if (isempty (method))
    method = ifelse (n < 400, "exact", "approximate");
  endif
  if (! any (strcmp (method, {"exact", "approximate"})))
    error ("dwtest: 'Method' must be 'exact' or 'approximate'.");
  endif
  if (! any (strcmp (tail, {"both", "right", "left"})))
    error ("dwtest: 'Tail' must be 'both', 'right', or 'left'.");
  endif

  ## Durbin-Watson statistic
  d = sum (diff (r) .^ 2) / sum (r .^ 2);

  ## Orthonormal basis of the column space of X (residual-maker M = I - Q*Q')
  [Q, R] = qr (x, 0);
  keep = abs (diag (R)) > max (size (x)) * eps (max (abs (diag (R))));
  Q = Q(:, keep);
  p = columns (Q);

  ## Durbin-Watson difference operator A (so that r'*A*r = sum (diff (r) .^ 2))
  A = 2 * eye (n);
  A(1,1) = 1;
  A(n,n) = 1;
  A -= diag (ones (n-1, 1), 1) + diag (ones (n-1, 1), -1);

  ## Lower-tail probability P(D <= d) under the null hypothesis
  if (strcmp (method, "exact"))
    M = eye (n) - Q * Q';
    MAM = M * A * M;
    nu = sort (eig ((MAM + MAM') / 2), "descend");
    nu = nu(1:n - p);                    ## n - p nonzero eigenvalues of M*A*M
    plow = dwtest_imhof_ (nu, d);
  else
    ## Mean and variance of D from traces, avoiding the eigendecomposition
    AQ = A * Q;
    QtAQ = Q' * A * Q;
    sumnu  = (2 * n - 2) - trace (QtAQ);              ## tr (M*A)
    sumnu2 = (6 * n - 8) - 2 * sum (AQ(:) .^ 2) + sum (QtAQ(:) .^ 2); ## tr((M*A)^2)
    m1 = sumnu / (n - p);
    s2 = 2 * (sumnu2 - sumnu ^ 2 / (n - p)) / ((n - p) * (n - p + 2));
    plow = normcdf ((d - m1) / sqrt (s2));
  endif
  plow = min (max (plow, 0), 1);

  ## Alternative hypothesis: 'right' -> positive autocorrelation (small D),
  ## 'left' -> negative autocorrelation (large D)
  switch (tail)
    case "right"
      pval = plow;
    case "left"
      pval = 1 - plow;
    case "both"
      pval = 2 * min (plow, 1 - plow);
  endswitch
  pval = min (max (pval, 0), 1);

endfunction

## Imhof's method for P(sum ((nu_i - d) * chi2_1) <= 0), which equals the
## lower-tail probability P(D <= d) of the Durbin-Watson statistic.
function plow = dwtest_imhof_ (nu, d)
  lam = nu - d;
  I = quadgk (@(u) dwtest_imhof_integrand_ (u, lam), 0, Inf, ...
              "AbsTol", 1e-10, "RelTol", 1e-9);
  plow = 0.5 - I / pi;
endfunction

function y = dwtest_imhof_integrand_ (u, lam)
  sz = size (u);
  uv = u(:).';                   ## row vector of quadrature nodes
  lam = lam(:);                  ## column vector of coefficients
  theta = 0.5 * sum (atan (lam * uv), 1);
  rho = prod ((1 + (lam .^ 2) * (uv .^ 2)) .^ 0.25, 1);
  y = sin (theta) ./ (uv .* rho);
  y(uv == 0) = 0.5 * sum (lam);  ## limit of the integrand as u -> 0
  y = reshape (y, sz);           ## match the shape quadgk expects
endfunction

function out = ifelse (cond, a, b)
  if (cond)
    out = a;
  else
    out = b;
  endif
endfunction

%!demo
%! ## Test regression residuals for autocorrelation
%! x = [ones(20, 1), (1:20)'];
%! y = x * [1; 0.5] + sin ((1:20)' / 2);   # add an autocorrelated component
%! b = x \ y;
%! r = y - x * b;
%! [p, d] = dwtest (r, x)

## Test the statistic value
%!test
%! x = [ones(6, 1), (1:6)'];
%! r = [1; -1; 1; -1; 1; -1];              # strong negative autocorrelation
%! [p, d] = dwtest (r, x);
%! assert_equal (d, sum (diff (r) .^ 2) / sum (r .^ 2), 1e-12);
%! assert_equal (d, 20 / 6, 1e-12);
%! assert_equal (p >= 0 && p <= 1, true);

%!test  # exact and approximate methods give similar p-values
%! x = [ones(30, 1), (1:30)', ((1:30)') .^ 2];
%! r = sin ((1:30)' / 3);
%! pe = dwtest (r, x, "Method", "exact");
%! pa = dwtest (r, x, "Method", "approximate");
%! assert_equal (pe, pa, 0.05);

%!test  # tail selection is consistent
%! x = [ones(15, 1), (1:15)'];
%! r = (-1) .^ (1:15)';                    # alternating -> D near 4
%! pr = dwtest (r, x, "Tail", "right");
%! pl = dwtest (r, x, "Tail", "left");
%! pb = dwtest (r, x, "Tail", "both");
%! assert_equal (pr + pl, 1, 1e-10);
%! assert_equal (pb, 2 * min (pr, pl), 1e-10);

%!test  # left tail rejects strong negative autocorrelation (D near 4)
%! x = [ones(20, 1), (1:20)'];
%! r = (-1) .^ (1:20)';
%! assert_equal (dwtest (r, x, "Tail", "left") < 0.05, true);

## Test input validation
%!error <Invalid call to dwtest> dwtest (1)
%!error <dwtest: R must be a real vector of residuals.> dwtest (ones (3, 3), ones (3, 2))
%!error <dwtest: X must be a real matrix with one row per residual.> ...
%! dwtest ([1;2;3], ones (2, 2))
%!error <dwtest: optional arguments must be given as Name-Value pairs.> ...
%! dwtest ([1;2;3], ones (3, 1), "Tail")
%!error <dwtest: unknown parameter name 'foo'.> ...
%! dwtest ([1;2;3], ones (3, 1), "foo", "bar")
%!error <dwtest: 'Method' must be 'exact' or 'approximate'.> ...
%! dwtest ([1;2;3], ones (3, 1), "Method", "fast")
%!error <dwtest: 'Tail' must be 'both', 'right', or 'left'.> ...
%! dwtest ([1;2;3], ones (3, 1), "Tail", "up")
