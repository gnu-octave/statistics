## Copyright (C) 2008 Arno Onken <asnelt@asnelt.org>
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
## @deftypefn  {statistics} {@var{p} =} mvtcdf (@var{x}, @var{rho}, @var{df})
## @deftypefnx {statistics} {@var{p} =} mvncdf (@var{x_lo}, @var{x_up}, @var{rho}, @var{df})
## @deftypefnx {statistics} {@var{p} =} mvncdf (@dots{}, @var{options})
## @deftypefnx {statistics} {[@var{p}, @var{err}] =} mvncdf (@dots{})
##
## Multivariate Student's t cumulative distribution function (CDF).
##
## @code{@var{p} = mvtcdf (@var{x}, @var{rho}, @var{df})} returns the cumulative
## probability of the multivariate student's t distribution with correlation
## parameters @var{rho} and degrees of freedom @var{df}, evaluated at each row
## of @var{x}.  The rows of the @math{NxD} matrix @var{x} correspond to sample
## observations and its columns correspond to variables or coordinates.  The
## return argument @var{p} is a column vector with the same number of rows as in
## @var{x}.
##
## @var{rho} is a symmetric, positive definite, @math{DxD} correlation matrix.
## @var{dF} is a scalar or a vector with @math{N} elements.
##
## Note: @code{mvtcdf} computes the CDF for the standard multivariate Student's
## t distribution, centered at the origin, with no scale parameters.  If
## @var{rho} is a covariance matrix, i.e. @code{diag(@var{rho})} is not all
## ones, @code{mvtcdf} rescales @var{rho} to transform it to a correlation
## matrix.  @code{mvtcdf} does not rescale @var{x}, though.
##
## The multivariate Student's t cumulative probability at @var{x} is defined as
## the probability that a random vector T, distributed as multivariate normal,
## will fall within the semi-infinite rectangle with upper limits defined by
## @var{x}.
## @itemize
## @item @math{Pr@{T(1)<=X(1), T(2)<=X(2), ... T(D)<=X(D)@}}.
## @end itemize
##
## @code{@var{p} = mvtcdf (@var{x_lo}, @var{x_hi}, @var{rho}, @var{df})} returns
## the multivariate Student's t cumulative probability evaluated over the
## rectangle (hyper-rectangle for multivariate data in @var{x}) with lower and
## upper limits defined by @var{x_lo} and @var{x_hi}, respectively.
##
## @code{[@var{p}, @var{err}] = mvtcdf (@dots{})} also returns an error estimate
## @var{err} in @var{p}.
##
## @code{@var{p} = mvtcdf (@dots{}, @var{options})} specifies the structure,
## which controls specific parameters for the numerical integration used to
## compute @var{p}. The required fields are:
##
## @multitable @columnfractions 0.2 0.05 0.75
## @item @qcode{"TolFun"} @tab @tab Maximum absolute error tolerance.  Default
## is 1e-8 for D < 4, or 1e-4 for D >= 4.
##
## @item @qcode{"MaxFunEvals"} @tab @tab Maximum number of integrand evaluations
## when @math{D >= 4}.  Default is 1e7.  Ignored when @math{D < 4}.
##
## @item @qcode{"Display"} @tab @tab Display options.  Choices are @qcode{"off"}
## (default), @qcode{"iter"}, which shows the probability and estimated error at
## each repetition, and @qcode{"final"}, which shows the final probability and
## related error after the integrand has converged successfully.  Ignored when
## @math{D < 4}.
## @end multitable
##
## @seealso{bvtcdf, mvtpdf, mvtrnd, mvtcdfqmc}
## @end deftypefn

function [p, err] = mvtcdf (varargin)

  ## Check for valid number on input and output arguments
  narginchk (3,5);

  ## Check for 'options' structure and parse parameters or add defaults
  if (isstruct (varargin{end}))
    if (isfield (varargin{end}, "TolFun"))
      TolFun = varargin{end}.TolFun;
    else
      error ("mvtcdf: options structure missing 'TolFun' field.");
    endif
    if (isempty (TolFun) && size (varargin{1}, 2) < 4)
      TolFun = 1e-8;
    elseif (isempty (TolFun) && size (varargin{1}, 2) < 26)
      TolFun = 1e-4;
    endif
    if (isfield (varargin{end}, "MaxFunEvals"))
      MaxFunEvals = varargin{end}.MaxFunEvals;
    else
      error ("mvtcdf: options structure missing 'MaxFunEvals' field.");
    endif
    if (isempty (MaxFunEvals))
      MaxFunEvals = 1e7;
    endif
    if (isfield (varargin{end}, "Display"))
      Display = varargin{end}.Display;
    else
      error ("mvtcdf: options structure missing 'Display' field.");
    endif
    DispOptions = {"off", "final", "iter"};
    if (sum (any (strcmpi (Display, DispOptions))) == 0)
      error ("mvtcdf: 'Display' field in 'options' has invalid value.");
    endif
    rem_nargin = nargin - 1;
  else
    if (size (varargin{1}, 2) < 4)
      TolFun = 1e-8;
    elseif (size (varargin{1}, 2) < 26)
      TolFun = 1e-4;
    endif
    MaxFunEvals = 1e7;
    Display = "off";
    rem_nargin = nargin;
  endif

  ## Check for X of X_lo and X_up
  if (rem_nargin < 4)   # MVTCDF(X_UP,SIGMA,DF)
    x_up_Only = true;
    x_up = varargin{1};

    ## Check for x being a matrix
    if (! ismatrix (x_up))
      error ("mvtcdf: X must be a matrix.");
    endif

    ## Create x_lo according to data type of x_lo
    x_lo = - Inf (size (x_up));
    if isa (x_up, "single")
      x_lo = single (x_lo);
    endif

    ## Get SIGMA and DF arguments
    rho = varargin{2};
    df = varargin{3};

  else                  # MVNCDF(X_LO,X_UP,SIGMA,DF)
    x_up_Only = false;
    x_lo = varargin{1};
    x_up = varargin{2};
    rho = varargin{3};
    df = varargin{4};

    ## Check for x_lo and x_up being matrices of the same size
    ## and that they define increasing limits
    if (! ismatrix (x_lo) || ! ismatrix (x_up))
      error ("mvtcdf: X_LO and X_UP must be matrices.");
    endif
    if (any (size (x_lo) != size (x_up)))
      error ("mvtcdf: X_LO and X_UP must be of the same size.");
    endif
    if (any (any (x_lo > x_up)))
      error ("mvtcdf: X_LO and X_UP must define increasing limits.");
    endif
  endif

  ## Check if data is single or double class
  is_type = "double";
  if (isa (x_up, "single") || isa (x_lo, "single") || ...
      isa (rho, "single") || isa (df, "single"))
    is_type = "single";
  endif

  ## Get size of data
  [n_x, d_x] = size (x_lo);
  if (d_x < 1)
    error ("mvtcdf: too few dimensions in data.");
  endif

  ## Force univariate column vector into a row vector
  if ((d_x == 1) && (size (rho, 1) == n_x))
    x_lo = x_lo';
    x_up = x_up';
    [n_x, d_x] = size (x_up);
  endif

  ## Check rho
  sz = size(rho);
  if (sz(1) != sz(2))
    error ("mvtcdf: correlation matrix RHO is not square.");
  elseif (! isequal (sz, [d_x, d_x]))
    error (strcat (["mvtcdf: correlation matrix RHO does not"], ...
                   [" match dimensions in data."]));
  endif

  ## Standardize rho to correlation if necessary (not the data)
  s = sqrt (diag (rho));
  if (any (s != 1))
    rho = rho ./ (s * s');
  endif

  ## Continue checking rho for being a valid correlation matrix
  [~, err] = cholcov (rho, 0);
  if (err != 0)
    error (strcat (["mvtcdf: correlation matrix RHO must be"], ...
                   [" positive semi-definite."]));
  endif

  ## Check df
  if (! isscalar (df) && ! (isvector (df) && length (df) == n_x))
    error (strcat (["mvtcdf: DF must be a scalar or a vector with"], ...
                   [" the same samples as in data."]));
  endif
  if (any (df <= 0) || ! isreal (df))
    error ("mvtcdf: DF must contain only positive real numbers.");
  endif

  ## Compute the cdf
  if (d_x == 1)
    p = tcdf (x_up, df) - tcdf (x_lo, df);
    if (nargout > 1)
      err = NaN (size (p), is_type);
    endif

  elseif (d_x < 4)
    if (x_up_Only)          # upper limit only
      if (d_x == 2)
        p = bvtcdf (x_up, rho(2), df, TolFun);
      else
        p = tvtcdf (x_up, rho([2 3 6]), df, TolFun);
      endif
    else                    # lower and upper limits present
      ## Compute the probability over the rectangle as sums and differences
      ## of integrals over semi-infinite half-rectangles.  For degenerate
      ## rectangles, force an exact zero by making each piece exactly zero.
      equalLimits = (x_lo == x_up);
      x_lo(equalLimits) = -Inf;
      x_up(equalLimits) = -Inf;
      p = zeros (n_x, 1, is_type);
      for i = 0:d_x
        k = nchoosek (1:d_x, i);
        for j = 1:size (k, 1)
          X = x_up;
          X(:,k(j,:)) = x_lo(:,k(j,:));
          if d_x == 2
            p = p + (-1)^i * bvtcdf (X, rho(2), df, TolFun/4);
          else
            p = p + (-1)^i * tvtcdf (X, rho([2 3 6]), df, TolFun/8);
          endif
        endfor
      endfor
    endif
    if (nargout > 1)
      err = repmat (cast (TolFun, is_type), size (p));
    endif

  elseif (d_x < 26)
    p = zeros (n_x, 1, is_type);
    err = zeros (n_x, 1, is_type);
    if (isscalar (df))
      df = repmat (df, n_x, 1);
    endif
    for i = 1:n_x
      [p(i), err(i)] = mvtcdfqmc (x_lo(i,:), x_up(i,:), rho, df(i), ...
                                  TolFun, MaxFunEvals, Display);
    endfor

  else
    error ("mvncdf: too many dimensions in data (limit = 25 columns).");
  endif

  ## Bound p in range [0, 1]
  p(p < 0) = 0;
  p(p > 1) = 1;

endfunction

## CDF for the trivariate Student's T
function p = tvtcdf (x, rho, df, TolFun)

  n_x = size (x, 1);
  if (isscalar (df))
    df = repmat (df, n_x, 1);
  endif

  ## Find a permutation that makes rho_23 == max(rho)
  [~,imax] = max (abs (rho));
  if (imax == 1)      # swap 1 and 3
    rho_12 = rho(3);
    rho_13 = rho(2);
    rho_23 = rho(1);
    x = x(:,[3 2 1]);
  elseif (imax == 2)  # swap 1 and 2
    rho_12 = rho(1); rho_13 = rho(3); rho_23 = rho(2);
    x = x(:,[2 1 3]);
  else                # x already in correct order
    rho_12 = rho(1);
    rho_13 = rho(2);
    rho_23 = rho(3);
  endif

  if (rho_23 >= 0)
    p1 = bvtcdf ([x(:,1) min(x(:,2:3), [], 2)], 0, df, TolFun / 4);
    p1(any (isnan (x), 2)) = NaN;
  else
    p1 = bvtcdf (x(:,1:2), 0, df, TolFun  /4) - ...
         bvtcdf ([x(:,1) -x(:,3)], 0, df, TolFun / 4);
    p1(p1 < 0) = 0;
  endif

  if (abs (rho_23) < 1)
    lo = asin (rho_23);
    hi = (sign (rho_23) + (rho_23 == 0)) .* pi ./ 2;
    p2 = zeros (size (p1), class (p1));
    for i = 1:n_x
      x1 = x(i,1);
      x2 = x(i,2);
      x3 = x(i,3);
      if (isfinite (x2) && isfinite (x3) && ~! isnan (x1))
        v = df(i);
        p2(i) = quadgk (@tvtIntegr1, lo, hi, "AbsTol", TolFun / 4, "RelTol", 0);
      endif
    endfor
  else
    p2 = zeros (class (p1));
  endif

  if (abs (rho_12) > 0)
    lo = 0;
    hi = asin (rho_12);
    rj = rho_12;
    rk = rho_13;
    p3 = zeros (size (p1), class (p1));
    for i = 1:n_x
      x1 = x(i,1);
      xj = x(i,2);
      xk = x(i,3);
      if (isfinite (x1) && isfinite (xj) && ! isnan (xk))
        v = df(i);
        p3(i) = quadgk (@tvtIntegr2, lo, hi, "AbsTol", TolFun / 4, "RelTol", 0);
      endif
    endfor
  else
    p3 = zeros (class (p1));
  endif

  if (abs (rho_13) > 0)
    lo = 0;
    hi = asin (rho_13);
    rj = rho_13;
    rk = rho_12;
    p4 = zeros (size (p1), class (p1));
    for i = 1:n_x
      x1 = x(i,1);
      xj = x(i,3);
      xk = x(i,2);
      if (isfinite (x1) && isfinite (xj) && ! isnan (xk))
        v = df(i);
        p4(i) = quadgk (@tvtIntegr2, lo, hi, "AbsTol", TolFun / 4, "RelTol", 0);
      endif
    endfor
  else
    p4 = zeros (class (p1));
  endif

  if (isa (x, "single") || isa (rho, "single") || isa (df, "single"))
    p = cast (p1 + (-p2 + p3 + p4) ./ (2 .* pi), "single");
  else
    p = cast (p1 + (-p2 + p3 + p4) ./ (2 .* pi), "double");
  endif

  ## Functions to compute the integrands
  function integrand = tvtIntegr1 (theta)
    st = sin(theta);
    c2t = cos(theta) .^ 2;
    w = sqrt (1 ./ (1 + ((x2 * st - x3) .^ 2 ./ c2t + x2 .^ 2) / v));
    integrand = w .^ v .* TCDF (x1 .* w, v);
  endfunction

  function integrand = tvtIntegr2 (theta)
    st = sin (theta);
    c2t = cos (theta) .^ 2;
    w = sqrt (1 ./ (1 + ((x1 *st - xj) .^ 2 ./ c2t + x1 .^ 2) / v));
    integrand = w .^ v .* TCDF (uk (st, c2t) .* w, v);
  endfunction

  function uk = uk (st, c2t)
    sinphi = st .* rk ./ rj;
    numeru = xk .* c2t - x1 .* (sinphi - rho_23 .* st) ...
                       - xj .* (rho_23 - st .* sinphi);
    denomu = sqrt (c2t .* (c2t - sinphi .* sinphi ...
                               - rho_23 .* (rho_23 - 2 .* st .* sinphi)));
    uk = numeru ./ denomu;
  endfunction

endfunction

## CDF for Student's T
function p = TCDF (x, df)
  p = betainc(df ./ (df + x .^ 2), df / 2, 0.5) / 2;
  reflect = (x > 0);
  p(reflect) = 1 - p(reflect);
endfunction

%!demo
%! ## Compute the cdf of a multivariate Student's t distribution with
%! ## correlation parameters rho = [1, 0.4; 0.4, 1] and 2 degrees of freedom.
%!
%! rho = [1, 0.4; 0.4, 1];
%! df = 2;
%! [X1, X2] = meshgrid (linspace (-2, 2, 25)', linspace (-2, 2, 25)');
%! X = [X1(:), X2(:)];
%! p = mvtcdf (X, rho, df);
%! surf (X1, X2, reshape (p, 25, 25));
%! title ("Bivariate Student's t cumulative distribution function");

## Test output against MATLAB R2018
%!test
%! x = [1, 2];
%! rho = [1, 0.5; 0.5, 1];
%! df = 4;
%! a = [-1, 0];
%! assert (mvtcdf(a, x, rho, df), 0.294196905339283, 1e-14);
%!test
%! x = [1, 2;2, 4;1, 5];
%! rho = [1, 0.5; 0.5, 1];
%! df = 4;
%! p =[0.790285178602166; 0.938703291727784; 0.81222737321336];
%! assert (mvtcdf(x, rho, df), p, 1e-14);
%!test
%! x = [1, 2, 2, 4, 1, 5];
%! rho = eye (6);
%! rho(rho == 0) = 0.5;
%! df = 4;
%! assert (mvtcdf(x, rho, df), 0.6874, 1e-4);

%!error mvtcdf (1)
%!error mvtcdf (1, 2)
%!error<mvtcdf: correlation matrix RHO does not match dimensions in data.> ...
%! mvtcdf (1, [2, 3; 3, 2], 1)
%!error<mvtcdf: correlation matrix RHO does not match dimensions in data.> ...
%! mvtcdf ([2, 3, 4], ones (2), 1)
%!error<mvtcdf: X_LO and X_UP must be of the same size.> ...
%! mvtcdf ([1, 2, 3], [2, 3], ones (2), 1)
%!error<mvtcdf: correlation matrix RHO must be positive semi-definite.> ...
%! mvtcdf ([2, 3], ones (2), [1, 2, 3])
%!error<mvtcdf: DF must be a scalar or a vector with the same samples as in> ...
%! mvtcdf ([2, 3], [1, 0.5; 0.5, 1], [1, 2, 3])

