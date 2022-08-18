## Copyright (C) 2022 Andreas Bertsatos <abertsatos@biol.uoa.gr>
## Copyright (C) 2008 Arno Onken <asnelt@asnelt.org>
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
## @deftypefn {Function File} @var{y} = mvncdf (@var{x})
## @deftypefnx {Function File} @var{y} = mvncdf (@var{x}, @var{mu}, @var{sigma})
## @deftypefnx {Function File} @var{y} = mvncdf (@var{x_lo}, @var{x_up}, @var{mu}, @var{sigma})
## @deftypefnx {Function File} @var{y} = mvncdf (@dots{}, @var{options})
## @deftypefnx {Function File} [@var{y}, @var{err}] = mvncdf (@dots{})
##
## Multivariate normal cumulative distribution function.
##
## @code{@var{y} = mvncdf (@var{x})} returns cumulative probability of the
## multivariate normal distribution evaluated at each row of @var{x} with zero
## mean and an identity covariance matrix. The row of matrix @var{x} correspond
## to observations and its columns to variables. The return argument @var{y} is
## a column vector with the same number of rows as in @var{x}.
##
## @code{@var{y} = mvncdf (@var{x}, @var{mu}, @var{sigma})} returns cumulative
## probability of the multivariate normal distribution evaluated at each row of
## @var{x} with mean @var{mu} and a covariance matrix @var{sigma}.  @var{mu} can
## be either a scalar (the same of every variable) or a row vector with the same
## number of elements as the number of variables in @var{x}.  @var{sigma}
## covariance matrix may be specified a row vector if it only contains variances
## along its diagonal and zero covariances of the diagonal.  In such a case, the
## diagonal vector @var{sigma} must have the same number of elements as the
## number of variables (columns) in @var{x}.  If you only want to specify sigma,
## you can pass an empty matrix for @var{mu}.
##
## The multivariate normal cumulative probability at @var{x} is defined as the
## probability that a random vector V, distributed as multivariate normal, will
## fall within the semi-infinite rectangle with upper limits defined by @var{x}.
## @itemize
## @item Pr@{V(1)<=X(1), V(2)<=X(2), ... V(D)<=X(D)@}.
## @end itemize
##
## @code{@var{y} = mvncdf (@var{x_lo}, @var{x_hi}, @var{mu}, @var{sigma})}
## returns the multivariate normal cumulative probability evaluated over the
## rectangle (hyper-rectangle for multivariate data in @var{x}) with lower and
## upper limits defined by @var{x_lo} and @var{x_hi} respectively.
##
## @code{[@var{y}, @var{err}] = mvncdf (@dots{})} also returns an error estimate
## @var{err} in @var{y}.
##
## @code{@var{y} = mvncdf (@dots{}, @var{options})} specifies the structure,
## which controls specific parameters for the numerical integration for
## numltivariate cases. The required fieds are:
##
## @multitable @columnfractions 0.2 0.8
## @item "TolFun" @tab --- Maximum absolute error tolerance.  Default is
## 1e-8 for D == 1 | 3, 1e-4 for D > 4. Note that for bivariate normal cdf, the
## Octave implementation has a presicion of more than 1e-10.
## @item "MaxFunEvals" @tab --- Maximum number of integrand evaluations.
## Default is 1e7 for D > 4.
## @item "Display" @tab --- Display options.  Choices are "off" (default),
## "iter", which shows the probability and estimated error at each repetition,
## and "final", which shows the final probability and related error after the
## integrand has converged successfully.
## @end multitable
##
## @end deftypefn

function [y, err] = mvncdf (varargin)
  ## Check for valid number on input and output arguments
  narginchk (1,5);
  nargoutchk (1,2);
  ## Check for 'options' structure and parse parameters or add defaults
  if (isstruct (varargin{end}))
    if (isfield (varargin{end}, "TolFun"))
      TolFun = varargin{end}.TolFun;
    else
      error ("mvncdf: options structure missing 'TolFun' field.");
    endif
    if (isempty (TolFun) && size (varargin{1}, 2) < 4)
      TolFun = 1e-8;
    elseif (isempty (TolFun) && size (varargin{1}, 2) < 26)
      TolFun = 1e-4;
    endif
    if (isfield (varargin{end}, "MaxFunEvals"))
      MaxFunEvals = varargin{end}.MaxFunEvals;
    else
      error ("mvncdf: options structure missing 'MaxFunEvals' field.");
    endif
    if (isempty (MaxFunEvals))
      MaxFunEvals = 1e7;
    endif
    if (isfield (varargin{end}, "Display"))
      Display = varargin{end}.Display;
    else
      error ("mvncdf: options structure missing 'Display' field.");
    endif
    DispOptions = {"off", "final", "iter"};
    if (sum (any (strcmpi (Display, DispOptions))) == 0)
      error ("mvncdf: 'Display' field in 'options' has invalid value.");
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
  if (rem_nargin < 4)   # MVNCDF(XU,MU,SIGMA)
    x_up_Only = true;
    x_up = varargin{1};
    ## Check for x being a matrix
    if (! ismatrix (x_up))
      error ("mvncdf: X must be a matrix.");
    endif
    ## Create x_lo according to data type of x_lo
    x_lo = - Inf (size (x_up));
    if isa (x_up, "single")
      x_lo = single (x_lo);
    endif
    ## Check for mu and sigma arguments
    if (rem_nargin > 1)
      mu = varargin{2};
    else
      mu = [];
    endif
    if (rem_nargin > 2)
      sigma = varargin{3};
    else
      sigma = [];
    endif
  else                  ## MVNCDF(XL,XU,MU,SIGMA)
    x_up_Only = false;
    x_lo = varargin{1};
    x_up = varargin{2};
    mu = varargin{3};
    sigma = varargin{4};
    ## Check for x_lo and x_up being matrices of the same size
    ## and that they define increasing limits
    if (! ismatrix (x_lo) || ! ismatrix (x_up))
      error ("mvncdf: X_lo and x_up must be matrices.");
    endif
    if (size (x_lo) != size (x_up))
      error ("mvncdf: X_lo and x_up must have the same size.");
    endif
    if (any (any (x_lo > x_up)))
      error ("mvncdf: X_lo and x_up must define increasing limits.");
    endif
  endif
  ## Check if data is single or double class
  is_type = "double";
  if (isa (x_lo, "single"))
    is_type = "single";
  endif
  ## Get size of data
  [n_x, d_x] = size (x_lo);
  ## Center data according to mu
  if (isempty (mu))         # already centered
    XLo0 = x_lo;
    XUp0 = x_up;
  elseif (isscalar (mu))    # mu is a scalar
    XLo0 = x_lo - mu;
    XUp0 = x_up - mu;
  elseif (isvector (mu))    # mu is a vector
    ## Get size of mu vector
    [n_mu, d_mu] = size (mu);
    if (d_mu != d_x)
      error ("mvncdf: wrong size of 'mu' vector.");
    endif
    if (n_mu == 1 || n_mu == n_x)
      XLo0 = x_lo - mu;
      XUp0 = x_up - mu;
    else
      error ("mvncdf: wrong size of 'mu' vector.");
    endif
  else
    error ("mvncdf: 'mu' must be either empty, a scalar, or a vector.");
  endif
  ## Check how sigma was parsed
  if (isempty (sigma))        # already standardized
    ## If x_lo and x_up are column vectors, transpose them to row vectors
    if (d_x == 1)
      XLo0 = XLo0';
      XUp0 = XUp0';
      [n_x, d_x] = size (XUp0);
    endif
    sigmaIsDiag = true;
    sigma = ones (1, d_x);
  else
    ## Check if sigma parsed as diagonal vector
    if (size (sigma, 1) == 1 && size (sigma, 2) > 1)
      sigmaIsDiag = true;
    else
      sigmaIsDiag = false;
    endif
    ## If x_lo and x_up are column vectors, transpose them to row vectors
    if (d_x == 1)
      if (isequal (size (sigma), [1, n_x]))
        XLo0 = XLo0';
        XUp0 = XUp0';
        [n_x, d_x] = size (XUp0);
      elseif (! isscalar (mu))
        error ("mvncdf: 'mu' must a scalar if sigma is a vector.");
      endif
    endif
    ## Check for sigma being a valid covariance matrix
    if (! sigmaIsDiag && (size (sigma, 1) != size (sigma, 2)))
      error ("mvncdf: covariance matrix is not symmetric.");
    elseif (! sigmaIsDiag && (! all (size (sigma) == [d_x, d_x])))
      error ("mvncdf: covariance matrix does not match x_up.");
    else
      ## If sigma is a covariance matrix check that it is positive semi-definite
      if (! sigmaIsDiag)
        [~, err] = chol (sigma);
        if (err != 0)
          error ("mvncdf: covariance matrix must be positive semi-definite.");
        endif
      else
        if (any (sigma) <= 0)
          error ("mvncdf: invalid sigma diagonal vector.");
        endif
      endif
    endif
  endif
  ## Standardize sigma and x data
  if (sigmaIsDiag)
    XLo0 = XLo0 ./ sqrt (sigma);
    XUp0 = XUp0 ./ sqrt (sigma);
  else
    s = sqrt (diag (sigma))';
    XLo0 = XLo0 ./ s;
    XUp0 = XUp0 ./ s;
    Rho = sigma ./ (s * s');
  endif
  ## Compute the cdf from standardized values.
  if (d_x == 1)
    y = normcdf (XUp0, 0, 1) - normcdf (XLo0, 0, 1);
    if (nargout > 1)
      err = NaN (size(y), is_type);
    endif
  elseif (sigmaIsDiag)
    y = prod (normcdf (XUp0, 0, 1) - normcdf (XLo0, 0, 1), 2);
    if (nargout > 1)
      err = NaN (size(y), is_type);
    endif
  elseif (d_x < 4)
    if (x_up_Only)          # upper limit only
      if (d_x == 2)
        y = bvncdf (x_up, mu, sigma);
      else
        y = tvncdf (XUp0, Rho([2 3 6]), TolFun);
      endif
    else                    # lower and upper limits present
      ## Compute the probability over the rectangle as sums and differences
      ## of integrals over semi-infinite half-rectangles.  For degenerate
      ## rectangles, force an exact zero by making each piece exactly zero.
      equalLimits = (XUp0 == XLo0);
      XUp0(equalLimits) = -Inf;
      XLo0(equalLimits) = -Inf;
      ## For bvncdf
      x_up(equalLimits) = -Inf;
      x_lo(equalLimits) = -Inf;
      y = zeros(n_x, 1, is_type);
      for i = 0:d_x
        k = nchoosek (1:d_x, i);
        for j = 1:size (k, 1)
          X = XUp0;
          X(:,k(j,:)) = XLo0(:,k(j,:));
          if d_x == 2
            x = x_up;
            x(:,k(j,:)) = x_lo(:,k(j,:));
            y = y + (-1) ^ i * bvncdf (x, mu, sigma);
          else
            y = y + (-1) ^ i * tvncdf (X, Rho([2 3 6]), TolFun / 8);
          endif
        endfor
      endfor
    endif
    if (nargout > 1)
      err = repmat (cast (TolFun, is_type), size (y));
    endif
  elseif (d_x < 26)
    y = zeros (n_x, 1, is_type);
    err = zeros (n_x, 1, is_type);
    for i = 1:n_x
      [y(i), err(i)] = mvtcdfqmc (XUp0(i,:), XUp0(i,:), Rho, Inf, ...
                                  TolFun, MaxFunEvals, Display);
    endfor
  else
    error ("mvncdf: too many dimensions in data (limit = 25 columns).");
  endif
  ## Bound y in range [0, 1]
  y(y < 0) = 0;
  y(y > 1) = 1;
endfunction


## function for computing a trivariate normal cdf
function p = tvncdf (x, rho, tol)
  ## Get size of data
  n = size(x,1);
  ## Check if data is single or double class
  is_type = "double";
  if (isa (x, "single") || isa (rho, "single"))
    is_type = "single";
  endif

  ## Find a permutation that makes rho_32 == max(rho)
  [dum,imax] = max(abs(rho)); %#ok<ASGLU>
  if imax == 1 % swap 1 and 3
    rho_21 = rho(3); rho_31 = rho(2); rho_32 = rho(1);
    x = x(:,[3 2 1]);
  elseif imax == 2 % swap 1 and 2
    rho_21 = rho(1); rho_31 = rho(3); rho_32 = rho(2);
    x = x(:,[2 1 3]);
  else % imax == 3
    rho_21 = rho(1); rho_31 = rho(2); rho_32 = rho(3);
  end

  phi = 0.5 * erfc (- x(:,1) / sqrt (2));
  p1 = phi .* bvncdf (x(:,2:3), [], rho_32);

  if abs(rho_21) > 0
    loLimit = 0;
    hiLimit = asin(rho_21);
    rho_j1 = rho_21;
    rho_k1 = rho_31;
    p2 = zeros (size (p1), is_type);
    for i = 1:n
      b1 = x(i,1); bj = x(i,2); bk = x(i,3);
      if isfinite(b1) && isfinite(bj) && ~isnan(bk)
        p2(i) = quadgk(@tvnIntegrand,loLimit,hiLimit,'AbsTol',tol/3,'RelTol',0);
      endif
    endfor
  else
    p2 = zeros (size (p1), is_type);
  endif

  if abs(rho_31) > 0
    loLimit = 0;
    hiLimit = asin(rho_31);
    rho_j1 = rho_31;
    rho_k1 = rho_21;
    p3 = zeros (size (p1), is_type);
    for i = 1:n
      b1 = x(i,1); bj = x(i,3); bk = x(i,2);
      if isfinite(b1) && isfinite(bj) && ~isnan(bk)
        p3(i) = quadgk(@tvnIntegrand,loLimit,hiLimit,'AbsTol',tol/3,'RelTol',0);
      endif
    endfor
  else
      p3 = zeros (size (p1), is_type);
  endif

  p = cast (p1 + (p2 + p3) ./ (2 .* pi), is_type);

  function integrand = tvnIntegrand(theta)
    # Integrand is exp( -(b1.^2 + bj.^2 - 2*b1*bj*sin(theta))/(2*cos(theta).^2))
    sintheta = sin (theta);
    cossqtheta = cos (theta) .^ 2;
    expon = ((b1 * sintheta - bj) .^ 2 ./ cossqtheta + b1 .^ 2) / 2;

    sinphi = sintheta .* rho_k1 ./ rho_j1;
    numeru = bk .* cossqtheta - b1 .* (sinphi - rho_32 .* sintheta) ...
                              - bj .* (rho_32 - sintheta .* sinphi);
    denomu = sqrt (cossqtheta .* (cossqtheta - sinphi .* sinphi ...
                             - rho_32 .* (rho_32 - 2 .* sintheta .* sinphi)));
    phi = 0.5 * erfc (- (numeru ./ denomu) / sqrt (2));
    integrand = exp (- expon) .* phi;
  endfunction
endfunction

%!demo
%! mu = [1, -1];
%! Sigma = [0.9, 0.4; 0.4, 0.3];
%! [X1, X2] = meshgrid (linspace (-1, 3, 25)', linspace (-3, 1, 25)');
%! X = [X1(:), X2(:)];
%! p = mvncdf (X, mu, Sigma);
%! Z = reshape (p, 25, 25);
%! surf (X1, X2, Z);
%! title ("Bivariate Normal Distribution");
%! ylabel "X1"
%! xlabel "X2"

%!demo
%! mu = [0, 0];
%! Sigma = [0.25, 0.3; 0.3, 1];
%! p = mvncdf ([0 0], [1 1], mu, Sigma);
%! x1 = -3:.2:3;
%! x2 = -3:.2:3;
%! [X1, X2] = meshgrid (x1, x2);
%! X = [X1(:), X2(:)];
%! y = mvnpdf (X, mu, Sigma);
%! y = reshape (y, length (x2), length (x1));
%! contour (x1, x2, y, [0.0001, 0.001, 0.01, 0.05, 0.15, 0.25, 0.35]);
%! xlabel ("x");
%! ylabel ("y");
%! title ("Probability over Rectangular Region");
%! line ([0, 0, 1, 1, 0], [1, 0, 0, 1, 1], "Linestyle", "--", "Color", "k");

%!test
%! fD = (-2:2)';
%! X = repmat (fD, 1, 4);
%! p = mvncdf (X);
%! assert (p, [0; 0.0006; 0.0625; 0.5011; 0.9121], ones (5, 1) * 1e-4);
%!test
%! mu = [1, -1];
%! Sigma = [0.9, 0.4; 0.4, 0.3];
%! [X1,X2] = meshgrid (linspace (-1, 3, 25)', linspace (-3, 1, 25)');
%! X = [X1(:), X2(:)];
%! p = mvncdf (X, mu, Sigma);
%! p_out = [0.00011878988774500, 0.00034404112322371, ...
%!          0.00087682502191813, 0.00195221905058185, ...
%!          0.00378235566873474, 0.00638175749734415, ...
%!          0.00943764224329656, 0.01239164888125426, ...
%!          0.01472750274376648, 0.01623228313374828]';
%! assert (p([1:10]), p_out, 1e-16);
%!test
%! mu = [1, -1];
%! Sigma = [0.9, 0.4; 0.4, 0.3];
%! [X1,X2] = meshgrid (linspace (-1, 3, 25)', linspace (-3, 1, 25)');
%! X = [X1(:), X2(:)];
%! p = mvncdf (X, mu, Sigma);
%! p_out = [0.8180695783608276, 0.8854485749482751, ...
%!          0.9308108777385832, 0.9579855743025508, ...
%!          0.9722897881414742, 0.9788150170059926, ...
%!          0.9813597788804785, 0.9821977956568989, ...
%!          0.9824283794464095, 0.9824809345614861]';
%! assert (p([616:625]), p_out, 1e-16);
%!test
%! mu = [0, 0];
%! Sigma = [0.25, 0.3; 0.3, 1];
%! [p, err] = mvncdf ([0 0], [1 1], mu, Sigma);
%! assert (p, 0.2097424404755626, 1e-16);
%! assert (err, 1e-08);
%!test
%! x = [1 2];
%! mu = [0.5 1.5];
%! sigma = [1.0 0.5; 0.5 1.0];
%! p = mvncdf (x, mu, sigma);
%! assert (p, 0.546244443857090, 1e-15);
%!test
%! x = [1 2];
%! mu = [0.5 1.5];
%! sigma = [1.0 0.5; 0.5 1.0];
%! a = [-inf 0];
%! p = mvncdf (a, x, mu, sigma);
%! assert (p, 0.482672935215631, 1e-15);
%!error p = mvncdf (randn (25,26), [], eye (26));
%!error p = mvncdf (randn (25,8), [], eye (9));
%!error p = mvncdf (randn (25,4), randn (25,5), [], eye (4));
%!error p = mvncdf (randn (25,4), randn (25,4), [2, 3; 2, 3], eye (4));
%!error p = mvncdf (randn (25,4), randn (25,4), ones (1, 5), eye (4));
%!error p = mvncdf ([-inf 0], [1, 2], [0.5 1.5], [1.0 0.5; 0.5 1.0], option);
