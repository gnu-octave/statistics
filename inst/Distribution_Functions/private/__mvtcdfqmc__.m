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
## @deftypefn  {statistics} {@var{p} =} __mvtcdfqmc__ (@var{fname}, @var{A}, @var{B}, @var{Rho}, @var{df})
## @deftypefnx {statistics} {@var{p} =} __mvtcdfqmc__ (@dots{}, @var{TolFun})
## @deftypefnx {statistics} {@var{p} =} __mvtcdfqmc__ (@dots{}, @var{TolFun}, @var{MaxFunEvals})
## @deftypefnx {statistics} {@var{p} =} __mvtcdfqmc__ (@dots{}, @var{TolFun}, @var{MaxFunEvals}, @var{Display})
## @deftypefnx {statistics} {[@var{p}, @var{err}] =} __mvtcdfqmc__ (@dots{})
## @deftypefnx {statistics} {[@var{p}, @var{err}, @var{FunEvals}] =} __mvtcdfqmc__ (@dots{})
##
## Quasi-Monte-Carlo computation of the multivariate Student's T CDF.
##
## The QMC multivariate Student's t distribution is evaluated between the lower
## limit @var{A} and upper limit @var{B} of the hyper-rectangle with a
## correlation matrix @var{Rho} and degrees of freedom @var{df}.
##
## @var{fname} is the name of the calling function, which prefixes every error
## and warning this helper emits so that a user sees the function they called
## rather than this private one.
##
## @multitable @columnfractions 0.2 0.8
## @item "TolFun" @tab --- Maximum absolute error tolerance.  Default is 1e-4.
## @item "MaxFunEvals" @tab --- Maximum number of integrand evaluations.
## Default is 1e7 for D > 4.
## @item "Display" @tab --- Display options.  Choices are "off" (default),
## "iter", which shows the probability and estimated error at each repetition,
## and "final", which shows the final probability and related error after the
## integrand has converged successfully.
## @end multitable
##
## @code{[@var{p}, @var{err}, @var{FunEvals}] = __mvtcdfqmc__ (@dots{})} returns
## the estimated probability, @var{p}, an estimate of the error, @var{err}, and
## the number of iterations until a successful convergence is met, unless the
## value in @var{MaxFunEvals} was reached.
##
## @seealso{mvtcdf, mvtpdf, mvtrnd}
## @end deftypefn

function [p, err, FunEvals] = __mvtcdfqmc__ (fname, A, B, Rho, df, varargin)

  ## Check for input arguments
  narginchk (5,8);
  ## Add defaults
  TolFun = 1e-4;
  MaxFunEvals = 1e7;
  Display = 'off';
  ## Parse optional arguments (TolFun, MaxFunEvals, Display)
  if (nargin > 5)
    TolFun = varargin{1};
    if (! isscalar (TolFun) || ! isreal (TolFun))
      error ("%s: TolFun must be a scalar.", fname);
    endif
  endif
  if (nargin > 6)
    MaxFunEvals = varargin{2};
    if (! isscalar (MaxFunEvals) || ! isreal (MaxFunEvals))
      error ("%s: MaxFunEvals must be a scalar.", fname);
    endif
    MaxFunEvals = floor (MaxFunEvals);
  endif
  if (nargin > 7)
    Display = varargin{3};
    DispOptions = {'off', 'final', 'iter'};
    if (sum (any (strcmpi (Display, DispOptions))) == 0)
      error ("%s: invalid value for 'Display' argument.", fname);
    endif
  endif
  ## Check if input is single or double class
  is_type = 'double';
  if (isa (A, 'single') || isa (B, 'single') || isa (Rho, 'single'))
    is_type = 'single';
  endif
  ## Check for appropriate lower upper limits and NaN values in data
  if (! all (A < B))
    if (any (A > B))
      error ("%s: inconsistent lower upper limits.", fname);
    elseif (any (isnan (A) | isnan (B)))
      warning ("%s: NaNs in data.", fname);
      p = NaN (is_type);
      err = NaN (is_type);
    else
      warning ("%s: zero distance between lower upper limits.", fname);
      p = zeros (is_type);
      err = zeros (is_type);
    endif
    FunEvals = 0;
    return;
  endif
  ## Ignore dimensions with infinite limits
  InfLim_idx = (A == -Inf) & (B == Inf);
  if (any (InfLim_idx))
    if (all (InfLim_idx))
      warning ("%s: infinite distance between lower upper limits.", fname);
      p = 1;
      err = 0;
      FunEvals = 0;
      return
    endif
    A(InfLim_idx) = [];
    B(InfLim_idx) = [];
    Rho(:,InfLim_idx) = [];
    Rho(InfLim_idx,:) = [];
  endif
  ## Get size of covariance matrix
  m = size (Rho, 1);
  ## Sort the order of integration according to increasing length of interval
  [~, ord] = sort (B - A);
  A = A(ord);
  B = B(ord);
  Rho = Rho(ord, ord);
  ## Check for highly correlated covariance matrix
  if any (any (abs (tril (Rho,-1)) > .999))
    warning ("%s: highly correlated covariance matrix Rho.", fname);
  endif
  ## Scale the integration limits and the Cholesky factor of Rho
  C = chol (Rho);
  c = diag (C);
  A = A(:) ./ c;
  B = B(:) ./ c;
  C = C ./ repmat (c',m,1);
  ## Set repetitions fof Monte Carlo
  MCreps = 25;
  MCdims = m - isinf (df);
  ## Set initial output
  p = zeros (is_type);
  sigsq = Inf (is_type);
  FunEvals = 0;
  err = NaN;
  ## Initialize vector
  P = [31, 47, 73, 113, 173, 263, 397, 593, 907, 1361, 2053, 3079, 4621, ...
       6947, 10427, 15641, 23473, 35221, 52837, 79259, 118891, 178349, ...
       267523, 401287, 601942, 902933, 1354471, 2031713];
  for i = 5:length (P);
    if ((FunEvals + 2*MCreps*P(i)) > MaxFunEvals)
      break;
    endif
    ## Compute the Niederreiter point set generator
    NRgen = 2 .^ ((1:MCdims) / (MCdims + 1));
    ## Compute randomized quasi-Monte Carlo estimate with P points
    [THat,sigsqTHat] = estimate_mvtqmc (MCreps, P(i), NRgen, C, df, ...
                                        A, B, is_type);
    FunEvals = FunEvals + 2 * MCreps *P(i);
    ## Recursively update the estimate and the error estimate
    p = p + (THat - p) ./ (1 + sigsqTHat ./ sigsq);
    sigsq = sigsqTHat ./ (1 + sigsqTHat ./ sigsq);
    ## Compute a conservative estimate of error
    err = 3.5 * sqrt (sigsq);
    ## Display output for every iteration
    if (strcmpi (Display, 'iter'))
      printf ("%s: Probability estimate: %0.4f ", fname, p);
      printf ("Error estimate: %0.4e Iterations: %d\n", err, FunEvals);
    endif
    if (err < TolFun)
      if (strcmpi (Display, 'final'))
        printf ("%s: Successfully converged!\n", fname);
        printf ("Final probability estimate: %0.4f ",p);
        printf ("Final error estimate: %0.4e Iterations: %d\n", err, FunEvals);
      endif
      return
    endif
  endfor
  warning ("%s: Error tolerance did NOT converge!", fname);
  printf ("Error tolerance: %0.4f Total Iterations: %d\n", TolFun, MaxFunEvals);
endfunction

## Randomized Quasi-Monte-Carlo estimate of the integral
function [THat, sigsqTHat] = estimate_mvtqmc (MCreps, P, NRgen, C, df, A, ...
                                              B, is_type)
  qq = (1:P)' * NRgen;
  THat = zeros (MCreps,1,is_type);
  for rep = 1:MCreps
    ## Generate A new random lattice of P points.  For MVT, this is in the
    ## m-dimensional unit hypercube, for MVN, in the (m-1)-dimensional unit
    ## hypercube.
    w = abs (2 * mod (qq + repmat (rand (size (NRgen), is_type), P, 1), 1) - 1);

    ## Compute the mean of the integrand over all P of the points, and all P
    ## of the antithetic points.
    THat(rep) = (F_qrsvn (A, B, C, df, w) + F_qrsvn (A, B, C, df, 1 - w)) ./ 2;
  endfor

  ## Return the MC mean and se^2
  sigsqTHat = var (THat) ./ MCreps;
  THat = mean (THat);
endfunction

## Integrand for computation of MVT probabilities
function TBar = F_qrsvn (A, B, C, df, w)
  N = size (w, 1);  # number of quasirandom points
  m = length (A);    # number of dimensions
  if isinf (df)
    rho = 1;
  else
    rho = chi_inv (w(:,m), df) ./ sqrt (df);
  endif
  rA = norm_cdf (rho .* A(1));       # A is already scaled by diag(C)
  rB = norm_cdf (rho .* B(1)) - rA;  # B is already scaled by diag(C)
  T = rB;
  Y = zeros (N, m, 'like', T);
  for i = 2:m
    z = min (max (rA + rB .* w(:,i-1), eps / 2), 1 - eps / 2);
    Y(:,i-1) = norm_inv (z);
    Ysum = Y * C(:,i);
    rA = norm_cdf (rho .* A(i) - Ysum);       # A is already scaled by diag(C)
    rB = norm_cdf (rho .* B(i) - Ysum) - rA;  # B is already scaled by diag(C)
    T = T .* rB;
  endfor
  TBar = sum (T, 1) ./ length (T);
endfunction

## Normal cumulative distribution function
function a = norm_cdf (b)
  a = 0.5 * erfc (- b ./ sqrt (2));
endfunction

## Inverse of normal cumulative distribution function
function a = norm_inv (b)
a = - sqrt (2) .* erfcinv (2 * b);
endfunction

## Inverse of chi cumulative distribution function
function a = chi_inv (b,df)
  a = sqrt (gammaincinv (b, df ./ 2) .* 2);
endfunction



