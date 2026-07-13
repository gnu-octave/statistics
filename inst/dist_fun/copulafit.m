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
## @deftypefn  {statistics} {@var{rho} =} copulafit (@qcode{"Gaussian"}, @var{u})
## @deftypefnx {statistics} {[@var{rho}, @var{nu}] =} copulafit (@qcode{"t"}, @var{u})
## @deftypefnx {statistics} {[@var{param}, @var{ci}] =} copulafit (@var{family}, @var{u})
## @deftypefnx {statistics} {[@dots{}] =} copulafit (@dots{}, @qcode{"alpha"}, @var{a})
##
## Fit a copula to data.
##
## @code{copulafit (@var{family}, @var{u})} returns the maximum-likelihood
## estimate of the parameter of a copula of the family @var{family}, fit to the
## data in @var{u}.  The rows of @var{u} are observations and its columns are
## variables; all entries must lie strictly inside the unit interval
## @math{(0,1)}, as produced for example by a probability-integral transform or
## by @code{ecdf}/@code{ksdensity}.
##
## @var{family} is the copula family name.  It can be @qcode{"Gaussian"} for the
## Gaussian family, @qcode{"t"} for the Student's t family, @qcode{"Clayton"}
## for the Clayton family, @qcode{"Gumbel"} for the Gumbel-Hougaard family, or
## @qcode{"Frank"} for the Frank family.
##
## The returned value depends on the family:
##
## @itemize @bullet
## @item For @qcode{"Gaussian"}, @code{@var{rho} = copulafit ("Gaussian",
## @var{u})} returns the estimated linear correlation matrix @var{rho}, computed
## as the sample correlation of the normal scores @code{norminv (@var{u})}.  The
## data may have two or more columns.
##
## @item For @qcode{"t"}, @code{copulafit ("t", @var{u})} returns the estimated
## correlation matrix @var{rho} and the degrees of freedom @var{nu} as
## @code{[@var{rho}, @var{nu}]}, obtained by maximizing the copula
## log-likelihood.  Only bivariate data (two columns) are supported.
##
## @item For the Archimedean families @qcode{"Clayton"}, @qcode{"Gumbel"}, and
## @qcode{"Frank"}, @code{[@var{param}, @var{ci}] = copulafit (@var{family},
## @var{u})} returns the scalar copula parameter @var{param} and, optionally, a
## two-element vector @var{ci} with the lower and upper confidence bounds.  Only
## bivariate data are supported.
## @end itemize
##
## @code{copulafit (@dots{}, @qcode{"alpha"}, @var{a})} sets the significance
## level for the confidence interval to @var{a}, so that @var{ci} has coverage
## @code{100 * (1 - @var{a})} percent.  The default is @code{@var{a} = 0.05}.
## The confidence interval is a Wald interval whose standard error is obtained
## from the outer-product-of-gradients estimate of the information.
##
## @seealso{copulastat, copulaparam, copulacdf, copulapdf, copularnd}
## @end deftypefn

function varargout = copulafit (family, u, varargin)

  ## Check arguments
  if (nargin < 2)
    print_usage ();
  endif

  if (! ischar (family))
    error (strcat ("copulafit: FAMILY must be one of 'Gaussian',", ...
                   " 't', 'Clayton', 'Gumbel', and 'Frank'."));
  endif

  if (! isnumeric (u) || ! isreal (u) || ! ismatrix (u) || isempty (u))
    error ("copulafit: U must be a nonempty numeric matrix.");
  endif

  if (any (u(:) <= 0) || any (u(:) >= 1))
    error ("copulafit: U must have all values in the open interval (0, 1).");
  endif

  [n, d] = size (u);
  if (d < 2)
    error ("copulafit: U must have at least two columns.");
  endif

  ## Parse the 'alpha' option
  alpha = 0.05;
  if (numel (varargin) > 0)
    if (numel (varargin) != 2 || ! ischar (varargin{1}) || ...
        ! strcmpi (varargin{1}, 'alpha'))
      error ("copulafit: invalid optional argument.");
    endif
    alpha = varargin{2};
    if (! (isnumeric (alpha) && isscalar (alpha) && isreal (alpha) ...
           && alpha > 0 && alpha < 1))
      error ("copulafit: ALPHA must be a scalar in the range (0, 1).");
    endif
  endif

  lower_family = lower (family);

  switch (lower_family)

    case 'gaussian'
      if (nargout > 1)
        error (strcat ("copulafit: the Gaussian copula fit returns only", ...
                       " the correlation matrix."));
      endif
      varargout{1} = correlation_of_scores (norminv (u));

    case 't'
      if (d != 2)
        error (strcat ("copulafit: the t copula fit supports bivariate", ...
                       " data only."));
      endif
      if (nargout > 2)
        error (strcat ("copulafit: confidence intervals are not supported", ...
                       " for the t copula fit."));
      endif
      [rho, nu] = fit_t (u);
      varargout{1} = [1, rho; rho, 1];
      varargout{2} = nu;

    case {'clayton', 'gumbel', 'frank'}
      if (d != 2)
        error (strcat ("copulafit: the %s copula fit supports bivariate", ...
                       " data only."), family);
      endif
      [param, ci] = fit_archimedean (lower_family, u, alpha);
      varargout{1} = param;
      if (nargout > 1)
        varargout{2} = ci;
      endif

    otherwise
      error ("copulafit: unknown copula family '%s'.", family);

  endswitch

endfunction

## Sample correlation matrix of the columns of z (mean-removed, unit-normalised)
function R = correlation_of_scores (z)
  R = corr (z);
  ## Guard the diagonal against round-off so it is exactly one
  R(1 : (rows (R) + 1) : end) = 1;
endfunction

## Maximum-likelihood fit of the bivariate Student's t copula.  Optimises the
## correlation and the degrees of freedom jointly on an unconstrained scale
## (rho = tanh (a), nu = exp (b)).
function [rho, nu] = fit_t (u)
  rho0 = corr (norminv (u))(1, 2);
  p0 = [atanh(rho0); log(4)];
  opts = optimset ("TolX", 1e-8, "TolFun", 1e-8, ...
                   "MaxFunEvals", 20000, "MaxIter", 20000);
  p = fminsearch (@(p) t_copula_nll (u, tanh (p(1)), exp (p(2))), p0, opts);
  rho = tanh (p(1));
  nu = exp (p(2));
endfunction

## Negative log-likelihood of the bivariate t copula
function v = t_copula_nll (u, rho, nu)
  if (abs (rho) >= 1 || nu <= 0)
    v = Inf;
    return;
  endif
  R = [1, rho; rho, 1];
  t = tinv (u, nu);
  c = mvtpdf (t, R, nu) ./ (tpdf (t(:, 1), nu) .* tpdf (t(:, 2), nu));
  if (any (! isfinite (c)) || any (c <= 0))
    v = Inf;
    return;
  endif
  v = -sum (log (c));
endfunction

## Maximum-likelihood fit of a bivariate Archimedean copula with a Wald
## confidence interval whose standard error uses the outer-product-of-gradients
## (BHHH) estimate of the information.
function [param, ci] = fit_archimedean (family, u, alpha)
  a0 = archimedean_start (family);
  opts = optimset ("TolX", 1e-10, "TolFun", 1e-10, ...
                   "MaxFunEvals", 10000, "MaxIter", 10000);
  param = fminsearch (@(a) -sum (log (copulapdf (family, u, a))), a0, opts);
  ## Per-observation score by central difference, then OPG information
  h = max (1e-6, abs (param) .* 1e-5);
  logc = @(a) log (copulapdf (family, u, a));
  score = (logc (param + h) - logc (param - h)) ./ (2 .* h);
  se = 1 ./ sqrt (sum (score .^ 2));
  z = norminv (1 - alpha ./ 2);
  ci = [param - z .* se, param + z .* se];
endfunction

## A robust starting value for the Archimedean maximum-likelihood search
function a0 = archimedean_start (family)
  switch (family)
    case 'clayton'
      a0 = 1;
    case 'gumbel'
      a0 = 1.5;
    case 'frank'
      a0 = 1;
  endswitch
endfunction

%!demo
%! ## Fit a Clayton copula to data and recover a confidence interval
%! u = copularnd ("Clayton", 2, 500);
%! [alpha, ci] = copulafit ("Clayton", u)

%!demo
%! ## Fit a Gaussian copula and report the correlation matrix
%! u = copularnd ("Gaussian", 0.6, 500);
%! rho = copulafit ("Gaussian", u)

## Shared probe data (near-Gaussian bivariate sample), reference values from
## MATLAB's copulafit.
%!shared u
%! u = [0.08,0.12; 0.17,0.25; 0.23,0.19; 0.31,0.42; 0.39,0.35; 0.46,0.51; ...
%!      0.52,0.48; 0.58,0.63; 0.64,0.59; 0.71,0.68; 0.77,0.82; 0.83,0.79; ...
%!      0.88,0.91; 0.93,0.87; 0.97,0.95];

## Gaussian: sample correlation of the normal scores (exact match to MATLAB)
%!test
%! rho = copulafit ("Gaussian", u);
%! assert_equal (rho, [1, 0.979591430658725; 0.979591430658725, 1], 1e-12);

## Archimedean point estimates and Wald confidence intervals (match MATLAB)
%!test
%! [a, ci] = copulafit ("Clayton", u);
%! assert_equal (a, 8.70842970823662, 1e-6);
%! assert_equal (ci, [4.48832722027934, 12.9285321961939], 1e-5);
%!test
%! [a, ci] = copulafit ("Frank", u);
%! assert_equal (a, 29.7092786222299, 1e-5);
%! assert_equal (ci, [7.05441865208815, 52.3641385923717], 1e-4);
%!test
%! [a, ci] = copulafit ("Gumbel", u);
%! assert_equal (a, 6.14597585483989, 1e-6);
%! assert_equal (ci, [2.94928782240414, 9.34266388727563], 1e-5);

## The t copula fit recovers the correlation (its nu is weakly identified for
## near-Gaussian data, so only rho is checked tightly).
%!test
%! [rho, nu] = copulafit ("t", u);
%! assert_equal (rho, [1, 0.98084; 0.98084, 1], 2e-3);
%! assert (nu > 10);

## The confidence level widens the interval as alpha shrinks
%!test
%! [~, ci95] = copulafit ("Clayton", u);
%! [~, ci99] = copulafit ("Clayton", u, "alpha", 0.01);
%! assert (ci99(1) < ci95(1) && ci99(2) > ci95(2));

## Round trip: fit recovers a parameter close to the generating one
%!test
%! rng (42);
%! u = copularnd ("Clayton", 2, 2000);
%! a = copulafit ("Clayton", u);
%! assert_equal (a, 2, 0.3);

## Test input validation
%!error <copulafit: FAMILY must be one of 'Gaussian', 't', 'Clayton', 'Gumbel', and 'Frank'.> ...
%! copulafit (5, [0.2, 0.3])
%!error <copulafit: U must be a nonempty numeric matrix.> ...
%! copulafit ("Gaussian", "foo")
%!error <copulafit: U must have all values in the open interval > ...
%! copulafit ("Gaussian", [0.2, 0.3; 1.2, 0.5])
%!error <copulafit: U must have at least two columns.> ...
%! copulafit ("Gaussian", [0.2; 0.3; 0.4])
%!error <copulafit: the Gaussian copula fit returns only the correlation matrix.> ...
%! [a, b] = copulafit ("Gaussian", [0.2, 0.3; 0.4, 0.5]);
%!error <copulafit: the t copula fit supports bivariate data only.> ...
%! copulafit ("t", [0.2, 0.3, 0.4; 0.5, 0.6, 0.7])
%!error <copulafit: the Clayton copula fit supports bivariate data only.> ...
%! copulafit ("Clayton", [0.2, 0.3, 0.4; 0.5, 0.6, 0.7])
%!error <copulafit: invalid optional argument.> ...
%! copulafit ("Clayton", [0.2, 0.3; 0.4, 0.5], "foo")
%!error <copulafit: ALPHA must be a scalar in the range > ...
%! copulafit ("Clayton", [0.2, 0.3; 0.4, 0.5], "alpha", 2)
%!error <copulafit: unknown copula family 'Foo'.> ...
%! copulafit ("Foo", [0.2, 0.3; 0.4, 0.5])
