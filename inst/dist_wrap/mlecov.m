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
## @deftypefn  {statistics} {@var{acov} =} mlecov (@var{params}, @var{data}, @var{Name}, @var{Value})
##
## Asymptotic covariance matrix of maximum likelihood estimators.
##
## @code{@var{acov} = mlecov (@var{params}, @var{data}, @dots{})}
## returns an approximation to the asymptotic covariance matrix of the maximum
## likelihood estimators of the parameters of a distribution, evaluated at the
## parameter values in @var{params} for the sample data in @var{data}.
## @var{params} is a numeric vector of parameter values (typically the estimates
## returned by @code{mle} or @code{fitdist}) and @var{data} is a numeric vector
## of the sample observations.  @var{acov} is a @math{p*p} matrix, where
## @math{p = numel (@var{params})}.
##
## The distribution is not identified by name; instead it is supplied through
## @qcode{Name-Value} paired arguments that give function handles to its
## density, its log density, or its negative log-likelihood.  Exactly
## @strong{one} of the following three arguments must be specified:
##
## @multitable @columnfractions 0.18 0.8
## @headitem @var{Name} @tab @var{Value}
##
## @item @qcode{'pdf'} @tab A function handle,
## @code{@var{f}(@var{data}, @var{p1}, @var{p2}, @dots{})}, that accepts the
## sample data as its first argument and the distribution parameters as
## subsequent scalar arguments, and returns a vector of probability density
## values, one per observation.
##
## @item @qcode{'logpdf'} @tab A function handle,
## @code{@var{f}(@var{data}, @var{p1}, @var{p2}, @dots{})}, with the same
## calling convention as @qcode{'pdf'} but returning the @emph{logarithm} of
## the density.
##
## @item @qcode{'nloglf'} @tab A function handle,
## @code{@var{nll}(@var{params}, @var{data}, @var{cens}, @var{freq})}, that
## returns the scalar negative log-likelihood of the whole sample.  It receives
## the current parameter vector, the data, the censoring vector, and the
## frequency vector, and is responsible for incorporating censoring and
## frequency itself.
##
## @item @qcode{'cdf'} @tab A function handle to the cumulative distribution
## function, with the same calling convention as @qcode{'pdf'}.  It is
## @strong{required} together with @qcode{'pdf'} when the data are censored, so
## that censored observations can contribute their survival probability.
##
## @item @qcode{'logsf'} @tab A function handle to the logarithm of the survivor
## function @math{log (1 - cdf)}, with the same calling convention as
## @qcode{'pdf'}.  It is @strong{required} together with @qcode{'logpdf'} when
## the data are censored.
##
## @item @qcode{'Censoring'} @tab A vector of the same size as @var{data}
## indicating censored observations (nonzero for right-censored).  By default no
## observation is censored.
##
## @item @qcode{'Frequency'} @tab A vector of nonnegative integer counts of the
## same size as @var{data}, giving the number of times each observation was
## observed.  By default it is @qcode{ones (size (@var{data}))}.
##
## @item @qcode{'Options'} @tab A structure that may contain a
## @qcode{'DerivStep'} field specifying the relative finite-difference step
## used to approximate the Hessian (a positive scalar or a vector the same size
## as @var{params}).  The default step is @qcode{eps ^ (1/4)}.
## @end multitable
##
## @strong{Computation and numerical behavior.}  @code{mlecov} approximates the
## covariance matrix as the inverse of the observed Fisher information, that is,
## the inverse of the Hessian of the @emph{aggregate} negative log-likelihood of
## the sample, evaluated by central finite differences at @var{params}.  The
## covariance is computed @emph{at} the supplied @var{params}; @code{mlecov}
## does not refit the parameters, so @var{params} should be the maximum
## likelihood estimates for the result to be meaningful.
##
## Whichever of @qcode{'pdf'}, @qcode{'logpdf'}, or @qcode{'nloglf'} is
## supplied, the Hessian is always formed by differencing the same aggregate
## negative log-likelihood rather than by differentiating the density itself.
## This makes the three input forms consistent with one another and is
## numerically far more stable than differentiating a density; as a consequence
## @var{acov} may differ from other implementations (including MATLAB) in
## ill-conditioned cases where those differentiate the density directly and
## return unreliable values or @code{NaN}.  If the computed Hessian is not
## positive definite (for example when @var{params} is not at a likelihood
## maximum), a warning is issued and @var{acov} is returned as an
## all-@code{NaN} matrix.
##
## @seealso{mle, fitdist, makedist}
## @end deftypefn

function acov = mlecov (params, data, varargin)

  ## Check number of input arguments
  if (nargin < 3)
    print_usage ();
  endif

  ## Check PARAMS and DATA
  if (! (isvector (params) && isnumeric (params) && isreal (params)
                           && ! isempty (params)))
    error ("mlecov: PARAMS must be a nonempty numeric vector of real values.");
  endif
  if (! (isvector (data) && isnumeric (data) && isreal (data)
                        && ! isempty (data)))
    error ("mlecov: DATA must be a nonempty numeric vector of real values.");
  endif

  ## Add defaults
  pdf = [];
  logpdf = [];
  nloglf = [];
  cdf = [];
  logsf = [];
  censor = [];
  freq = [];
  derivstep = eps ^ (1/4);

  ## Parse optional arguments as NAME-VALUE pairs
  if (mod (numel (varargin), 2) != 0)
    error ("mlecov: optional arguments must be in NAME-VALUE pairs.");
  endif
  while (numel (varargin) > 0)
    name = varargin{1};
    value = varargin{2};
    if (! (ischar (name) && isrow (name)))
      error ("mlecov: NAME arguments must be character vectors.");
    endif
    switch (tolower (name))
      case 'pdf'
        pdf = value;
      case 'logpdf'
        logpdf = value;
      case 'nloglf'
        nloglf = value;
      case 'cdf'
        cdf = value;
      case 'logsf'
        logsf = value;
      case 'censoring'
        censor = value;
      case 'frequency'
        freq = value;
      case 'options'
        if (! isstruct (value))
          error ("mlecov: 'Options' argument must be a structure.");
        endif
        if (isfield (value, 'DerivStep') && ! isempty (value.DerivStep))
          derivstep = value.DerivStep;
        endif
      otherwise
        error ("mlecov: unknown parameter name '%s'.", name);
    endswitch
    varargin([1:2]) = [];
  endwhile

  ## Exactly one distribution function must be specified
  nfun = (! isempty (pdf)) + (! isempty (logpdf)) + (! isempty (nloglf));
  if (nfun == 0)
    error (strcat ("mlecov: a distribution must be specified with one of", ...
                   " the 'pdf', 'logpdf', or 'nloglf' arguments."));
  elseif (nfun > 1)
    error (strcat ("mlecov: only one of the 'pdf', 'logpdf', or 'nloglf'", ...
                   " arguments can be specified."));
  endif

  ## Check that supplied handles are function handles
  if (! isempty (pdf) && ! is_function_handle (pdf))
    error ("mlecov: 'pdf' argument must be a function handle.");
  endif
  if (! isempty (logpdf) && ! is_function_handle (logpdf))
    error ("mlecov: 'logpdf' argument must be a function handle.");
  endif
  if (! isempty (nloglf) && ! is_function_handle (nloglf))
    error ("mlecov: 'nloglf' argument must be a function handle.");
  endif
  if (! isempty (cdf) && ! is_function_handle (cdf))
    error ("mlecov: 'cdf' argument must be a function handle.");
  endif
  if (! isempty (logsf) && ! is_function_handle (logsf))
    error ("mlecov: 'logsf' argument must be a function handle.");
  endif

  ## Add defaults and validate FREQUENCY and CENSORING
  if (isempty (freq))
    freq = ones (size (data));
  elseif (! isequal (size (data), size (freq)))
    error (strcat ("mlecov: 'Frequency' argument must have the same size", ...
                   " as the input data in DATA."));
  elseif (any (freq(:) < 0) || any (freq(:) != round (freq(:))))
    error (strcat ("mlecov: 'Frequency' argument must contain non-negative", ...
                   " integer values."));
  endif
  if (isempty (censor))
    censor = zeros (size (data));
  elseif (! isequal (size (data), size (censor)))
    error (strcat ("mlecov: 'Censoring' argument must have the same size", ...
                   " as the input data in DATA."));
  endif
  docens = any (censor(:) != 0);

  ## Censored data need a survivor function for the 'pdf'/'logpdf' forms
  if (docens && ! isempty (pdf) && isempty (cdf))
    error (strcat ("mlecov: a 'cdf' function handle is required for", ...
                   " censored data when using the 'pdf' argument."));
  endif
  if (docens && ! isempty (logpdf) && isempty (logsf))
    error (strcat ("mlecov: a 'logsf' function handle is required for", ...
                   " censored data when using the 'logpdf' argument."));
  endif

  ## Validate DerivStep and build the per-parameter finite-difference step
  theta = params(:).';
  p = numel (theta);
  if (! (isnumeric (derivstep) && isreal (derivstep) && all (derivstep(:) > 0)
                               && (isscalar (derivstep)
                                   || numel (derivstep) == p)))
    error (strcat ("mlecov: 'DerivStep' must be a positive real scalar or", ...
                   " a vector the same size as PARAMS."));
  endif
  hstep = derivstep(:).' .* max (abs (theta), 1);

  ## Assemble the aggregate negative log-likelihood as a function of the
  ## parameter vector only, capturing data, censoring, and frequency
  if (! isempty (pdf))
    form = 'pdf';
  elseif (! isempty (logpdf))
    form = 'logpdf';
  else
    form = 'nloglf';
  endif
  nllfun = @(t) aggregate_nll (t, form, pdf, cdf, logpdf, logsf, nloglf, ...
                               data, censor, freq, docens);

  ## Central finite-difference Hessian of the negative log-likelihood
  H = num_hessian (nllfun, theta, hstep);
  H = (H + H') / 2;

  ## Invert only if the Hessian is positive definite; otherwise warn and NaN
  [~, notpd] = chol (H);
  if (notpd != 0)
    warning (strcat ("mlecov: unable to compute a covariance matrix", ...
                     " because the computed Hessian matrix is not positive", ...
                     " definite."));
    acov = NaN (p);
  else
    acov = inv (H);
    acov = (acov + acov') / 2;
  endif

endfunction

## Aggregate negative log-likelihood at parameter vector T
function nll = aggregate_nll (t, form, pdf, cdf, logpdf, logsf, nloglf, ...
                              data, censor, freq, docens)
  switch (form)
    case 'pdf'
      pc = num2cell (t);
      dens = pdf (data, pc{:});
      if (docens)
        surv = 1 - cdf (data, pc{:});
        terms = (censor == 0) .* log (dens) + (censor != 0) .* log (surv);
      else
        terms = log (dens);
      endif
      nll = -sum (freq(:) .* terms(:));
    case 'logpdf'
      pc = num2cell (t);
      lpd = logpdf (data, pc{:});
      if (docens)
        lsf = logsf (data, pc{:});
        terms = (censor == 0) .* lpd + (censor != 0) .* lsf;
      else
        terms = lpd;
      endif
      nll = -sum (freq(:) .* terms(:));
    case 'nloglf'
      nll = nloglf (t, data, censor, freq);
  endswitch
endfunction

## Central finite-difference approximation of the Hessian of NLLFUN at THETA
function H = num_hessian (nllfun, theta, hstep)
  p = numel (theta);
  H = zeros (p);
  f0 = nllfun (theta);
  for i = 1:p
    ei = zeros (1, p);
    ei(i) = hstep(i);
    fpi = nllfun (theta + ei);
    fmi = nllfun (theta - ei);
    H(i,i) = (fpi - 2 * f0 + fmi) / (hstep(i) ^ 2);
    for j = (i + 1):p
      ej = zeros (1, p);
      ej(j) = hstep(j);
      fpp = nllfun (theta + ei + ej);
      fpm = nllfun (theta + ei - ej);
      fmp = nllfun (theta - ei + ej);
      fmm = nllfun (theta - ei - ej);
      H(i,j) = (fpp - fpm - fmp + fmm) / (4 * hstep(i) * hstep(j));
      H(j,i) = H(i,j);
    endfor
  endfor
endfunction

%!demo
%! ## Asymptotic covariance matrix of the ML estimates for a normal fit.
%! x = [2.1, 3.4, 1.9, 5.2, 4.1, 2.8, 3.3, 4.7, 2.2, 3.9, 3.0, 4.5];
%! phat = mle (x);
%! acov = mlecov (phat, x, 'pdf', @(x, mu, sigma) normpdf (x, mu, sigma))

## Reference values below were computed with MATLAB's mlecov using the reliable
## 'nloglf' input form (its 'pdf'/'logpdf' paths differentiate the density
## directly and return NaN or unreliable values on these inputs).

%!test
%! x = [2.1, 3.4, 1.9, 5.2, 4.1, 2.8, 3.3, 4.7, 2.2, 3.9, 3.0, 4.5];
%! phat = [mean(x), std(x, 1)];
%! nll = @(p, data, cens, freq) -sum (log (normpdf (data, p(1), p(2))));
%! acov = mlecov (phat, x, 'nloglf', nll);
%! assert (acov, [0.0887673606073251, 0; 0, 0.0443836789388774], 1e-6);

%!test
%! x = [2.1, 3.4, 1.9, 5.2, 4.1, 2.8, 3.3, 4.7, 2.2, 3.9, 3.0, 4.5];
%! phat = [3.5, 1.0];
%! nll = @(p, data, cens, freq) -sum (log (normpdf (data, p(1), p(2))));
%! acov = mlecov (phat, x, 'nloglf', nll);
%! ref = [0.0841894975321553, 0.00570776282157731; ...
%!        0.00570776282157731, 0.0380517511049636];
%! assert (acov, ref, 1e-6);

%!test
%! x = [2.1, 3.4, 1.9, 5.2, 4.1, 2.8, 3.3, 4.7, 2.2, 3.9, 3.0, 4.5];
%! phat = gamfit (x);
%! nll = @(p, data, cens, freq) -sum (log (gampdf (data, p(1), p(2))));
%! acov = mlecov (phat, x, 'nloglf', nll);
%! ref = [17.6711714212941, -0.553235139352906; ...
%!        -0.553235139352906, 0.0181745610037496];
%! assert (acov, ref, 5e-4);

%!test
%! x = [2.1, 3.4, 1.9, 5.2, 4.1, 2.8, 3.3, 4.7, 2.2, 3.9, 3.0, 4.5];
%! f = [1, 2, 1, 1, 3, 1, 2, 1, 1, 1, 2, 1];
%! phat = [mean(x), std(x, 1)];
%! nll = @(p, data, cens, freq) -sum (freq .* log (normpdf (data, p(1), p(2))));
%! acov = mlecov (phat, x, 'nloglf', nll, 'Frequency', f);
%! ref = [0.0630373869314514, -0.00427967204059965; ...
%!        -0.00427967204059965, 0.0484445550668851];
%! assert (acov, ref, 1e-6);

%!test
%! x = [2.1, 3.4, 1.9, 5.2, 4.1, 2.8, 3.3, 4.7, 2.2, 3.9, 3.0, 4.5];
%! c = [0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0];
%! phat = [mean(x), std(x, 1)];
%! nll = @(p, data, cens, freq) -sum ((1 - cens) .* ...
%!         log (normpdf (data, p(1), p(2))) + ...
%!         cens .* log (1 - normcdf (data, p(1), p(2))));
%! acov = mlecov (phat, x, 'nloglf', nll, 'Censoring', c);
%! ref = [0.100361972376949, -0.0148199117431384; ...
%!        -0.0148199117431384, 0.054254671056408];
%! assert (acov, ref, 1e-6);

## The 'pdf', 'logpdf', and 'nloglf' forms differentiate the same aggregate
## negative log-likelihood, so they agree (unlike MATLAB, whose density-based
## 'pdf' path returns NaN for the normal and 51.84 for the exponential here).
%!test
%! x = [2.1, 3.4, 1.9, 5.2, 4.1, 2.8, 3.3, 4.7, 2.2, 3.9, 3.0, 4.5];
%! phat = [mean(x), std(x, 1)];
%! a_pdf = mlecov (phat, x, 'pdf', @(x, mu, s) normpdf (x, mu, s));
%! a_logpdf = mlecov (phat, x, 'logpdf', @(x, mu, s) log (normpdf (x, mu, s)));
%! a_nloglf = mlecov (phat, x, 'nloglf', ...
%!                    @(p, d, c, f) -sum (log (normpdf (d, p(1), p(2)))));
%! assert (a_pdf, a_nloglf, 1e-10);
%! assert (a_logpdf, a_nloglf, 1e-10);
%! assert (a_pdf, [0.0887673606073251, 0; 0, 0.0443836789388774], 1e-6);

%!test
%! x = [2.1, 3.4, 1.9, 5.2, 4.1, 2.8, 3.3, 4.7, 2.2, 3.9, 3.0, 4.5];
%! acov = mlecov (mean (x), x, 'pdf', @(x, mu) exppdf (x, mu));
%! assert (acov, 0.977552156166645, 1e-6);

## Censored 'pdf' form (with a 'cdf' handle) matches the hand-built censored
## 'nloglf' likelihood.
%!test
%! x = [2.1, 3.4, 1.9, 5.2, 4.1, 2.8, 3.3, 4.7, 2.2, 3.9, 3.0, 4.5];
%! c = [0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0];
%! phat = [mean(x), std(x, 1)];
%! acov = mlecov (phat, x, 'pdf', @(x, mu, s) normpdf (x, mu, s), ...
%!                'cdf', @(x, mu, s) normcdf (x, mu, s), 'Censoring', c);
%! ref = [0.100361972376949, -0.0148199117431384; ...
%!        -0.0148199117431384, 0.054254671056408];
%! assert (acov, ref, 1e-6);

## A user-supplied 'DerivStep' still yields the correct covariance.
%!test
%! x = [2.1, 3.4, 1.9, 5.2, 4.1, 2.8, 3.3, 4.7, 2.2, 3.9, 3.0, 4.5];
%! phat = [mean(x), std(x, 1)];
%! nll = @(p, data, cens, freq) -sum (log (normpdf (data, p(1), p(2))));
%! acov = mlecov (phat, x, 'nloglf', nll, ...
%!                'Options', struct ('DerivStep', 1e-4));
%! assert (acov, [0.0887673606073251, 0; 0, 0.0443836789388774], 1e-5);

## A non-positive-definite Hessian yields a warning and an all-NaN matrix.
%!warning <mlecov: unable to compute a covariance matrix because the computed Hessian matrix is not positive definite.> ...
%! mlecov (1, [1, 2, 3, 4, 5], 'nloglf', @(p, d, c, f) -sum (p(1) .* d));
%!test
%! warning ("off", "all", "local");
%! acov = mlecov (1, [1, 2, 3, 4, 5], 'nloglf', @(p, d, c, f) -sum (p(1) .* d));
%! assert (isnan (acov));

## Test input validation
%!error <Invalid call to mlecov> mlecov (1, [1, 2, 3])
%!error <mlecov: PARAMS must be a nonempty numeric vector of real values.> ...
%! mlecov ([1, 2; 3, 4], [1, 2, 3], 'pdf', @(x, a) x)
%!error <mlecov: PARAMS must be a nonempty numeric vector of real values.> ...
%! mlecov ([1, 2i], [1, 2, 3], 'pdf', @(x, a, b) x)
%!error <mlecov: DATA must be a nonempty numeric vector of real values.> ...
%! mlecov ([1, 2], ones (2, 2), 'pdf', @(x, a, b) x)
%!error <mlecov: optional arguments must be in NAME-VALUE pairs.> ...
%! mlecov ([1, 2], [1, 2, 3], 'pdf')
%!error <mlecov: NAME arguments must be character vectors.> ...
%! mlecov ([1, 2], [1, 2, 3], 5, @sin)
%!error <mlecov: a distribution must be specified with one of the 'pdf', 'logpdf', or 'nloglf' arguments.> ...
%! mlecov ([1, 2], [1, 2, 3], 'Frequency', [1, 1, 1])
%!error <mlecov: only one of the 'pdf', 'logpdf', or 'nloglf' arguments can be specified.> ...
%! mlecov ([1, 2], [1, 2, 3], 'pdf', @sin, 'nloglf', @cos)
%!error <mlecov: 'pdf' argument must be a function handle.> ...
%! mlecov ([1, 2], [1, 2, 3], 'pdf', 5)
%!error <mlecov: 'nloglf' argument must be a function handle.> ...
%! mlecov ([1, 2], [1, 2, 3], 'nloglf', 'text')
%!error <mlecov: unknown parameter name 'bogus'.> ...
%! mlecov ([1, 2], [1, 2, 3], 'pdf', @sin, 'bogus', 1)
%!error <mlecov: 'Frequency' argument must have the same size as the input data in DATA.> ...
%! mlecov ([1, 2], [1, 2, 3], 'nloglf', @(varargin) 1, 'Frequency', [1, 1])
%!error <mlecov: 'Frequency' argument must contain non-negative integer values.> ...
%! mlecov ([1, 2], [1, 2, 3], 'nloglf', @(varargin) 1, 'Frequency', [1, 0.5, 1])
%!error <mlecov: 'Censoring' argument must have the same size as the input data in DATA.> ...
%! mlecov ([1, 2], [1, 2, 3], 'nloglf', @(varargin) 1, 'Censoring', [1, 0])
%!error <mlecov: a 'cdf' function handle is required for censored data when using the 'pdf' argument.> ...
%! mlecov ([1, 2], [1, 2, 3], 'pdf', @(x, a, b) x, 'Censoring', [1, 0, 0])
%!error <mlecov: a 'logsf' function handle is required for censored data when using the 'logpdf' argument.> ...
%! mlecov ([1, 2], [1, 2, 3], 'logpdf', @(x, a, b) x, 'Censoring', [1, 0, 0])
%!error <mlecov: 'Options' argument must be a structure.> ...
%! mlecov ([1, 2], [1, 2, 3], 'nloglf', @(varargin) 1, 'Options', 5)
%!error <mlecov: 'DerivStep' must be a positive real scalar or a vector the same size as PARAMS.> ...
%! mlecov ([1, 2], [1, 2, 3], 'nloglf', @(varargin) 1, ...
%!         'Options', struct ('DerivStep', -1))
