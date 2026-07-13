## Copyright (C) 2024 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {statistics} {@var{phat} =} mle (@var{x})
## @deftypefnx {statistics} {@var{phat} =} mle (@var{x}, @var{Name}, @var{Value})
## @deftypefnx {statistics} {[@var{phat}, @var{pci}] =} mle (@dots{})
##
## Compute maximum likelihood estimates.
##
## @code{@var{phat} = mle (@var{x})} returns the maximum likelihood estimates
## (MLEs) for the parameters of a normal distribution using the sample data in
## @var{x}, which must be a numeric vector of real values.
##
## @code{@var{phat} = mle (@var{x}, @var{Name}, @var{Value})} returns the MLEs
## with additional options specified by @qcode{Name-Value} pair arguments listed
## below.
##
## @multitable @columnfractions 0.18 0.8
## @headitem @var{Name} @tab @var{Value}
##
## @item @qcode{'distribution'} @tab A character vector specifying the
## distribution type for which to estimate parameters.
##
## @item @qcode{'Ntrials'} @tab A scalar specifying the number of trials
## for the corresponding element of @var{x} for the binomial distribution.
##
## @item @qcode{'theta'} @tab A scalar specifying the location parameter
## for the generalized Pareto distribution.
##
## @item @qcode{'mu'} @tab A scalar specifying the location parameter
## for the half-normal distribution.
##
## @item @qcode{'censoring'} @tab A vector of the same size as @var{x}
## indicating censored data in @var{x}.  By default it is
## @qcode{@var{censor} = zeros (size (@var{x}))}.
##
## @item @qcode{'frequency'} @tab A vector of nonnegative integer counts of
## the same size as @var{x} used as frequency observations.  By default it is
## @qcode{@var{freq} = ones (size (@var{x}))}.
##
## @item @qcode{'alpha'} @tab A scalar in the range @math{(0,1)}, as the
## significance level for the confidence interval @var{pci}.  By default it is
## 0.05 corresponding to 95% confidence intervals.
##
## @item @qcode{'options'} @tab A structure specifying the control
## parameters for the iterative algorithm used to compute ML estimates with the
## @code{fminsearch} function.
##
## @item @qcode{'pdf'} @tab A function handle
## @code{@@(@var{data}, @var{p1}, @var{p2}, @dots{})} to the probability density
## of a @strong{custom} distribution, whose parameters are then estimated by
## maximum likelihood.  Requires @qcode{'start'}.  It is mutually exclusive with
## @qcode{'distribution'} and with @qcode{'logpdf'}/@qcode{'nloglf'}.
##
## @item @qcode{'cdf'} @tab A function handle to the cumulative distribution
## function of the custom distribution, with the same calling convention as
## @qcode{'pdf'}.  Required together with @qcode{'pdf'} for censored or
## truncated data.
##
## @item @qcode{'logpdf'} @tab A function handle to the log probability density
## of a custom distribution, with the same calling convention as @qcode{'pdf'}.
## Requires @qcode{'start'}.
##
## @item @qcode{'logsf'} @tab A function handle to the log survivor function
## @math{log (1 - cdf)} of the custom distribution, with the same calling
## convention as @qcode{'pdf'}.  Required together with @qcode{'logpdf'} for
## censored data.
##
## @item @qcode{'nloglf'} @tab A function handle
## @code{@@(@var{params}, @var{data}, @var{cens}, @var{freq})} returning the
## scalar negative log-likelihood of a custom distribution.  Requires
## @qcode{'start'}.
##
## @item @qcode{'start'} @tab A vector of initial parameter values for a
## custom-distribution fit.  Required with @qcode{'pdf'}, @qcode{'logpdf'}, or
## @qcode{'nloglf'}.
##
## @item @qcode{'lowerbound'} @tab A scalar or vector of lower bounds for the
## custom-distribution parameters.  By default they are unbounded below.
##
## @item @qcode{'upperbound'} @tab A scalar or vector of upper bounds for the
## custom-distribution parameters.  By default they are unbounded above.
##
## @item @qcode{'truncationbounds'} @tab A two-element vector @qcode{[L U]}
## giving the truncation interval of a custom distribution.  Requires a
## @qcode{'cdf'} function.
##
## @item @qcode{'optimfun'} @tab The optimizer for a custom-distribution fit.
## Only @qcode{'fminsearch'} is supported; bounded fits are handled by internal
## reparameterization of the constrained parameters.
## @end multitable
##
## When a custom distribution is specified through @qcode{'pdf'},
## @qcode{'logpdf'}, or @qcode{'nloglf'}, the parameters are estimated by
## maximizing the likelihood with @code{fminsearch}, and the second output
## @var{pci} gives asymptotic normal (Wald) confidence intervals computed from
## the observed Fisher information at @var{phat} (see @code{mlecov}).  Bounded
## parameters are estimated on an internally reparameterized unconstrained
## scale.
##
## @seealso{mlecov, fitdist, makedist}
## @end deftypefn

function [phat, pci] = mle (x, varargin)

  ## Check data
  if (! (isvector (x) && isnumeric (x) && isreal (x)))
    error ("mle: X must be a numeric vector of real values.");
  endif

  ## Add defaults
  censor = [];
  freq = ones (size (x));
  alpha = 0.05;
  ntrials = [];
  mu = 0;
  theta = 1;
  options.Display = 'off';
  options.MaxFunEvals = 400;
  options.MaxIter = 200;
  options.TolX = 1e-6;
  distname = 'normal';
  userdist = false;
  custpdf = [];
  custlogpdf = [];
  custnloglf = [];
  custcdf = [];
  custlogsf = [];
  start = [];
  lowbnd = [];
  uppbnd = [];
  truncbnd = [];
  optimfun = 'fminsearch';

  ## Parse extra arguments
  if (mod (numel (varargin), 2) != 0)
    error ("mle: optional arguments must be in NAME-VALUE pairs.");
  endif
  while (numel (varargin) > 0)
    switch (tolower (varargin{1}))
      case 'distribution'
        distname = varargin{2};
        userdist = true;
      case 'censoring'
        censor = varargin{2};
        if (! isequal (size (x), size (censor)) && ! isempty (censor))
          error (strcat ("mle: 'censoring' argument must have the same", ...
                         " size as the input data in X."));
        endif
      case 'frequency'
        freq = varargin{2};
        if (! isequal (size (x), size (freq)))
          error (strcat ("mle: 'frequency' argument must have the same", ...
                         " size as the input data in X."));
        endif
        if (any (freq != round (freq)) || any (freq < 0))
          error (strcat ("mle: 'frequency' argument must contain", ...
                         " non-negative integer values."));
        endif
      case 'alpha'
        alpha = varargin{2};
        if (! isscalar (alpha) || ! isreal (alpha) || alpha <= 0 || alpha >= 1)
          error ("mle: invalid value for 'alpha' argument.");
        endif
      case 'ntrials'
        ntrials = varargin{2};
        if (! (isscalar (ntrials) && isreal (ntrials) && ntrials > 0
                                  && fix (ntrials) == ntrials))
          error (strcat ("mle: 'ntrials' argument must be a positive", ...
                         " integer scalar value."));
        endif
      case {'mu'}
        mu = varargin{2};
      case {'theta'}
        theta = varargin{2};
      case 'options'
        options = varargin{2};
        if (! isstruct (options) || ! isfield (options, 'Display') ||
            ! isfield (options, 'MaxFunEvals') || ! isfield (options, 'MaxIter')
                                               || ! isfield (options, 'TolX'))
          error (strcat ("mle: 'options' argument must be a structure", ...
                         " compatible for 'fminsearch'."));
        endif

      case 'pdf'
        custpdf = varargin{2};
      case 'logpdf'
        custlogpdf = varargin{2};
      case 'nloglf'
        custnloglf = varargin{2};
      case 'cdf'
        custcdf = varargin{2};
      case 'logsf'
        custlogsf = varargin{2};
      case 'start'
        start = varargin{2};
      case 'lowerbound'
        lowbnd = varargin{2};
      case 'upperbound'
        uppbnd = varargin{2};
      case 'truncationbounds'
        truncbnd = varargin{2};
      case 'optimfun'
        optimfun = varargin{2};
      otherwise
        error ("mle: unknown parameter name.");
    endswitch
    varargin([1:2]) = [];
  endwhile

  ## Custom-distribution fit through user-supplied likelihood functions.  When a
  ## 'pdf', 'logpdf', or 'nloglf' handle is given, the parameters are estimated
  ## by maximum likelihood and the named-distribution path below is bypassed.
  ncust = (! isempty (custpdf)) + (! isempty (custlogpdf)) ...
                                + (! isempty (custnloglf));
  if (ncust > 0)
    if (userdist)
      error (strcat ("mle: the 'distribution' argument cannot be combined", ...
                     " with a custom 'pdf', 'logpdf', or 'nloglf' function."));
    endif
    if (ncust > 1)
      error (strcat ("mle: only one of the 'pdf', 'logpdf', or 'nloglf'", ...
                     " arguments can be specified."));
    endif
    [phat, pci] = mle_custom (x, custpdf, custlogpdf, custnloglf, custcdf, ...
                              custlogsf, start, lowbnd, uppbnd, truncbnd, ...
                              censor, freq, alpha, options, optimfun, nargout);
    return;
  endif

  ## Switch to known distributions
  switch (tolower (distname))

    case 'bernoulli'
      if (! isempty (censor))
        error (strcat ("mle: censoring is not supported for", ...
                       " the Bernoulli distribution."));
      elseif (any (x != 0 & x != 1))
        error ("mle: invalid data for the Bernoulli distribution.");
      endif
      if (! isempty (freq))
        x = expandFreq (x, freq);
      endif
      if (nargout < 2)
        phat = binofit (sum (x), numel (x));
      else
        [phat, pci] = binofit (sum (x), numel (x), alpha);
      endif

    case 'beta'
      if (! isempty (censor))
        error ("mle: censoring is not supported for the Beta distribution.");
      endif
      if (nargout < 2)
        phat = betafit (x, alpha, freq, options);
      else
        [phat, pci] = betafit (x, alpha, freq, options);
      endif

    case {'binomial', 'bino'}
      if (! isempty (censor))
        error (strcat ("mle: censoring is not supported for", ...
                       " the Binomial distribution."));
      elseif (isempty (ntrials))
        error (strcat ("mle: 'Ntrials' parameter is required", ...
                       " for the Binomial distribution."));
      endif
      if (nargout < 2)
        phat = binofit (sum (x .* freq), sum (freq) .* ntrials);
      else
        [phat, pci] = binofit (sum (x .* freq), sum (freq) .* ntrials, alpha);
      endif

    case {'bisa', 'BirnbaumSaunders'}
      if (nargout < 2)
        phat = bisafit (x, alpha, censor, freq, options);
      else
        [phat, pci] = bisafit (x, alpha, censor, freq, options);
      endif

    case 'burr'
      if (nargout < 2)
        phat = burrfit (x, alpha, censor, freq, options);
      else
        [phat, pci] = burrfit (x, alpha, censor, freq, options);
      endif

    case {'ev', 'extreme value'}
      if (nargout < 2)
        phat = evfit (x, alpha, censor, freq, options);
      else
        [phat, pci] = evfit (x, alpha, censor, freq, options);
      endif

    case {'exp', 'exponential'}
      if (nargout < 2)
        phat = expfit (x, alpha, censor, freq);
      else
        [phat, pci] = expfit (x, alpha, censor, freq);
      endif

    case {'gam', 'gamma'}
      if (nargout < 2)
        phat = gamfit (x, alpha, censor, freq, options);
      else
        [phat, pci] = gamfit (x, alpha, censor, freq, options);
      endif

    case {'geo', 'geometric'}
      if (! isempty (censor))
        error (strcat ("mle: censoring is not supported for the", ...
                       " Geometric distribution."));
      endif
      if (nargout < 2)
        phat = geofit (x, alpha, freq);
      else
        [phat, pci] = geofit (x, alpha, freq);
      endif

    case {'gev', 'generalized extreme value'}
      if (! isempty (censor))
        error (strcat ("mle: censoring is not supported for the", ...
                       " Generalized Extreme Value distribution."));
      endif
      if (nargout < 2)
        phat = gevfit (x, alpha, freq, options);
      else
        [phat, pci] = gevfit (x, alpha, freq, options);
      endif

    case {'gp', 'generalized pareto'}
      if (! isempty (censor))
        error (strcat ("mle: censoring is not supported for", ...
                       " the Generalized Pareto distribution."));
      endif
      if (any (x < theta))
        error (strcat ("mle: invalid 'theta' location parameter", ...
                       " for the Generalized Pareto distribution."));
      endif
      if (nargout < 2)
        phat = gpfit (x, theta, alpha, freq, options);
      else
        [phat, pci] = gpfit (x, theta, alpha, freq, options);
      endif

    case 'gumbel'
      if (nargout < 2)
        phat = gumbelfit (x, alpha, censor, freq, options);
      else
        [phat, pci] = gumbelfit (x, alpha, censor, freq, options);
      endif

    case {'hn', 'half normal', 'halfnormal'}
      if (! isempty (censor))
        error (strcat ("mle: censoring is not supported for", ...
                       " the Half Normal distribution."));
      endif
      if (any (x < mu))
        error (strcat ("mle: invalid 'mu' location parameter", ...
                       " for the Half Normal distribution."));
      endif
      if (nargout < 2)
        phat = hnfit (x, mu, alpha, freq);
      else
        [phat, pci] = hnfit (x, mu, alpha, freq);
      endif

    case {'invg', 'inversegaussian', 'inverse gaussian'}
      if (nargout < 2)
        phat = invgfit (x, alpha, censor, freq, options);
      else
        [phat, pci] = invgfit (x, alpha, censor, freq, options);
      endif

    case {'logi', 'logistic'}
      if (nargout < 2)
        phat = logifit (x, alpha, censor, freq, options);
      else
        [phat, pci] = logifit (x, alpha, censor, freq, options);
      endif

    case {'logl', 'loglogistic'}
      if (nargout < 2)
        phat = loglfit (x, alpha, censor, freq, options);
      else
        [phat, pci] = loglfit (x, alpha, censor, freq, options);
      endif

    case {'logn', 'lognormal'}
      if (nargout < 2)
        phat = lognfit (x, alpha, censor, freq, options);
      else
        [phat, pci] = lognfit (x, alpha, censor, freq, options);
      endif

    case {'naka', 'nakagami'}
      if (nargout < 2)
        phat = nakafit (x, alpha, censor, freq, options);
      else
        [phat, pci] = nakafit (x, alpha, censor, freq, options);
      endif

    case {'nbin', 'negativebinomial', 'negative binomial'}
      if (! isempty (censor))
        error (strcat ("mle: censoring is not supported for", ...
                       " the Negative Binomial distribution."));
      endif
      if (nargout < 2)
        phat = nbinfit (x, alpha, freq, options);
      else
        [phat, pci] = nbinfit (x, alpha, freq, options);
      endif

    case {'norm', 'normal'}
      if (nargout < 2)
        [muhat, sigmahat] = normfit (x, alpha, censor, freq, options);
        phat = [muhat, sigmahat];
      else
        [muhat, sigmahat, muci, sigmaci] = normfit (x, alpha, censor, ...
                                                    freq, options);
        phat = [muhat, sigmahat];
        pci = [muci, sigmaci];
      endif

    case {'poiss', 'poisson'}
      if (! isempty (censor))
        error (strcat ("mle: censoring is not supported for", ...
                       " the Poisson distribution."));
      endif
      if (nargout < 2)
        phat = poissfit (x, alpha, freq);
      else
        [phat, pci] = poissfit (x, alpha, freq);
      endif

    case {'rayl', 'rayleigh'}
      if (nargout < 2)
        phat = raylfit (x, alpha, censor, freq);
      else
        [phat, pci] = raylfit (x, alpha, censor, freq);
      endif

    case {'rice', 'rician'}
      if (nargout < 2)
        phat = ricefit (x, alpha, censor, freq, options);
      else
        [phat, pci] = ricefit (x, alpha, censor, freq, options);
      endif

    case {'tls', 'tlocationscale'}
      if (nargout < 2)
        phat = tlsfit (x, alpha, censor, freq, options);
      else
        [phat, pci] = tlsfit (x, alpha, censor, freq, options);
      endif

    case {'unid', 'uniform discrete', 'discrete uniform', 'discrete'}
      if (! isempty (censor))
        error (strcat ("mle: censoring is not supported for", ...
                       " the Discrete Uniform distribution."));
      endif
      if (nargout < 2)
        phat = unidfit (x, alpha, freq);
      else
        [phat, pci] = unidfit (x, alpha, freq);
      endif

    case {'unif', 'uniform', 'continuous uniform'}
      if (! isempty (censor))
        error (strcat ("mle: censoring is not supported for", ...
                       " the Continuous Uniform distribution."));
      endif
      if (nargout < 2)
        phat = uniffit (x, alpha, freq);
      else
        [phat, pci] = uniffit (x, alpha, freq);
      endif

    case {'wbl', 'weibull'}
      if (nargout < 2)
        phat = wblfit (x, alpha, censor, freq, options);
      else
        [phat, pci] = wblfit (x, alpha, censor, freq, options);
      endif

    otherwise
      error ("mle: unrecognized distribution name.");

  endswitch

endfunction

## Helper function for expanding data according to frequency vector
function [x, freq] = expandFreq (x, freq)
  ## Remove NaNs and zeros
  remove = isnan (freq);
  x(remove) = [];
  freq(remove) = [];
  if (! all (freq == 1))
    xf = [];
    for i = 1:numel (freq)
      xf = [xf, repmat(x(i), 1, freq(i))];
    endfor
  endif
  x = xf;
endfunction

## Maximum likelihood fit of a user-supplied custom distribution
function [phat, pci] = mle_custom (x, custpdf, custlogpdf, custnloglf, ...
                        custcdf, custlogsf, start, lowbnd, uppbnd, truncbnd, ...
                        censor, freq, alpha, options, optimfun, nout)

  ## Only fminsearch is available in core Octave; bounded fits are handled by
  ## reparameterization, so a separate constrained optimizer is not required.
  if (! (ischar (optimfun) && strcmpi (optimfun, 'fminsearch')))
    error (strcat ("mle: 'optimfun' only supports 'fminsearch'; bounded", ...
                   " fits are handled by internal reparameterization."));
  endif

  ## 'start' is required and sets the number of parameters
  if (isempty (start))
    error (strcat ("mle: a 'start' vector of initial parameter values is", ...
                   " required for a custom distribution fit."));
  endif
  if (! (isvector (start) && isnumeric (start) && isreal (start)))
    error ("mle: 'start' must be a numeric vector of real values.");
  endif
  start = start(:).';
  k = numel (start);

  ## Identify the input form and validate the supplied handles
  if (! isempty (custpdf))
    form = 'pdf';
    if (! is_function_handle (custpdf))
      error ("mle: 'pdf' argument must be a function handle.");
    endif
  elseif (! isempty (custlogpdf))
    form = 'logpdf';
    if (! is_function_handle (custlogpdf))
      error ("mle: 'logpdf' argument must be a function handle.");
    endif
  else
    form = 'nloglf';
    if (! is_function_handle (custnloglf))
      error ("mle: 'nloglf' argument must be a function handle.");
    endif
  endif
  if (! isempty (custcdf) && ! is_function_handle (custcdf))
    error ("mle: 'cdf' argument must be a function handle.");
  endif
  if (! isempty (custlogsf) && ! is_function_handle (custlogsf))
    error ("mle: 'logsf' argument must be a function handle.");
  endif

  ## Data, censoring, and frequency as columns
  x = x(:);
  n = numel (x);
  if (isempty (censor))
    cens = false (n, 1);
  else
    cens = logical (censor(:));
  endif
  freq = freq(:);
  docens = any (cens);

  ## Truncation bounds
  dotrunc = ! isempty (truncbnd);
  if (dotrunc && ! (isnumeric (truncbnd) && isreal (truncbnd)
                    && numel (truncbnd) == 2 && truncbnd(1) < truncbnd(2)))
    error (strcat ("mle: 'truncationbounds' must be a two-element vector", ...
                   " [L U] with L < U."));
  endif

  ## Censoring and truncation need the complementary functions
  if (strcmp (form, 'pdf') && (docens || dotrunc) && isempty (custcdf))
    error (strcat ("mle: a 'cdf' function handle is required for censored", ...
                   " or truncated data when using the 'pdf' argument."));
  endif
  if (strcmp (form, 'logpdf') && docens && isempty (custlogsf))
    error (strcat ("mle: a 'logsf' function handle is required for", ...
                   " censored data when using the 'logpdf' argument."));
  endif
  if (strcmp (form, 'logpdf') && dotrunc && isempty (custcdf))
    error (strcat ("mle: a 'cdf' function handle is required for truncated", ...
                   " data when using the 'logpdf' argument."));
  endif

  ## Parameter bounds, expanded to per-parameter row vectors
  if (isempty (lowbnd))
    lb = -Inf (1, k);
  elseif (isscalar (lowbnd))
    lb = lowbnd * ones (1, k);
  elseif (numel (lowbnd) == k)
    lb = lowbnd(:).';
  else
    error ("mle: 'lowerbound' must be a scalar or match the size of 'start'.");
  endif
  if (isempty (uppbnd))
    ub = Inf (1, k);
  elseif (isscalar (uppbnd))
    ub = uppbnd * ones (1, k);
  elseif (numel (uppbnd) == k)
    ub = uppbnd(:).';
  else
    error ("mle: 'upperbound' must be a scalar or match the size of 'start'.");
  endif
  if (any (lb >= ub))
    error (strcat ("mle: each 'lowerbound' must be strictly less than its", ...
                   " corresponding 'upperbound'."));
  endif
  if (any (start <= lb) || any (start >= ub))
    error ("mle: 'start' values must lie strictly within the given bounds.");
  endif

  ## Aggregate negative log-likelihood of the sample at a parameter row vector
  nllfun = @(th) custom_nll (th, form, custpdf, custlogpdf, custnloglf, ...
                             custcdf, custlogsf, x, cens, freq, docens, ...
                             dotrunc, truncbnd);

  ## Maximize the likelihood.  Bounded parameters are optimized on an
  ## unconstrained transformed scale and mapped back.
  if (any (isfinite (lb)) || any (isfinite (ub)))
    obj = @(u) nllfun (to_con (u, lb, ub));
    uopt = fminsearch (obj, to_uncon (start, lb, ub), options);
    phat = to_con (uopt, lb, ub);
  else
    phat = fminsearch (nllfun, start, options);
  endif

  ## Asymptotic (Wald) confidence intervals from the observed information
  if (nout > 1)
    acov = mlecov (phat, x, 'nloglf', @(pp, dd, cc, ff) nllfun (pp));
    se = sqrt (diag (acov)).';
    z = norminv (1 - alpha / 2);
    pci = [phat - z .* se; phat + z .* se];
  else
    pci = [];
  endif

endfunction

## Aggregate negative log-likelihood for a custom distribution
function nll = custom_nll (th, form, cpdf, clogpdf, cnloglf, ccdf, clogsf, ...
                           x, cens, freq, docens, dotrunc, tb)
  if (strcmp (form, 'nloglf'))
    nll = cnloglf (th, x, double (cens), freq);
    return;
  endif
  pc = num2cell (th);
  unc = ! cens;
  terms = zeros (size (x));
  if (strcmp (form, 'pdf'))
    terms(unc) = log (cpdf (x(unc), pc{:}));
    if (docens)
      terms(cens) = log (1 - ccdf (x(cens), pc{:}));
    endif
  else
    terms(unc) = clogpdf (x(unc), pc{:});
    if (docens)
      terms(cens) = clogsf (x(cens), pc{:});
    endif
  endif
  if (dotrunc)
    logZ = log (ccdf (tb(2), pc{:}) - ccdf (tb(1), pc{:}));
    ## Right-censored survival is renormalized to the truncation interval
    if (docens)
      terms(cens) = log (ccdf (tb(2), pc{:}) - ccdf (x(cens), pc{:}));
    endif
    terms = terms - logZ;
  endif
  nll = -sum (freq .* terms);
endfunction

## Map constrained parameters to an unconstrained scale for optimization
function u = to_uncon (th, lb, ub)
  u = th;
  for i = 1:numel (th)
    if (isfinite (lb(i)) && isfinite (ub(i)))
      u(i) = log ((th(i) - lb(i)) / (ub(i) - th(i)));
    elseif (isfinite (lb(i)))
      u(i) = log (th(i) - lb(i));
    elseif (isfinite (ub(i)))
      u(i) = log (ub(i) - th(i));
    endif
  endfor
endfunction

## Map unconstrained parameters back to the constrained scale
function th = to_con (u, lb, ub)
  th = u;
  for i = 1:numel (u)
    if (isfinite (lb(i)) && isfinite (ub(i)))
      th(i) = lb(i) + (ub(i) - lb(i)) / (1 + exp (-u(i)));
    elseif (isfinite (lb(i)))
      th(i) = lb(i) + exp (u(i));
    elseif (isfinite (ub(i)))
      th(i) = ub(i) - exp (u(i));
    endif
  endfor
endfunction

%!demo
%! ## Fit a custom (normal) distribution by maximum likelihood and return the
%! ## asymptotic 95% confidence intervals of the estimates.
%! x = [2.1, 3.4, 1.9, 5.2, 4.1, 2.8, 3.3, 4.7, 2.2, 3.9, 3.0, 4.5];
%! pdf = @(x, mu, sigma) normpdf (x, mu, sigma);
%! [phat, pci] = mle (x, 'pdf', pdf, 'start', [mean(x), std(x)])

## Test custom-distribution fitting (values verified against MATLAB)
%!test
%! x = [2.1, 3.4, 1.9, 5.2, 4.1, 2.8, 3.3, 4.7, 2.2, 3.9, 3.0, 4.5];
%! [phat, pci] = mle (x, 'pdf', @(x, mu, s) normpdf (x, mu, s), ...
%!                    'start', [mean(x), std(x)]);
%! assert (phat, [3.42499970800201, 1.03208912390818], 1e-4);
%! assert (pci, [2.8410510441517, 0.619175174501268; ...
%!               4.00894837185232, 1.44500307331509], 1e-4);

%!test
%! x = [2.1, 3.4, 1.9, 5.2, 4.1, 2.8, 3.3, 4.7, 2.2, 3.9, 3.0, 4.5];
%! nll = @(p, d, c, f) -sum (log (normpdf (d, p(1), p(2))));
%! phat = mle (x, 'nloglf', nll, 'start', [3, 1]);
%! assert (phat, [3.42499959294639, 1.03208933560634], 1e-4);

%!test
%! x = [2.1, 3.4, 1.9, 5.2, 4.1, 2.8, 3.3, 4.7, 2.2, 3.9, 3.0, 4.5];
%! phat = mle (x, 'logpdf', @(x, mu, s) log (normpdf (x, mu, s)), ...
%!             'start', [3, 1]);
%! assert (phat, [3.42499959294639, 1.03208933560634], 1e-4);

%!test
%! ## Alpha propagates into the Wald interval
%! x = [2.1, 3.4, 1.9, 5.2, 4.1, 2.8, 3.3, 4.7, 2.2, 3.9, 3.0, 4.5];
%! [~, pci] = mle (x, 'pdf', @(x, mu, s) normpdf (x, mu, s), ...
%!                 'start', [mean(x), std(x)], 'alpha', 0.10);
%! assert (pci, [2.93493454085403, 0.685560813868748; ...
%!               3.91506487514999, 1.37861743394761], 1e-4);

%!test
%! ## Frequency-weighted fit
%! x = [2.1, 3.4, 1.9, 5.2, 4.1, 2.8, 3.3, 4.7, 2.2, 3.9, 3.0, 4.5];
%! f = [1, 2, 1, 1, 3, 1, 2, 1, 1, 1, 2, 1];
%! phat = mle (x, 'pdf', @(x, mu, s) normpdf (x, mu, s), ...
%!             'start', [mean(x), std(x)], 'frequency', f);
%! assert (phat, [3.47058837030174, 0.902783454993327], 1e-4);

%!test
%! ## Right-censored fit (needs a cdf)
%! x = [2.1, 3.4, 1.9, 5.2, 4.1, 2.8, 3.3, 4.7, 2.2, 3.9, 3.0, 4.5];
%! c = [0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0];
%! phat = mle (x, 'pdf', @(x, mu, s) normpdf (x, mu, s), ...
%!             'cdf', @(x, mu, s) normcdf (x, mu, s), ...
%!             'start', [mean(x), std(x)], 'censoring', c);
%! assert (phat, [3.62195934926527, 1.03307264832117], 1e-4);

%!test
%! ## Bounded fit via reparameterization (lower bound inactive at the optimum)
%! x = [2.1, 3.4, 1.9, 5.2, 4.1, 2.8, 3.3, 4.7, 2.2, 3.9, 3.0, 4.5];
%! phat = mle (x, 'pdf', @(x, mu) exppdf (x, mu), 'start', 3, 'lowerbound', 0);
%! assert (phat, 3.42499980926514, 1e-4);

%!test
%! ## Truncated fit on [1, 6]
%! x = [2.1, 3.4, 1.9, 5.2, 4.1, 2.8, 3.3, 4.7, 2.2, 3.9, 3.0, 4.5];
%! phat = mle (x, 'pdf', @(x, mu, s) normpdf (x, mu, s), ...
%!             'cdf', @(x, mu, s) normcdf (x, mu, s), ...
%!             'start', [mean(x), std(x)], 'truncationbounds', [1, 6]);
%! assert (phat, [3.41148048015442, 1.12176348358835], 1e-4);

## Test input validation
%!error <mle: X must be a numeric vector of real values.> mle (ones (2))
%!error <mle: X must be a numeric vector of real values.> mle ('text')
%!error <mle: X must be a numeric vector of real values.> mle ([1, 2, 3, i, 5])
%!error <mle: optional arguments must be in NAME-VALUE pairs.> ...
%! mle ([1:50], 'distribution')
%!error <mle: 'censoring' argument must have the same size as the input data in X.> ...
%! mle ([1:50], 'censoring', logical ([1,0,1,0]))
%!error <mle: 'frequency' argument must have the same size as the input data in X.> ...
%! mle ([1:50], 'frequency', [1,0,1,0])
%!error <mle: 'frequency' argument must contain non-negative integer values.> ...
%! mle ([1 0 1 0], 'frequency', [-1 1 0 0])
%!error <mle: 'frequency' argument must contain non-negative integer values.> ...
%! mle ([1 0 1 0], 'distribution', 'nbin', 'frequency', [-1 1 0 0])
%!error <mle: invalid value for 'alpha' argument.> mle ([1:50], 'alpha', [0.05, 0.01])
%!error <mle: invalid value for 'alpha' argument.> mle ([1:50], 'alpha', 1)
%!error <mle: invalid value for 'alpha' argument.> mle ([1:50], 'alpha', -1)
%!error <mle: invalid value for 'alpha' argument.> mle ([1:50], 'alpha', i)
%!error <mle: 'ntrials' argument must be a positive integer scalar value.> ...
%! mle ([1:50], 'ntrials', -1)
%!error <mle: 'ntrials' argument must be a positive integer scalar value.> ...
%! mle ([1:50], 'ntrials', [20, 50])
%!error <mle: 'ntrials' argument must be a positive integer scalar value.> ...
%! mle ([1:50], 'ntrials', [20.3])
%!error <mle: 'ntrials' argument must be a positive integer scalar value.> ...
%! mle ([1:50], 'ntrials', 3i)
%!error <mle: 'options' argument must be a structure compatible for 'fminsearch'.> ...
%! mle ([1:50], 'options', 4)
%!error <mle: 'options' argument must be a structure compatible for 'fminsearch'.> ...
%! mle ([1:50], 'options', struct ('x', 3))
%!error <mle: unknown parameter name.> mle ([1:50], 'NAME', 'value')
%!error <mle: censoring is not supported for the Bernoulli distribution.> ...
%! mle ([1 0 1 0], 'distribution', 'bernoulli', 'censoring', [1 1 0 0])
%!error <mle: invalid data for the Bernoulli distribution.> ...
%! mle ([1 2 1 0], 'distribution', 'bernoulli')
%!error <mle: censoring is not supported for the Beta distribution.> ...
%! mle ([1 0 1 0], 'distribution', 'beta', 'censoring', [1 1 0 0])
%!error <mle: censoring is not supported for the Binomial distribution.> ...
%! mle ([1 0 1 0], 'distribution', 'bino', 'censoring', [1 1 0 0])
%!error <mle: 'Ntrials' parameter is required for the Binomial distribution.> ...
%! mle ([1 0 1 0], 'distribution', 'bino')
%!error <mle: censoring is not supported for the Geometric distribution.> ...
%! mle ([1 0 1 0], 'distribution', 'geo', 'censoring', [1 1 0 0])
%!error <mle: censoring is not supported for the Generalized Extreme Value distribution.> ...
%! mle ([1 0 1 0], 'distribution', 'gev', 'censoring', [1 1 0 0])
%!error <mle: censoring is not supported for the Generalized Pareto distribution.> ...
%! mle ([1 0 1 0], 'distribution', 'gp', 'censoring', [1 1 0 0])
%!error <mle: invalid 'theta' location parameter for the Generalized Pareto distribution.> ...
%! mle ([1 0 -1 0], 'distribution', 'gp')
%!error <mle: censoring is not supported for the Half Normal distribution.> ...
%! mle ([1 0 1 0], 'distribution', 'hn', 'censoring', [1 1 0 0])
%!error <mle: invalid 'mu' location parameter for the Half Normal distribution.> ...
%! mle ([1 0 -1 0], 'distribution', 'hn')
%!error <mle: censoring is not supported for the Negative Binomial distribution.> ...
%! mle ([1 0 1 0], 'distribution', 'nbin', 'censoring', [1 1 0 0])
%!error <mle: censoring is not supported for the Poisson distribution.> ...
%! mle ([1 0 1 0], 'distribution', 'poisson', 'censoring', [1 1 0 0])
%!error <mle: censoring is not supported for the Discrete Uniform distribution.> ...
%! mle ([1 0 1 0], 'distribution', 'unid', 'censoring', [1 1 0 0])
%!error <mle: censoring is not supported for the Continuous Uniform distribution.> ...
%! mle ([1 0 1 0], 'distribution', 'unif', 'censoring', [1 1 0 0])
%!error <mle: unrecognized distribution name.> mle ([1:50], 'distribution', 'value')
%!error <mle: censoring is not supported for the Continuous Uniform distribution.> ...
%! mle ([1 0 1 0], 'distribution', 'unif', 'censoring', [1 1 0 0])
%!error <mle: the 'distribution' argument cannot be combined with a custom 'pdf', 'logpdf', or 'nloglf' function.> ...
%! mle ([1:50], 'distribution', 'normal', 'pdf', @(x, a, b) normpdf (x, a, b))
%!error <mle: only one of the 'pdf', 'logpdf', or 'nloglf' arguments can be specified.> ...
%! mle ([1:50], 'pdf', @sin, 'nloglf', @cos)
%!error <mle: a 'start' vector of initial parameter values is required for a custom distribution fit.> ...
%! mle ([1:50], 'pdf', @(x, a, b) normpdf (x, a, b))
%!error <mle: 'start' must be a numeric vector of real values.> ...
%! mle ([1:50], 'pdf', @(x, a, b) normpdf (x, a, b), 'start', 'text')
%!error <mle: 'pdf' argument must be a function handle.> ...
%! mle ([1:50], 'pdf', 5, 'start', [0, 1])
%!error <mle: 'nloglf' argument must be a function handle.> ...
%! mle ([1:50], 'nloglf', 5, 'start', [0, 1])
%!error <mle: 'cdf' argument must be a function handle.> ...
%! mle ([1:50], 'pdf', @(x, a, b) normpdf (x, a, b), 'cdf', 5, 'start', [0, 1])
%!error <mle: a 'cdf' function handle is required for censored or truncated data when using the 'pdf' argument.> ...
%! mle ([1:50], 'pdf', @(x, a, b) normpdf (x, a, b), 'start', [0, 1], ...
%!      'censoring', [1, zeros(1, 49)])
%!error <mle: 'truncationbounds' must be a two-element vector \[L U\] with L < U.> ...
%! mle ([1:50], 'pdf', @(x, a, b) normpdf (x, a, b), ...
%!      'cdf', @(x, a, b) normcdf (x, a, b), 'start', [0, 1], ...
%!      'truncationbounds', [6, 1])
%!error <mle: each 'lowerbound' must be strictly less than its corresponding 'upperbound'.> ...
%! mle ([1:50], 'pdf', @(x, a, b) normpdf (x, a, b), 'start', [0, 1], ...
%!      'lowerbound', [0, 2], 'upperbound', [0, 5])
%!error <mle: 'start' values must lie strictly within the given bounds.> ...
%! mle ([1:50], 'pdf', @(x, a, b) normpdf (x, a, b), 'start', [0, 1], ...
%!      'lowerbound', [1, 0])
%!error <mle: 'optimfun' only supports 'fminsearch'; bounded fits are handled by internal reparameterization.> ...
%! mle ([1:50], 'pdf', @(x, a, b) normpdf (x, a, b), 'start', [0, 1], ...
%!      'optimfun', 'fmincon')
