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

function [phat, pci] = mle (x, varargin)

  ## Check data
  if (! (isvector (x) && isnumeric (x) && isreal (x)))
    error ("mle: X must be a numeric vector of real values.");
  endif

  ## Add defaults
  censor = [];
  freq = [];
  alpha = 0.05;
  ntrials = [];
  mu = 0;
  options.Display = "off";
  options.MaxFunEvals = 400;
  options.MaxIter = 200;
  options.TolX = 1e-6;
  distname = "normal";

  ## Parse extra arguments
  if (mod (numel (varargin), 2) != 0)
    error ("mle: optional arguments must be in NAME-VALUE pairs.");
  endif
  while (numel (varargin) > 0)
    switch (tolower (varargin{1}))
      case "distribution"
        distname = varargin{2};
      case "censoring"
        censor = varargin{2};
        if (! isequal (size (x), size (censor)) && ! isempty (censor))
          error (strcat (["mle: 'censoring' argument must have the same"], ...
                         [" size as the input data in X."]));
        endif
      case "frequency"
        freq = varargin{2};
        if (! isequal (size (x), size (freq)) && ! isempty (freq))
          error (strcat (["mle: 'frequency' argument must have the same"], ...
                         [" size as the input data in X."]));
        endif
        if (any (freq != round (freq)) || any (freq < 0))
          error (strcat (["mle: 'frequency' argument must contain"], ...
                         [" non-negative integer values."]));
        endif
      case "alpha"
        alpha = varargin{2};
        if (! isscalar (alpha) || ! isreal (alpha) || alpha <= 0 || alpha >= 1)
          error ("mle: invalid value for 'alpha' argument.");
        endif
      case "ntrials"
        ntrials = varargin{2};
        if (! (isscalar (ntrials) && isreal (ntrials) && ntrials > 0
                                  && fix (ntrials) == ntrials))
          error (strcat (["mle: 'ntrials' argument must be a positive"], ...
                         [" integer scalar value."]));
        endif
      case {"theta", "mu"}
        mu = varargin{2};
      case "options"
        options = varargin{2};
        if (! isstruct (options) || ! isfield (options, "Display") ||
            ! isfield (options, "MaxFunEvals") || ! isfield (options, "MaxIter")
                                               || ! isfield (options, "TolX"))
          error (strcat (["mle: 'options' argument must be a structure"], ...
                         [" compatible for 'fminsearch'."]));
        endif

      case {"pdf", "cdf", "logpdf", "logsf", "nloglf", "truncationbounds", ...
            "start", "lowerbound", "upperbound", "optimfun"}
        printf ("mle: parameter not supported yet.");
      otherwise
        error ("mle: unknown parameter name.");
    endswitch
    varargin([1:2]) = [];
  endwhile

  ## Switch to known distributions
  switch (tolower (distname))

    case "bernoulli"
      if (! isempty (censor))
        error (strcat (["mle: censoring is not supported for"], ...
                       [" the Bernoulli distribution."]));
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

    case "beta"
      if (! isempty (censor))
        error ("mle: censoring is not supported for the Beta distribution.");
      endif
      if (! isempty (freq))
        x = expandFreq (x, freq);
      endif
      if (nargout < 2)
        phat = betafit (x, alpha, options);
      else
        [phat, pci] = betafit (x, alpha, options);
      endif

    case {"binomial", "bino"}
      if (! isempty (censor))
        error (strcat (["mle: censoring is not supported for"], ...
                       [" the Binomial distribution."]));
      elseif (isempty (ntrials))
        error (strcat (["mle: 'Ntrials' parameter is required"], ...
                       [" for the Binomial distribution."]));
      endif
      if (! isempty (freq))
        x = expandFreq (x, freq);
      endif
      if (nargout < 2)
        phat = binofit (sum (x), ntrials);
      else
        [phat, pci] = binofit (sum (x), ntrials, alpha);
      endif

    case {"bisa", "BirnbaumSaunders"}
      if (nargout < 2)
        phat = bisafit (x, alpha, censor, freq, options);
      else
        [phat, pci] = bisafit (x, alpha, censor, freq, options);
      endif

    case "burr"
      if (nargout < 2)
        phat = burrfit (x, alpha, censor, freq, options);
      else
        [phat, pci] = burrfit (x, alpha, censor, freq, options);
      endif

    case {"ev", "extreme value"}
      if (nargout < 2)
        phat = evfit (x, alpha, censor, freq, options);
      else
        [phat, pci] = evfit (x, alpha, censor, freq, options);
      endif

    case {"exp", "exponential"}
      if (nargout < 2)
        phat = expfit (x, alpha, censor, freq, options);
      else
        [phat, pci] = expfit (x, alpha, censor, freq, options);
      endif

    case {"gam", "gamma"}
      if (nargout < 2)
        phat = gamfit (x, alpha, censor, freq, options);
      else
        [phat, pci] = gamfit (x, alpha, censor, freq, options);
      endif

    case {"geo", "geometric"}
      if (! isempty (censor))
        error (strcat (["mle: censoring is not supported for the"], ...
                       [" Geometric distribution."]));
      endif
      if (nargout < 2)
        phat = gevfit (x, alpha, freq);
      else
        [phat, pci] = gevfit (x, alpha, freq);
      endif

    case {"gev", "generalized extreme value"}
      if (! isempty (censor))
        error (strcat (["mle: censoring is not supported for the"], ...
                       [" Generalized Extreme Value distribution."]));
      endif
      if (! isempty (freq))
        x = expandFreq (x, freq);
      endif
      if (nargout < 2)
        phat = gevfit (x, alpha, options);
      else
        [phat, pci] = gevfit (x, alpha, options);
      endif

    case {"gp", "generalized pareto"}
      if (! isempty (censor))
        error (strcat (["mle: censoring is not supported for"], ...
                       [" the Generalized Pareto distribution."]));
      endif
      if (! isempty (freq))
        x = expandFreq (x, freq);
      endif
      if (any (x < mu))
        error (strcat (["mle: invalid 'theta' location parameter"], ...
                       [" for the Generalized Pareto distribution."]));
      endif
      if (nargout < 2)
        phat = gevfit (x - mu, alpha, options);
      else
        [phat, pci] = gevfit (x - mu, alpha, options);
      endif

    case "gumbel"
      if (nargout < 2)
        phat = gumbelfit (x, alpha, censor, freq, options);
      else
        [phat, pci] = gumbelfit (x, alpha, censor, freq, options);
      endif

    case {"hn", "half normal", "halfnormal"}
      if (! isempty (censor))
        error (strcat (["mle: censoring is not supported for"], ...
                       [" the Half Normal distribution."]));
      endif
      if (any (x < mu))
        error (strcat (["mle: invalid 'mu' location parameter"], ...
                       [" for the Half Normal distribution."]));
      endif
      if (nargout < 2)
        phat = hnfit (x, mu, alpha, freq);
      else
        [phat, pci] = hnfit (x, mu, alpha, freq);
      endif

    case {"invg", "inversegaussian", "inverse gaussian"}
      if (nargout < 2)
        phat = invgfit (x, alpha, censor, freq, options);
      else
        [phat, pci] = invgfit (x, alpha, censor, freq, options);
      endif

    case {"logi", "logistic"}
      if (nargout < 2)
        phat = logifit (x, alpha, censor, freq, options);
      else
        [phat, pci] = logifit (x, alpha, censor, freq, options);
      endif

    case {"logl", "loglogistic"}
      if (nargout < 2)
        phat = loglfit (x, alpha, censor, freq, options);
      else
        [phat, pci] = loglfit (x, alpha, censor, freq, options);
      endif

    case {"logn", "lognormal"}
      if (nargout < 2)
        phat = lognfit (x, alpha, censor, freq, options);
      else
        [phat, pci] = lognfit (x, alpha, censor, freq, options);
      endif

    case {"naka", "nakagami"}
      if (nargout < 2)
        phat = nakafit (x, alpha, censor, freq, options);
      else
        [phat, pci] = nakafit (x, alpha, censor, freq, options);
      endif

    case {"nbin", "negativebinomial", "negative binomial"}
      if (! isempty (censor))
        error (strcat (["mle: censoring is not supported for"], ...
                       [" the Negative Binomial distribution."]));
      endif
      if (nargout < 2)
        phat = nbinfit (x, alpha, freq, options);
      else
        [phat, pci] = nbinfit (x, alpha, freq, options);
      endif

    case {"norm", "normal"}
      if (nargout < 2)
        [muhat, sigmahat] = normfit (x, alpha, censor, freq, options);
        phat = [muhat, sigmahat];
      else
        [muhat, sigmahat, muci, sigmaci] = normfit (x, alpha, censor, ...
                                                    freq, options);
        phat = [muhat, sigmahat];
        pci = [muci, sigmaci];
      endif

    case {"poiss", "poisson"}
      if (! isempty (censor))
        error (strcat (["mle: censoring is not supported for"], ...
                       [" the Poisson distribution."]));
      endif
      if (nargout < 2)
        phat = poissfit (x, alpha, freq);
      else
        [phat, pci] = poissfit (x, alpha, freq);
      endif

    case {"rayl", "rayleigh"}
      if (nargout < 2)
        phat = raylfit (x, alpha, censor, freq);
      else
        [phat, pci] = raylfit (x, alpha, censor, freq);
      endif

    case {"rice", "rician"}
      if (nargout < 2)
        phat = ricefit (x, alpha, censor, freq, options);
      else
        [phat, pci] = ricefit (x, alpha, censor, freq, options);
      endif

    case {"tls", "tlocationscale"}
      if (nargout < 2)
        phat = tlsfit (x, alpha, censor, freq, options);
      else
        [phat, pci] = tlsfit (x, alpha, censor, freq, options);
      endif

    case {"unid", "uniform discrete", "discrete uniform", "discrete"}
      if (! isempty (censor))
        error (strcat (["mle: censoring is not supported for"], ...
                       [" the Discrete Uniform distribution."]));
      endif
      if (nargout < 2)
        phat = unidfit (x, alpha, freq);
      else
        [phat, pci] = unidfit (x, alpha, freq);
      endif

    case {"unif", "uniform", "continuous uniform"}
      if (! isempty (censor))
        error (strcat (["mle: censoring is not supported for"], ...
                       [" the Continuous Uniform distribution."]));
      endif
      if (nargout < 2)
        phat = uniffit (x, alpha, freq);
      else
        [phat, pci] = uniffit (x, alpha, freq);
      endif

    case {"wbl", "weibull"}
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
function xf = expandFreq (x, freq)
  xf = [];
  for i = 1:numel (freq)
    xf = [xf, repmat(x(i), 1, freq(i))];
  endfor
endfunction

## Test input validation
%!error <mle: X must be a numeric vector of real values.> mle (ones (2))
%!error <mle: X must be a numeric vector of real values.> mle ("text")
%!error <mle: X must be a numeric vector of real values.> mle ([1, 2, 3, i, 5])
%!error <mle: optional arguments must be in NAME-VALUE pairs.> ...
%! mle ([1:50], "distribution")
%!error <mle: 'censoring' argument must have the same size as the input data in X.> ...
%! mle ([1:50], "censoring", logical ([1,0,1,0]))
%!error <mle: 'frequency' argument must have the same size as the input data in X.> ...
%! mle ([1:50], "frequency", [1,0,1,0])
%!error <mle: 'frequency' argument must contain non-negative integer values.> ...
%! mle ([1 0 1 0], "frequency", [-1 1 0 0])
%!error <mle: 'frequency' argument must contain non-negative integer values.> ...
%! mle ([1 0 1 0], "distribution", "nbin", "frequency", [-1 1 0 0])
%!error <mle: invalid value for 'alpha' argument.> mle ([1:50], "alpha", [0.05, 0.01])
%!error <mle: invalid value for 'alpha' argument.> mle ([1:50], "alpha", 1)
%!error <mle: invalid value for 'alpha' argument.> mle ([1:50], "alpha", -1)
%!error <mle: invalid value for 'alpha' argument.> mle ([1:50], "alpha", i)
%!error <mle: 'ntrials' argument must be a positive integer scalar value.> ...
%! mle ([1:50], "ntrials", -1)
%!error <mle: 'ntrials' argument must be a positive integer scalar value.> ...
%! mle ([1:50], "ntrials", [20, 50])
%!error <mle: 'ntrials' argument must be a positive integer scalar value.> ...
%! mle ([1:50], "ntrials", [20.3])
%!error <mle: 'ntrials' argument must be a positive integer scalar value.> ...
%! mle ([1:50], "ntrials", 3i)
%!error <mle: 'options' argument must be a structure compatible for 'fminsearch'.> ...
%! mle ([1:50], "options", 4)
%!error <mle: 'options' argument must be a structure compatible for 'fminsearch'.> ...
%! mle ([1:50], "options", struct ("x", 3))
%!error <mle: unknown parameter name.> mle ([1:50], "NAME", "value")
%!error <mle: censoring is not supported for the Bernoulli distribution.> ...
%! mle ([1 0 1 0], "distribution", "bernoulli", "censoring", [1 1 0 0])
%!error <mle: invalid data for the Bernoulli distribution.> ...
%! mle ([1 2 1 0], "distribution", "bernoulli")
%!error <mle: censoring is not supported for the Beta distribution.> ...
%! mle ([1 0 1 0], "distribution", "beta", "censoring", [1 1 0 0])
%!error <mle: censoring is not supported for the Binomial distribution.> ...
%! mle ([1 0 1 0], "distribution", "bino", "censoring", [1 1 0 0])
%!error <mle: 'Ntrials' parameter is required for the Binomial distribution.> ...
%! mle ([1 0 1 0], "distribution", "bino")
%!error <mle: censoring is not supported for the Geometric distribution.> ...
%! mle ([1 0 1 0], "distribution", "geo", "censoring", [1 1 0 0])
%!error <mle: censoring is not supported for the Generalized Extreme Value distribution.> ...
%! mle ([1 0 1 0], "distribution", "gev", "censoring", [1 1 0 0])
%!error <mle: censoring is not supported for the Generalized Pareto distribution.> ...
%! mle ([1 0 1 0], "distribution", "gp", "censoring", [1 1 0 0])
%!error <mle: invalid 'theta' location parameter for the Generalized Pareto distribution.> ...
%! mle ([1 0 -1 0], "distribution", "gp")
%!error <mle: censoring is not supported for the Half Normal distribution.> ...
%! mle ([1 0 1 0], "distribution", "hn", "censoring", [1 1 0 0])
%!error <mle: invalid 'mu' location parameter for the Half Normal distribution.> ...
%! mle ([1 0 -1 0], "distribution", "hn")
%!error <mle: censoring is not supported for the Negative Binomial distribution.> ...
%! mle ([1 0 1 0], "distribution", "nbin", "censoring", [1 1 0 0])
%!error <mle: censoring is not supported for the Poisson distribution.> ...
%! mle ([1 0 1 0], "distribution", "poisson", "censoring", [1 1 0 0])
%!error <mle: censoring is not supported for the Discrete Uniform distribution.> ...
%! mle ([1 0 1 0], "distribution", "unid", "censoring", [1 1 0 0])
%!error <mle: censoring is not supported for the Continuous Uniform distribution.> ...
%! mle ([1 0 1 0], "distribution", "unif", "censoring", [1 1 0 0])
%!error <mle: unrecognized distribution name.> mle ([1:50], "distribution", "value")
%!error <mle: censoring is not supported for the Continuous Uniform distribution.> ...
%! mle ([1 0 1 0], "distribution", "unif", "censoring", [1 1 0 0])
