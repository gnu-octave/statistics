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
## @deftypefn  {statistics} {@var{pd} =} fitdist (@var{x}, @var{distname})
## @deftypefnx {statistics} {@var{pd} =} fitdist (@var{x}, @var{distname}, @var{Name}, @var{Value})
## @deftypefnx {statistics} {[@var{pdca}, @var{gn}, @var{gl}] =} fitdist (@var{x}, @var{distname}, @qcode{"By"}, @var{groupvar})
## @deftypefnx {statistics} {[@var{pdca}, @var{gn}, @var{gl}] =} fitdist (@var{x}, @var{distname}, @qcode{"By"}, @var{groupvar}, @var{Name}, @var{Value})
##
## Create probability distribution object.
##
## @code{@var{pd} = fitdist (@var{x}, @var{distname})} creates a probability
## distribution distribution object by fitting the distribution specified by
## @var{distname} to the data in vector @var{x}.
##
## @code{@var{pd} = fitdist (@var{x}, @var{distname}, @var{Name}, @var{Value})}
## creates the probability distribution object with additional options specified
## by one or more @qcode{Name-Value} pair arguments listed below.
##
## @multitable @columnfractions 0.18 0.02 0.8
## @headitem @var{Name} @tab @tab @var{Value}
##
## @item @qcode{"distribution"} @tab @tab A character vector specifying the
## distribution type for which to estimate parameters.
##
## @item @qcode{"Ntrials"} @tab @tab A scalar specifying the number of trials
## for the corresponding element of @var{x} for the binomial distribution.
##
## @item @qcode{"theta"} @tab @tab A scalar specifying the location parameter
## for the generalized Pareto distribution.
##
## @item @qcode{"mu"} @tab @tab A scalar specifying the location parameter
## for the half-normal distribution.
##
## @item @qcode{"censoring"} @tab @tab A vector of the same size as @var{x}
## indicating censored data in @var{x}.  By default it is
## @qcode{@var{censor} = zeros (size (@var{x}))}.
##
## @item @qcode{"frequency"} @tab @tab A vector of nonnegative integer counts of
## the same size as @var{x} used as frequency observations.  By default it is
## @qcode{@var{freq} = ones (size (@var{x}))}.
##
## @item @qcode{"alpha"} @tab @tab A scalar in the range @math{(0,1)}, as the
## significance level for the confidence interval @var{pci}.  By default it is
## 0.05 corresponding to 95% confidence intervals.
##
## @item @qcode{"options"} @tab @tab A structure specifying the control
## parameters for the iterative algorithm used to compute ML estimates with the
## @code{fminsearch} function.
## @end multitable
##
## @code{[@var{pdca}, @var{gn}, @var{gl}] = fitdist (@var{x}, @var{distname},
## @qcode{"By"}, @var{groupvar})} creates probability distribution objects by
## fitting the distribution specified by @var{distname} to the data in @var{x}
## based on the grouping variable @var{groupvar}. It returns a cell array of
## fitted probability distribution object, @var{pdca}, a cell array of group
## labels, @var{gn}, and a cell array of grouping variable levels, @var{gl}.
##
## @code{[@var{pdca}, @var{gn}, @var{gl}] = fitdist (@var{x}, @var{distname},
## @qcode{"By"}, @var{groupvar}, @var{Name}, @var{Value})} returns the same
## output arguments using additional options specified by one or more
## @qcode{Name-Value} pair arguments mentioned above.
##
## Note: calling @code{fitdist} without any input arguments will return a cell
## array of character vectors listing all supported distributions.
##
## @seealso{makedist}
## @end deftypefn

function [varargout] = fitdist (varargin)

  ## Add list of supported probability distribution objects
  PDO = {'Beta'; 'Binomial'; 'BirnbaumSaunders'; 'Burr'; 'Exponential'; ...
         'ExtremeValue'; 'Gamma'; 'GeneralizedExtremeValue'; ...
         'GeneralizedPareto'; 'HalfNormal'; 'InverseGaussian'; ...
         'Kernel'; 'Logistic'; 'Loglogistic'; 'Lognormal'; 'Nakagami'; ...
         'NegativeBinomial'; 'Normal'; 'Poisson'; 'Rayleigh'; 'Rician'; ...
         'Stable'; 'tLocationScale'; 'Weibull'};

  ABBR = {"bisa", "ev", "gev", "gp", "hn", "invg", "nbin", "tls"};

  ## Check for input arguments
  if (nargin == 0)
    varargout{1} = PDO;
    return
  elseif (nargin == 1)
    error ("fitdist: DISTNAME is required.");
  else
    x = varargin{1};
    distname = varargin{2};
    varargin([1:2]) = [];
  endif

  ## Check distribution name
  if (! (ischar (distname) && size (distname, 1) == 1))
    error ("fitdist: DISTNAME must be a character vector.");
  elseif (! (any (strcmpi (distname, PDO)) || any (strcmpi (distname, ABBR))))
    error ("fitdist: unrecognized distribution name.");
  endif

  ## Check data in X being a real vector
  if (! (isvector (x) && isnumeric (x) && isreal (x)))
    error ("fitdist: X must be a numeric vector of real values.");
  endif

  ## Add defaults
  groupvar = [];
  censor = zeros (size (x));
  freq = ones (size (x));
  alpha = 0.05;
  ntrials = 1;
  mu = 0;
  theta = 1;
  options.Display = "off";
  options.MaxFunEvals = 400;
  options.MaxIter = 200;
  options.TolX = 1e-6;

  ## Parse extra arguments
  if (mod (numel (varargin), 2) != 0)
    error ("fitdist: optional arguments must be in NAME-VALUE pairs.");
  endif
  while (numel (varargin) > 0)
    switch (tolower (varargin{1}))
      case "by"
        groupvar = varargin{2};
        if (! isequal (size (x), size (groupvar)) && ! isempty (groupvar))
          error (strcat (["fitdist: GROUPVAR argument must have the same"], ...
                         [" size as the input data in X."]));
        endif
      case "censoring"
        censor = varargin{2};
        if (! isequal (size (x), size (censor)))
          error (strcat (["fitdist: 'censoring' argument must have the"], ...
                         [" same size as the input data in X."]));
        endif
      case "frequency"
        freq = varargin{2};
        if (! isequal (size (x), size (freq)))
          error (strcat (["fitdist: 'frequency' argument must have the"], ...
                         [" same size as the input data in X."]));
        endif
        if (any (freq != round (freq)) || any (freq < 0))
          error (strcat (["fitdist: 'frequency' argument must contain"], ...
                         [" non-negative integer values."]));
        endif
      case "alpha"
        alpha = varargin{2};
        if (! isscalar (alpha) || ! isreal (alpha) || alpha <= 0 || alpha >= 1)
          error ("fitdist: invalid value for 'alpha' argument.");
        endif
      case "ntrials"
        ntrials = varargin{2};
        if (! (isscalar (ntrials) && isreal (ntrials) && ntrials > 0
                                  && fix (ntrials) == ntrials))
          error (strcat (["fitdist: 'ntrials' argument must be a positive"], ...
                         [" integer scalar value."]));
        endif
      case {"mu"}
        mu = varargin{2};
      case {"theta"}
        theta = varargin{2};
      case "options"
        options = varargin{2};
        if (! isstruct (options) || ! isfield (options, "Display") ||
            ! isfield (options, "MaxFunEvals") || ! isfield (options, "MaxIter")
                                               || ! isfield (options, "TolX"))
          error (strcat (["fitdist: 'options' argument must be a"], ...
                         [" structure compatible for 'fminsearch'."]));
        endif

      case {"kernel", "support", "width"}
        warning ("fitdist: parameter not supported yet.");
      otherwise
        error ("fitdist: unknown parameter name.");
    endswitch
    varargin([1:2]) = [];
  endwhile

  ## Handle group variable
  if (isempty (groupvar) && nargout > 1)
    error ("fitdist: must define GROUPVAR for more than one output arguments.");
  endif
  if (! isempty (groupvar))
    [g, gn, gl] = grp2idx (groupvar);
    groups = numel (gn);
  endif

  ## Switch to selected distribution
  switch (tolower (distname))

    case "beta"
      if (isempty (groupvar))
        varargout{1} = BetaDistribution.fit (x, alpha, freq, options);
      else
        pd = BetaDistribution.fit (x(g==1), alpha, freq(g==1), options);
        for i = 2:groups
          pd(i) = BetaDistribution.fit (x(g==i), alpha, freq(g==i), options);
        endfor
        varargout{1} = pd;
        varargout{2} = gn;
        varargout{3} = gl;
      endif

    case "binomial"
      if (any (x > ntrials))
        error ("fitdist: invalid NTRIALS value for Binomial distribution.")
      endif
      if (isempty (groupvar))
        varargout{1} = BinomialDistribution.fit (x, ntrials, alpha, freq);
      else
        pd = BinomialDistribution.fit (x(g==1), ntrials, alpha, freq(g==1));
        for i = 2:groups
          pd(i) = BinomialDistribution.fit (x(g==i), ntrials, alpha, freq(g==i));
        endfor
        varargout{1} = pd;
        varargout{2} = gn;
        varargout{3} = gl;
      endif

    case {"birnbaumsaunders", "bisa"}
      if (isempty (groupvar))
        varargout{1} = BirnbaumSaundersDistribution.fit ...
                       (x, alpha, censor, freq, options);
      else
        pd = BirnbaumSaundersDistribution.fit ...
             (x(g==1), alpha, censor(g==1), freq(g==1), options);
        for i = 2:groups
          pd(i) = BirnbaumSaundersDistribution.fit ...
                  (x(g==i), alpha, censor(g==i), freq(g==i), options);
        endfor
        varargout{1} = pd;
        varargout{2} = gn;
        varargout{3} = gl;
      endif

    case "burr"
      if (isempty (groupvar))
        varargout{1} = BurrDistribution.fit ...
                       (x, alpha, censor, freq, options);
      else
        pd = BurrDistribution.fit ...
             (x(g==1), alpha, censor(g==1), freq(g==1), options);
        for i = 2:groups
          pd(i) = BurrDistribution.fit ...
                  (x(g==i), alpha, censor(g==i), freq(g==i), options);
        endfor
        varargout{1} = pd;
        varargout{2} = gn;
        varargout{3} = gl;
      endif

    case "exponential"
      if (isempty (groupvar))
        varargout{1} = ExponentialDistribution.fit (x, alpha, censor, freq);
      else
        pd = ExponentialDistribution.fit ...
             (x(g==1), alpha, censor(g==1), freq(g==1));
        for i = 2:groups
          pd(i) = ExponentialDistribution.fit ...
                  (x(g==i), alpha, censor(g==i), freq(g==i));
        endfor
        varargout{1} = pd;
        varargout{2} = gn;
        varargout{3} = gl;
      endif

    case {"extremevalue", "ev"}
      if (isempty (groupvar))
        varargout{1} = ExtremeValueDistribution.fit ...
                       (x, alpha, censor, freq, options);
      else
        pd = ExtremeValueDistribution.fit ...
             (x(g==1), alpha, censor(g==1), freq(g==1), options);
        for i = 2:groups
          pd(i) = ExtremeValueDistribution.fit ...
                  (x(g==i), alpha, censor(g==i), freq(g==i), options);
        endfor
        varargout{1} = pd;
        varargout{2} = gn;
        varargout{3} = gl;
      endif

    case "gamma"
      if (isempty (groupvar))
        varargout{1} = GammaDistribution.fit (x, alpha, censor, freq, options);
      else
        pd = GammaDistribution.fit ...
             (x(g==1), alpha, censor(g==1), freq(g==1), options);
        for i = 2:groups
          pd(i) = GammaDistribution.fit ...
                  (x(g==i), alpha, censor(g==i), freq(g==i), options);
        endfor
        varargout{1} = pd;
        varargout{2} = gn;
        varargout{3} = gl;
      endif

    case {"generalizedextremevalue", "gev"}
      if (isempty (groupvar))
        varargout{1} = GeneralizedExtremeValueDistribution.fit ...
                       (x, alpha, freq, options);
      else
        pd = GeneralizedExtremeValueDistribution.fit ...
             (x(g==1), alpha, freq(g==1), options);
        for i = 2:groups
          pd(i) = GeneralizedExtremeValueDistribution.fit ...
                  (x(g==i), alpha, freq(g==i), options);
        endfor
        varargout{1} = pd;
        varargout{2} = gn;
        varargout{3} = gl;
      endif

    case {"generalizedpareto", "gp"}
      if (any (x - theta < 0))
        error (strcat (["fitdist: invalid THETA value for generalized"], ...
                       [" Pareto distribution."]));
      endif
      if (isempty (groupvar))
        varargout{1} = GeneralizedParetoDistribution.fit ...
                       (x, theta, alpha, freq, options);
      else
        pd = GeneralizedParetoDistribution.fit ...
             (x(g==1), theta, alpha, freq(g==1), options);
        for i = 2:groups
          pd(i) = GeneralizedParetoDistribution.fit ...
                  (x(g==i), theta, alpha, freq(g==i), options);
        endfor
        varargout{1} = pd;
        varargout{2} = gn;
        varargout{3} = gl;
      endif

    case {"halfnormal", "hn"}
      if (any (x - mu < 0))
        error ("fitdist: invalid MU value for half-normal distribution.");
      endif
      if (isempty (groupvar))
        varargout{1} = HalfNormalDistribution.fit (x, mu, alpha, freq);
      else
        pd = HalfNormalDistribution.fit (x(g==1), mu, alpha, freq(g==1));
        for i = 2:groups
          pd(i) = HalfNormalDistribution.fit (x(g==i), mu, alpha, freq(g==i));
        endfor
        varargout{1} = pd;
        varargout{2} = gn;
        varargout{3} = gl;
      endif

    case {"inversegaussian", "invg"}
      if (isempty (groupvar))
        varargout{1} = InverseGaussianDistribution.fit ...
                       (x, alpha, censor, freq, options);
      else
        pd = InverseGaussianDistribution.fit ...
             (x(g==1), alpha, censor(g==1), freq(g==1), options);
        for i = 2:groups
          pd(i) = InverseGaussianDistribution.fit ...
                  (x(g==i), alpha, censor(g==i), freq(g==i), options);
        endfor
        varargout{1} = pd;
        varargout{2} = gn;
        varargout{3} = gl;
      endif

    case "kernel"
      warning ("fitdist: 'Kernel' distribution not supported yet.");
      if (isempty (groupvar))
        varargout{1} = [];
      else
        varargout{1} = [];
        varargout{2} = gn;
        varargout{3} = gl;
      endif

    case "logistic"
      if (isempty (groupvar))
        varargout{1} = LogisticDistribution.fit ...
                       (x, alpha, censor, freq, options);
      else
        pd = LogisticDistribution.fit ...
             (x(g==1), alpha, censor(g==1), freq(g==1), options);
        for i = 2:groups
          pd(i) = LogisticDistribution.fit ...
                  (x(g==i), alpha, censor(g==i), freq(g==i), options);
        endfor
        varargout{1} = pd;
        varargout{2} = gn;
        varargout{3} = gl;
      endif

    case "loglogistic"
      if (isempty (groupvar))
        varargout{1} = LoglogisticDistribution.fit ...
                       (x, alpha, censor, freq, options);
      else
        pd = LoglogisticDistribution.fit ...
             (x(g==1), alpha, censor(g==1), freq(g==1), options);
        for i = 2:groups
          pd(i) = LoglogisticDistribution.fit ...
                  (x(g==i), alpha, censor(g==i), freq(g==i), options);
        endfor
        varargout{1} = pd;
        varargout{2} = gn;
        varargout{3} = gl;
      endif

    case "lognormal"
      if (isempty (groupvar))
        varargout{1} = LognormalDistribution.fit ...
                       (x, alpha, censor, freq, options);
      else
        pd = LognormalDistribution.fit ...
             (x(g==1), alpha, censor(g==1), freq(g==1), options);
        for i = 2:groups
          pd(i) = LognormalDistribution.fit ...
                  (x(g==i), alpha, censor(g==i), freq(g==i), options);
        endfor
        varargout{1} = pd;
        varargout{2} = gn;
        varargout{3} = gl;
      endif

    case "nakagami"
      if (isempty (groupvar))
        varargout{1} = NakagamiDistribution.fit ...
                       (x, alpha, censor, freq, options);
      else
        pd = NakagamiDistribution.fit ...
             (x(g==1), alpha, censor(g==1), freq(g==1), options);
        for i = 2:groups
          pd(i) = NakagamiDistribution.fit ...
                  (x(g==i), alpha, censor(g==i), freq(g==i), options);
        endfor
        varargout{1} = pd;
        varargout{2} = gn;
        varargout{3} = gl;
      endif

    case {"negativebinomial", "nbin"}
      if (isempty (groupvar))
        varargout{1} = NegativeBinomialDistribution.fit ...
                       (x, alpha, freq, options);
      else
        pd = NegativeBinomialDistribution.fit ...
             (x(g==1), alpha, freq(g==1), options);
        for i = 2:groups
          pd(i) = NegativeBinomialDistribution.fit ...
                  (x(g==i), alpha, freq(g==i), options);
        endfor
        varargout{1} = pd;
        varargout{2} = gn;
        varargout{3} = gl;
      endif

    case "normal"
      if (isempty (groupvar))
        varargout{1} = NormalDistribution.fit (x, alpha, censor, freq, options);
      else
        pd = NormalDistribution.fit ...
             (x(g==1), alpha, censor(g==1), freq(g==1), options);
        for i = 2:groups
          pd(i) = NormalDistribution.fit ...
                  (x(g==i), alpha, censor(g==i), freq(g==i), options);
        endfor
        varargout{1} = pd;
        varargout{2} = gn;
        varargout{3} = gl;
      endif

    case "poisson"
      if (isempty (groupvar))
        varargout{1} = PoissonDistribution.fit (x, alpha, freq);
      else
        pd = PoissonDistribution.fit (x(g==1), alpha, freq(g==1));
        for i = 2:groups
          pd(i) = PoissonDistribution.fit (x(g==i), alpha, freq(g==i));
        endfor
        varargout{1} = pd;
        varargout{2} = gn;
        varargout{3} = gl;
      endif

    case "rayleigh"
      if (isempty (groupvar))
        varargout{1} = RayleighDistribution.fit (x, alpha, censor, freq);
      else
        pd = RayleighDistribution.fit (x(g==1), alpha, censor(g==1), freq(g==1));
        for i = 2:groups
          pd(i) = RayleighDistribution.fit ...
                  (x(g==i), alpha, censor(g==i), freq(g==i));
        endfor
        varargout{1} = pd;
        varargout{2} = gn;
        varargout{3} = gl;
      endif

    case "rician"
      if (isempty (groupvar))
        varargout{1} = RicianDistribution.fit (x, alpha, censor, freq, options);
      else
        pd = RicianDistribution.fit ...
             (x(g==1), alpha, censor(g==1), freq(g==1), options);
        for i = 2:groups
          pd(i) = RicianDistribution.fit ...
                  (x(g==i), alpha, censor(g==i), freq(g==i), options);
        endfor
        varargout{1} = pd;
        varargout{2} = gn;
        varargout{3} = gl;
      endif

    case "stable"
      warning ("fitdist: 'Stable' distribution not supported yet.");
      if (isempty (groupvar))
        varargout{1} = [];
      else
        varargout{1} = [];
        varargout{2} = gn;
        varargout{3} = gl;
      endif

    case {"tlocationscale", "tls"}
      if (isempty (groupvar))
        varargout{1} = tLocationScaleDistribution.fit ...
                       (x, alpha, censor, freq, options);
      else
        pd = tLocationScaleDistribution.fit ...
             (x(g==1), alpha, censor(g==1), freq(g==1), options);
        for i = 2:groups
          pd(i) = tLocationScaleDistribution.fit ...
                  (x(g==i), alpha, censor(g==i), freq(g==i), options);
        endfor
        varargout{1} = pd;
        varargout{2} = gn;
        varargout{3} = gl;
      endif

    case "weibull"
      if (isempty (groupvar))
        varargout{1} = WeibullDistribution.fit (x, alpha, censor, freq, options);
      else
        pd = WeibullDistribution.fit ...
             (x(g==1), alpha, censor(g==1), freq(g==1), options);
        for i = 2:groups
          pd(i) = WeibullDistribution.fit ...
                  (x(g==i), alpha, censor(g==i), freq(g==i), options);
        endfor
        varargout{1} = pd;
        varargout{2} = gn;
        varargout{3} = gl;
      endif

  endswitch

endfunction

## Test output
%!test
%! x = betarnd (1, 1, 100, 1);
%! pd = fitdist (x, "Beta");
%! [phat, pci] = betafit (x);
%! assert ([pd.a, pd.b], phat);
%! assert (paramci (pd), pci);
%!test
%! x1 = betarnd (1, 1, 100, 1);
%! x2 = betarnd (5, 2, 100, 1);
%! pd = fitdist ([x1; x2], "Beta", "By", [ones(100,1); 2*ones(100,1)]);
%! [phat, pci] = betafit (x1);
%! assert ([pd(1).a, pd(1).b], phat);
%! assert (paramci (pd(1)), pci);
%! [phat, pci] = betafit (x2);
%! assert ([pd(2).a, pd(2).b], phat);
%! assert (paramci (pd(2)), pci);
%!test
%! N = 1;
%! x = binornd (N, 0.5, 100, 1);
%! pd = fitdist (x, "binomial");
%! [phat, pci] = binofit (sum (x), numel (x));
%! assert ([pd.N, pd.p], [N, phat]);
%! assert (paramci (pd), pci);
%!test
%! N = 3;
%! x = binornd (N, 0.4, 100, 1);
%! pd = fitdist (x, "binomial", "ntrials", N);
%! [phat, pci] = binofit (sum (x), numel (x) * N);
%! assert ([pd.N, pd.p], [N, phat]);
%! assert (paramci (pd), pci);
%!test
%! N = 1;
%! x1 = binornd (N, 0.5, 100, 1);
%! x2 = binornd (N, 0.7, 100, 1);
%! pd = fitdist ([x1; x2], "binomial", "By", [ones(100,1); 2*ones(100,1)]);
%! [phat, pci] = binofit (sum (x1), numel (x1));
%! assert ([pd(1).N, pd(1).p], [N, phat]);
%! assert (paramci (pd(1)), pci);
%! [phat, pci] = binofit (sum (x2), numel (x2));
%! assert ([pd(2).N, pd(2).p], [N, phat]);
%! assert (paramci (pd(2)), pci);
%!test
%! N = 5;
%! x1 = binornd (N, 0.5, 100, 1);
%! x2 = binornd (N, 0.8, 100, 1);
%! pd = fitdist ([x1; x2], "binomial", "ntrials", N, ...
%!               "By", [ones(100,1); 2*ones(100,1)]);
%! [phat, pci] = binofit (sum (x1), numel (x1) * N);
%! assert ([pd(1).N, pd(1).p], [N, phat]);
%! assert (paramci (pd(1)), pci);
%! [phat, pci] = binofit (sum (x2), numel (x2) * N);
%! assert ([pd(2).N, pd(2).p], [N, phat]);
%! assert (paramci (pd(2)), pci);
%!test
%! x = bisarnd (1, 1, 100, 1);
%! pd = fitdist (x, "BirnbaumSaunders");
%! [phat, pci] = bisafit (x);
%! assert ([pd.beta, pd.gamma], phat);
%! assert (paramci (pd), pci);
%!test
%! x1 = bisarnd (1, 1, 100, 1);
%! x2 = bisarnd (5, 2, 100, 1);
%! pd = fitdist ([x1; x2], "bisa", "By", [ones(100,1); 2*ones(100,1)]);
%! [phat, pci] = bisafit (x1);
%! assert ([pd(1).beta, pd(1).gamma], phat);
%! assert (paramci (pd(1)), pci);
%! [phat, pci] = bisafit (x2);
%! assert ([pd(2).beta, pd(2).gamma], phat);
%! assert (paramci (pd(2)), pci);
%!test
%! x = burrrnd (1, 2, 1, 100, 1);
%! pd = fitdist (x, "Burr");
%! [phat, pci] = burrfit (x);
%! assert ([pd.alpha, pd.c, pd.k], phat);
%! assert (paramci (pd), pci);
%!test
%! rand ("seed", 4);   # for reproducibility
%! x1 = burrrnd (1, 2, 1, 100, 1);
%! rand ("seed", 3);   # for reproducibility
%! x2 = burrrnd (1, 0.5, 2, 100, 1);
%! pd = fitdist ([x1; x2], "burr", "By", [ones(100,1); 2*ones(100,1)]);
%! [phat, pci] = burrfit (x1);
%! assert ([pd(1).alpha, pd(1).c, pd(1).k], phat);
%! assert (paramci (pd(1)), pci);
%! [phat, pci] = burrfit (x2);
%! assert ([pd(2).alpha, pd(2).c, pd(2).k], phat);
%! assert (paramci (pd(2)), pci);
%!test
%! x = exprnd (1, 100, 1);
%! pd = fitdist (x, "exponential");
%! [muhat, muci] = expfit (x);
%! assert ([pd.mu], muhat);
%! assert (paramci (pd), muci);
%!test
%! x1 = exprnd (1, 100, 1);
%! x2 = exprnd (5, 100, 1);
%! pd = fitdist ([x1; x2], "exponential", "By", [ones(100,1); 2*ones(100,1)]);
%! [muhat, muci] = expfit (x1);
%! assert ([pd(1).mu], muhat);
%! assert (paramci (pd(1)), muci);
%! [muhat, muci] = expfit (x2);
%! assert ([pd(2).mu], muhat);
%! assert (paramci (pd(2)), muci);
%!test
%! x = evrnd (1, 1, 100, 1);
%! pd = fitdist (x, "ev");
%! [phat, pci] = evfit (x);
%! assert ([pd.mu, pd.sigma], phat);
%! assert (paramci (pd), pci);
%!test
%! x1 = evrnd (1, 1, 100, 1);
%! x2 = evrnd (5, 2, 100, 1);
%! pd = fitdist ([x1; x2], "extremevalue", "By", [ones(100,1); 2*ones(100,1)]);
%! [phat, pci] = evfit (x1);
%! assert ([pd(1).mu, pd(1).sigma], phat);
%! assert (paramci (pd(1)), pci);
%! [phat, pci] = evfit (x2);
%! assert ([pd(2).mu, pd(2).sigma], phat);
%! assert (paramci (pd(2)), pci);
%!test
%! x = gamrnd (1, 1, 100, 1);
%! pd = fitdist (x, "Gamma");
%! [phat, pci] = gamfit (x);
%! assert ([pd.a, pd.b], phat);
%! assert (paramci (pd), pci);
%!test
%! x1 = gamrnd (1, 1, 100, 1);
%! x2 = gamrnd (5, 2, 100, 1);
%! pd = fitdist ([x1; x2], "Gamma", "By", [ones(100,1); 2*ones(100,1)]);
%! [phat, pci] = gamfit (x1);
%! assert ([pd(1).a, pd(1).b], phat);
%! assert (paramci (pd(1)), pci);
%! [phat, pci] = gamfit (x2);
%! assert ([pd(2).a, pd(2).b], phat);
%! assert (paramci (pd(2)), pci);
%!test
%! rand ("seed", 4);   # for reproducibility
%! x = gevrnd (-0.5, 1, 2, 1000, 1);
%! pd = fitdist (x, "generalizedextremevalue");
%! [phat, pci] = gevfit (x);
%! assert ([pd.k, pd.sigma, pd.mu], phat);
%! assert (paramci (pd), pci);
%!test
%! rand ("seed", 5);   # for reproducibility
%! x1 = gevrnd (-0.5, 1, 2, 1000, 1);
%! rand ("seed", 9);   # for reproducibility
%! x2 = gevrnd (0, 1, -4, 1000, 1);
%! pd = fitdist ([x1; x2], "gev", "By", [ones(1000,1); 2*ones(1000,1)]);
%! [phat, pci] = gevfit (x1);
%! assert ([pd(1).k, pd(1).sigma, pd(1).mu], phat);
%! assert (paramci (pd(1)), pci);
%! [phat, pci] = gevfit (x2);
%! assert ([pd(2).k, pd(2).sigma, pd(2).mu], phat);
%! assert (paramci (pd(2)), pci);
%!test
%! x = gprnd (1, 1, 1, 100, 1);
%! pd = fitdist (x, "GeneralizedPareto");
%! [phat, pci] = gpfit (x, 1);
%! assert ([pd.k, pd.sigma, pd.theta], phat);
%! assert (paramci (pd), pci);
%!test
%! x = gprnd (1, 1, 2, 100, 1);
%! pd = fitdist (x, "GeneralizedPareto", "theta", 2);
%! [phat, pci] = gpfit (x, 2);
%! assert ([pd.k, pd.sigma, pd.theta], phat);
%! assert (paramci (pd), pci);
%!test
%! x1 = gprnd (1, 1, 1, 100, 1);
%! x2 = gprnd (0, 2, 1, 100, 1);
%! pd = fitdist ([x1; x2], "gp", "By", [ones(100,1); 2*ones(100,1)]);
%! [phat, pci] = gpfit (x1, 1);
%! assert ([pd(1).k, pd(1).sigma, pd(1).theta], phat);
%! assert (paramci (pd(1)), pci);
%! [phat, pci] = gpfit (x2, 1);
%! assert ([pd(2).k, pd(2).sigma, pd(2).theta], phat);
%! assert (paramci (pd(2)), pci);
%!test
%! x1 = gprnd (3, 2, 2, 100, 1);
%! x2 = gprnd (2, 3, 2, 100, 1);
%! pd = fitdist ([x1; x2], "GeneralizedPareto", "theta", 2, ...
%!               "By", [ones(100,1); 2*ones(100,1)]);
%! [phat, pci] = gpfit (x1, 2);
%! assert ([pd(1).k, pd(1).sigma, pd(1).theta], phat);
%! assert (paramci (pd(1)), pci);
%! [phat, pci] = gpfit (x2, 2);
%! assert ([pd(2).k, pd(2).sigma, pd(2).theta], phat);
%! assert (paramci (pd(2)), pci);
%!test
%! x = hnrnd (0, 1, 100, 1);
%! pd = fitdist (x, "HalfNormal");
%! [phat, pci] = hnfit (x, 0);
%! assert ([pd.mu, pd.sigma], phat);
%! assert (paramci (pd), pci);
%!test
%! x = hnrnd (1, 1, 100, 1);
%! pd = fitdist (x, "HalfNormal", "mu", 1);
%! [phat, pci] = hnfit (x, 1);
%! assert ([pd.mu, pd.sigma], phat);
%! assert (paramci (pd), pci);
%!test
%! x1 = hnrnd (0, 1, 100, 1);
%! x2 = hnrnd (0, 2, 100, 1);
%! pd = fitdist ([x1; x2], "HalfNormal", "By", [ones(100,1); 2*ones(100,1)]);
%! [phat, pci] = hnfit (x1, 0);
%! assert ([pd(1).mu, pd(1).sigma], phat);
%! assert (paramci (pd(1)), pci);
%! [phat, pci] = hnfit (x2, 0);
%! assert ([pd(2).mu, pd(2).sigma], phat);
%! assert (paramci (pd(2)), pci);
%!test
%! x1 = hnrnd (2, 1, 100, 1);
%! x2 = hnrnd (2, 2, 100, 1);
%! pd = fitdist ([x1; x2], "HalfNormal", "mu", 2, ...
%!               "By", [ones(100,1); 2*ones(100,1)]);
%! [phat, pci] = hnfit (x1, 2);
%! assert ([pd(1).mu, pd(1).sigma], phat);
%! assert (paramci (pd(1)), pci);
%! [phat, pci] = hnfit (x2, 2);
%! assert ([pd(2).mu, pd(2).sigma], phat);
%! assert (paramci (pd(2)), pci);
%!test
%! x = invgrnd (1, 1, 100, 1);
%! pd = fitdist (x, "InverseGaussian");
%! [phat, pci] = invgfit (x);
%! assert ([pd.mu, pd.lambda], phat);
%! assert (paramci (pd), pci);
%!test
%! x1 = invgrnd (1, 1, 100, 1);
%! x2 = invgrnd (5, 2, 100, 1);
%! pd = fitdist ([x1; x2], "InverseGaussian", "By", [ones(100,1); 2*ones(100,1)]);
%! [phat, pci] = invgfit (x1);
%! assert ([pd(1).mu, pd(1).lambda], phat);
%! assert (paramci (pd(1)), pci);
%! [phat, pci] = invgfit (x2);
%! assert ([pd(2).mu, pd(2).lambda], phat);
%! assert (paramci (pd(2)), pci);
%!test
%! x = logirnd (1, 1, 100, 1);
%! pd = fitdist (x, "logistic");
%! [phat, pci] = logifit (x);
%! assert ([pd.mu, pd.sigma], phat);
%! assert (paramci (pd), pci);
%!test
%! x1 = logirnd (1, 1, 100, 1);
%! x2 = logirnd (5, 2, 100, 1);
%! pd = fitdist ([x1; x2], "logistic", "By", [ones(100,1); 2*ones(100,1)]);
%! [phat, pci] = logifit (x1);
%! assert ([pd(1).mu, pd(1).sigma], phat);
%! assert (paramci (pd(1)), pci);
%! [phat, pci] = logifit (x2);
%! assert ([pd(2).mu, pd(2).sigma], phat);
%! assert (paramci (pd(2)), pci);
%!test
%! x = loglrnd (1, 1, 100, 1);
%! pd = fitdist (x, "loglogistic");
%! [phat, pci] = loglfit (x);
%! assert ([pd.mu, pd.sigma], phat);
%! assert (paramci (pd), pci);
%!test
%! x1 = loglrnd (1, 1, 100, 1);
%! x2 = loglrnd (5, 2, 100, 1);
%! pd = fitdist ([x1; x2], "loglogistic", "By", [ones(100,1); 2*ones(100,1)]);
%! [phat, pci] = loglfit (x1);
%! assert ([pd(1).mu, pd(1).sigma], phat);
%! assert (paramci (pd(1)), pci);
%! [phat, pci] = loglfit (x2);
%! assert ([pd(2).mu, pd(2).sigma], phat);
%! assert (paramci (pd(2)), pci);
%!test
%! x = lognrnd (1, 1, 100, 1);
%! pd = fitdist (x, "lognormal");
%! [phat, pci] = lognfit (x);
%! assert ([pd.mu, pd.sigma], phat);
%! assert (paramci (pd), pci);
%!test
%! x1 = lognrnd (1, 1, 100, 1);
%! x2 = lognrnd (5, 2, 100, 1);
%! pd = fitdist ([x1; x2], "lognormal", "By", [ones(100,1); 2*ones(100,1)]);
%! [phat, pci] = lognfit (x1);
%! assert ([pd(1).mu, pd(1).sigma], phat);
%! assert (paramci (pd(1)), pci);
%! [phat, pci] = lognfit (x2);
%! assert ([pd(2).mu, pd(2).sigma], phat);
%! assert (paramci (pd(2)), pci);
%!test
%! x = nakarnd (2, 0.5, 100, 1);
%! pd = fitdist (x, "Nakagami");
%! [phat, pci] = nakafit (x);
%! assert ([pd.mu, pd.omega], phat);
%! assert (paramci (pd), pci);
%!test
%! x1 = nakarnd (2, 0.5, 100, 1);
%! x2 = nakarnd (5, 0.8, 100, 1);
%! pd = fitdist ([x1; x2], "Nakagami", "By", [ones(100,1); 2*ones(100,1)]);
%! [phat, pci] = nakafit (x1);
%! assert ([pd(1).mu, pd(1).omega], phat);
%! assert (paramci (pd(1)), pci);
%! [phat, pci] = nakafit (x2);
%! assert ([pd(2).mu, pd(2).omega], phat);
%! assert (paramci (pd(2)), pci);
%!test
%! randp ("seed", 123);
%! randg ("seed", 321);
%! x = nbinrnd (2, 0.5, 100, 1);
%! pd = fitdist (x, "negativebinomial");
%! [phat, pci] = nbinfit (x);
%! assert ([pd.R, pd.P], phat);
%! assert (paramci (pd), pci);
%!test
%! randp ("seed", 345);
%! randg ("seed", 543);
%! x1 = nbinrnd (2, 0.5, 100, 1);
%! randp ("seed", 432);
%! randg ("seed", 234);
%! x2 = nbinrnd (5, 0.8, 100, 1);
%! pd = fitdist ([x1; x2], "nbin", "By", [ones(100,1); 2*ones(100,1)]);
%! [phat, pci] = nbinfit (x1);
%! assert ([pd(1).R, pd(1).P], phat);
%! assert (paramci (pd(1)), pci);
%! [phat, pci] = nbinfit (x2);
%! assert ([pd(2).R, pd(2).P], phat);
%! assert (paramci (pd(2)), pci);
%!test
%! x = normrnd (1, 1, 100, 1);
%! pd = fitdist (x, "normal");
%! [muhat, sigmahat, muci, sigmaci] = normfit (x);
%! assert ([pd.mu, pd.sigma], [muhat, sigmahat]);
%! assert (paramci (pd), [muci, sigmaci]);
%!test
%! x1 = normrnd (1, 1, 100, 1);
%! x2 = normrnd (5, 2, 100, 1);
%! pd = fitdist ([x1; x2], "normal", "By", [ones(100,1); 2*ones(100,1)]);
%! [muhat, sigmahat, muci, sigmaci] = normfit (x1);
%! assert ([pd(1).mu, pd(1).sigma], [muhat, sigmahat]);
%! assert (paramci (pd(1)), [muci, sigmaci]);
%! [muhat, sigmahat, muci, sigmaci] = normfit (x2);
%! assert ([pd(2).mu, pd(2).sigma], [muhat, sigmahat]);
%! assert (paramci (pd(2)), [muci, sigmaci]);
%!test
%! x = poissrnd (1, 100, 1);
%! pd = fitdist (x, "poisson");
%! [phat, pci] = poissfit (x);
%! assert (pd.lambda, phat);
%! assert (paramci (pd), pci);
%!test
%! x1 = poissrnd (1, 100, 1);
%! x2 = poissrnd (5, 100, 1);
%! pd = fitdist ([x1; x2], "poisson", "By", [ones(100,1); 2*ones(100,1)]);
%! [phat, pci] = poissfit (x1);
%! assert (pd(1).lambda, phat);
%! assert (paramci (pd(1)), pci);
%! [phat, pci] = poissfit (x2);
%! assert (pd(2).lambda, phat);
%! assert (paramci (pd(2)), pci);
%!test
%! x = raylrnd (1, 100, 1);
%! pd = fitdist (x, "rayleigh");
%! [phat, pci] = raylfit (x);
%! assert (pd.sigma, phat);
%! assert (paramci (pd), pci);
%!test
%! x1 = raylrnd (1, 100, 1);
%! x2 = raylrnd (5, 100, 1);
%! pd = fitdist ([x1; x2], "rayleigh", "By", [ones(100,1); 2*ones(100,1)]);
%! [phat, pci] = raylfit (x1);
%! assert ( pd(1).sigma, phat);
%! assert (paramci (pd(1)), pci);
%! [phat, pci] = raylfit (x2);
%! assert (pd(2).sigma, phat);
%! assert (paramci (pd(2)), pci);
%!test
%! x = ricernd (1, 1, 100, 1);
%! pd = fitdist (x, "rician");
%! [phat, pci] = ricefit (x);
%! assert ([pd.s, pd.sigma], phat);
%! assert (paramci (pd), pci);
%!test
%! x1 = ricernd (1, 1, 100, 1);
%! x2 = ricernd (5, 2, 100, 1);
%! pd = fitdist ([x1; x2], "rician", "By", [ones(100,1); 2*ones(100,1)]);
%! [phat, pci] = ricefit (x1);
%! assert ([pd(1).s, pd(1).sigma], phat);
%! assert (paramci (pd(1)), pci);
%! [phat, pci] = ricefit (x2);
%! assert ([pd(2).s, pd(2).sigma], phat);
%! assert (paramci (pd(2)), pci);
%!warning <fitdist: 'Stable' distribution not supported yet.> ...
%! fitdist ([1 2 3 4 5], "Stable");
%!test
%! x = tlsrnd (0, 1, 1, 100, 1);
%! pd = fitdist (x, "tlocationscale");
%! [phat, pci] = tlsfit (x);
%! assert ([pd.mu, pd.sigma, pd.nu], phat);
%! assert (paramci (pd), pci);
%!test
%! x1 = tlsrnd (0, 1, 1, 100, 1);
%! x2 = tlsrnd (5, 2, 1, 100, 1);
%! pd = fitdist ([x1; x2], "tlocationscale", "By", [ones(100,1); 2*ones(100,1)]);
%! [phat, pci] = tlsfit (x1);
%! assert ([pd(1).mu, pd(1).sigma, pd(1).nu], phat);
%! assert (paramci (pd(1)), pci);
%! [phat, pci] = tlsfit (x2);
%! assert ([pd(2).mu, pd(2).sigma, pd(2).nu], phat);
%! assert (paramci (pd(2)), pci);
%!test
%! x = [1 2 3 4 5];
%! pd = fitdist (x, "weibull");
%! [phat, pci] = wblfit (x);
%! assert ([pd.lambda, pd.k], phat);
%! assert (paramci (pd), pci);
%!test
%! x = [1 2 3 4 5 6 7 8 9 10];
%! pd = fitdist (x, "weibull", "By", [1 1 1 1 1 2 2 2 2 2]);
%! [phat, pci] = wblfit (x(1:5));
%! assert ([pd(1).lambda, pd(1).k], phat);
%! assert (paramci (pd(1)), pci);
%! [phat, pci] = wblfit (x(6:10));
%! assert ([pd(2).lambda, pd(2).k], phat);
%! assert (paramci (pd(2)), pci);

## Test input validation
%!error <fitdist: DISTNAME is required.> fitdist (1)
%!error <fitdist: DISTNAME must be a character vector.> fitdist (1, ["as";"sd"])
%!error <fitdist: unrecognized distribution name.> fitdist (1, "some")
%!error <fitdist: X must be a numeric vector of real values.> ...
%! fitdist (ones (2), "normal")
%!error <fitdist: X must be a numeric vector of real values.> ...
%! fitdist ([i, 2, 3], "normal")
%!error <fitdist: X must be a numeric vector of real values.> ...
%! fitdist (["a", "s", "d"], "normal")
%!error <fitdist: optional arguments must be in NAME-VALUE pairs.> ...
%! fitdist ([1, 2, 3], "normal", "By")
%!error <fitdist: GROUPVAR argument must have the same size as the input data in X.> ...
%! fitdist ([1, 2, 3], "normal", "By", [1, 2])
%!error <fitdist: 'censoring' argument must have the same size as the input data in X.> ...
%! fitdist ([1, 2, 3], "normal", "Censoring", [1, 2])
%!error <fitdist: 'frequency' argument must have the same size as the input data in X.> ...
%! fitdist ([1, 2, 3], "normal", "frequency", [1, 2])
%!error <fitdist: 'frequency' argument must contain non-negative integer values.> ...
%! fitdist ([1, 2, 3], "negativebinomial", "frequency", [1, -2, 3])
%!error <fitdist: invalid value for 'alpha' argument.> ...
%! fitdist ([1, 2, 3], "normal", "alpha", [1, 2])
%!error <fitdist: invalid value for 'alpha' argument.> ...
%! fitdist ([1, 2, 3], "normal", "alpha", i)
%!error <fitdist: invalid value for 'alpha' argument.> ...
%! fitdist ([1, 2, 3], "normal", "alpha", -0.5)
%!error <fitdist: invalid value for 'alpha' argument.> ...
%! fitdist ([1, 2, 3], "normal", "alpha", 1.5)
%!error <fitdist: 'ntrials' argument must be a positive integer scalar value.> ...
%! fitdist ([1, 2, 3], "normal", "ntrials", [1, 2])
%!error <fitdist: 'ntrials' argument must be a positive integer scalar value.> ...
%! fitdist ([1, 2, 3], "normal", "ntrials", 0)
%!error <fitdist: 'options' argument must be a structure compatible for 'fminsearch'.> ...
%! fitdist ([1, 2, 3], "normal", "options", 0)
%!error <fitdist: 'options' argument must be a structure compatible for 'fminsearch'.> ...
%! fitdist ([1, 2, 3], "normal", "options", struct ("options", 1))
%!warning fitdist ([1, 2, 3], "kernel", "kernel", "normal");
%!warning fitdist ([1, 2, 3], "kernel", "support", "positive");
%!warning fitdist ([1, 2, 3], "kernel", "width", 1);
%!error <fitdist: unknown parameter name.> ...
%! fitdist ([1, 2, 3], "normal", "param", struct ("options", 1))
%!error <fitdist: must define GROUPVAR for more than one output arguments.> ...
%! [pdca, gn, gl] = fitdist ([1, 2, 3], "normal");
%!error <fitdist: invalid THETA value for generalized Pareto distribution.> ...
%! fitdist ([1, 2, 3], "generalizedpareto", "theta", 2);
%!error <fitdist: invalid MU value for half-normal distribution.> ...
%! fitdist ([1, 2, 3], "halfnormal", "mu", 2);
