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
## @deftypefn  {statistics} {@var{pd} =} makedist (@var{distname})
## @deftypefnx {statistics} {@var{pd} =} makedist (@var{distname}, @var{Name}, @var{Value})
## @deftypefnx {statistics} {@var{list} =} makedist
##
## Create probability distribution object.
##
## @code{@var{pd} = makedist (@var{distname})} creates a probability
## distribution object for the distribution specified in @var{distname}, using
## the default parameter values.
##
## @code{@var{pd} = makedist (@var{distname}, @var{Name}, @var{Value})} also
## creates a probability distribution object with one or more distribution
## parameter values specified by @qcode{Name-Value} pair arguments.
##
## @code{@var{list} = makedist}  returns a cell array, @var{list}, containing a
## list of the probability distributions that makedist can create.
##
## @seealso{fitdist}
## @end deftypefn

function pd = makedist (varargin)

  ## Add list of supported probability distribution objects
  PDO = {'Beta'; 'Binomial'; 'BirnbaumSaunders'; 'Burr'; 'Exponential'; ...
         'ExtremeValue'; 'Gamma'; 'GeneralizedExtremeValue'; ...
         'GeneralizedPareto'; 'HalfNormal'; 'InverseGaussian'; ...
         'Logistic'; 'Loglogistic'; 'Lognormal'; 'Loguniform'; ...
         'Multinomial'; 'Nakagami'; 'NegativeBinomial'; 'Normal'; ...
         'PiecewiseLinear'; 'Poisson'; 'Rayleigh'; 'Rician'; ...
         'Stable'; 'tLocationScale'; 'Triangular'; 'Uniform'; 'Weibull'};

  ABBR = {'bisa', 'ev', 'gev', 'gp', 'hn', 'invg', 'nbin', 'tls'};

  ## Check for input arguments
  if (nargin == 0)
    pd = PDO;
    return
  else
    distname = varargin{1};
    varargin(1) = [];
  endif

  ## Check distribution name
  if (! (ischar (distname) && size (distname, 1) == 1))
    error ("makedist: DISTNAME must be a character vector.");
  elseif (! (any (strcmpi (distname, PDO)) || any (strcmpi (distname, ABBR))))
    error ("makedist: unrecognized distribution name.");
  endif

  ## Check for additional arguments being in pairs
  if (mod (numel (varargin), 2) != 0)
    error ("makedist: optional arguments must be in NAME-VALUE pairs.");
  endif

  ## Switch to selected distribution
  switch (tolower (distname))

    case 'beta'
      ## Add default parameters
      a = 1;
      b = 1;
      while (numel (varargin) > 0)
        switch (tolower (varargin{1}))
          case 'a'
            a = varargin{2};
          case 'b'
            b = varargin{2};
          otherwise
            error ("makedist: unknown parameter for 'Beta' distribution.");
        endswitch
        varargin([1:2]) = [];
      endwhile
      pd = BetaDistribution (a, b);

    case 'binomial'
      N = 1;
      p = 0.5;
      while (numel (varargin) > 0)
        switch (tolower (varargin{1}))
          case 'n'
            N = varargin{2};
          case 'p'
            p = varargin{2};
          otherwise
            error ("makedist: unknown parameter for 'Binomial' distribution.");
        endswitch
        varargin([1:2]) = [];
      endwhile
      pd = BinomialDistribution (N, p);

    case {'birnbaumsaunders', 'bisa'}
      beta = 1;
      gamma = 1;
      while (numel (varargin) > 0)
        switch (tolower (varargin{1}))
          case 'beta'
            beta = varargin{2};
          case 'gamma'
            gamma = varargin{2};
          otherwise
            error (strcat ("makedist: unknown parameter for", ...
                           " 'BirnbaumSaunders' distribution."));
        endswitch
        varargin([1:2]) = [];
      endwhile
      pd = BirnbaumSaundersDistribution (beta, gamma);

    case 'burr'
      alpha = 1;
      c = 1;
      k = 1;
      while (numel (varargin) > 0)
        switch (tolower (varargin{1}))
          case {'lambda', 'alpha'}
            alpha = varargin{2};
          case 'c'
            c = varargin{2};
          case 'k'
            k = varargin{2};
          otherwise
            error ("makedist: unknown parameter for 'Burr' distribution.");
        endswitch
        varargin([1:2]) = [];
      endwhile
      pd = BurrDistribution (alpha, c, k);

    case 'exponential'
      mu = 1;
      while (numel (varargin) > 0)
        switch (tolower (varargin{1}))
          case 'mu'
            mu = varargin{2};
          otherwise
            error (strcat ("makedist: unknown parameter for", ...
                           " 'Exponential' distribution."));
        endswitch
        varargin([1:2]) = [];
      endwhile
      pd = ExponentialDistribution (mu);

    case {'extremevalue', 'ev'}
      mu = 0;
      sigma = 1;
      while (numel (varargin) > 0)
        switch (tolower (varargin{1}))
          case 'mu'
            mu = varargin{2};
          case 'sigma'
            sigma = varargin{2};
          otherwise
            error (strcat ("makedist: unknown parameter for", ...
                           " 'ExtremeValue' distribution."));
        endswitch
        varargin([1:2]) = [];
      endwhile
      pd = ExtremeValueDistribution (mu, sigma);

    case 'gamma'
      a = 1;
      b = 1;
      while (numel (varargin) > 0)
        switch (tolower (varargin{1}))
          case 'a'
            a = varargin{2};
          case 'b'
            b = varargin{2};
          otherwise
            error ("makedist: unknown parameter for 'Gamma' distribution.");
        endswitch
        varargin([1:2]) = [];
      endwhile
      pd = GammaDistribution (a, b);

    case {'generalizedextremevalue', 'gev'}
      k = 0;
      sigma = 1;
      mu = 0;
      while (numel (varargin) > 0)
        switch (tolower (varargin{1}))
          case 'k'
            k = varargin{2};
          case 'sigma'
            sigma = varargin{2};
          case 'mu'
            mu = varargin{2};
          otherwise
            error (strcat ("makedist: unknown parameter for", ...
                           " 'GeneralizedExtremeValue' distribution."));
        endswitch
        varargin([1:2]) = [];
      endwhile
      pd = GeneralizedExtremeValueDistribution (k, sigma, mu);

    case {'generalizedpareto', 'gp'}
    k = 1;
    sigma = 1;
    theta = 1;
      while (numel (varargin) > 0)
        switch (tolower (varargin{1}))
          case 'k'
            k = varargin{2};
          case 'sigma'
            sigma = varargin{2};
          case 'theta'
            theta = varargin{2};
          otherwise
            error (strcat ("makedist: unknown parameter for", ...
                           " 'GeneralizedPareto' distribution."));
        endswitch
        varargin([1:2]) = [];
      endwhile
      pd = GeneralizedParetoDistribution (k, sigma, theta);

    case {'halfnormal', 'hn'}
      mu = 0;
      sigma = 1;
      while (numel (varargin) > 0)
        switch (tolower (varargin{1}))
          case 'mu'
            mu = varargin{2};
          case 'sigma'
            sigma = varargin{2};
          otherwise
            error (strcat ("makedist: unknown parameter for", ...
                           " 'HalfNormal' distribution."));
        endswitch
        varargin([1:2]) = [];
      endwhile
      pd = HalfNormalDistribution (mu, sigma);

    case {'inversegaussian', 'invg'}
      mu = 1;
      lambda = 1;
      while (numel (varargin) > 0)
        switch (tolower (varargin{1}))
          case 'mu'
            mu = varargin{2};
          case 'lambda'
            lambda = varargin{2};
          otherwise
            error (strcat ("makedist: unknown parameter for", ...
                           " 'InverseGaussian' distribution."));
        endswitch
        varargin([1:2]) = [];
      endwhile
      pd = InverseGaussianDistribution (mu, lambda);

    case 'logistic'
      mu = 0;
      sigma = 1;
      while (numel (varargin) > 0)
        switch (tolower (varargin{1}))
          case 'mu'
            mu = varargin{2};
          case 'sigma'
            sigma = varargin{2};
          otherwise
            error (strcat ("makedist: unknown parameter for", ...
                           " 'Logistic' distribution."));
        endswitch
        varargin([1:2]) = [];
      endwhile
      pd = LogisticDistribution (mu, sigma);

    case 'loglogistic'
      mu = 0;
      sigma = 1;
      while (numel (varargin) > 0)
        switch (tolower (varargin{1}))
          case 'mu'
            mu = varargin{2};
          case 'sigma'
            sigma = varargin{2};
          otherwise
            error (strcat ("makedist: unknown parameter for", ...
                           " 'Loglogistic' distribution."));
        endswitch
        varargin([1:2]) = [];
      endwhile
      pd = LoglogisticDistribution (mu, sigma);

    case 'lognormal'
      mu = 0;
      sigma = 1;
      while (numel (varargin) > 0)
        switch (tolower (varargin{1}))
          case 'mu'
            mu = varargin{2};
          case 'sigma'
            sigma = varargin{2};
          otherwise
            error (strcat ("makedist: unknown parameter for", ...
                           " 'Lognormal' distribution."));
        endswitch
        varargin([1:2]) = [];
      endwhile
      pd = LognormalDistribution (mu, sigma);

    case 'loguniform'
      lower = 1;
      upper = 4;
      while (numel (varargin) > 0)
        switch (tolower (varargin{1}))
          case 'lower'
            lower = varargin{2};
          case 'upper'
            upper = varargin{2};
          otherwise
            error (strcat ("makedist: unknown parameter for", ...
                           " 'Loguniform' distribution."));
        endswitch
        varargin([1:2]) = [];
      endwhile
      pd = LoguniformDistribution (lower, upper);

    case 'multinomial'
      probs = [0.5, 0.5];
      while (numel (varargin) > 0)
        switch (tolower (varargin{1}))
          case 'probabilities'
            probs = varargin{2};
          otherwise
            error (strcat ("makedist: unknown parameter for", ...
                           " 'Multinomial' distribution."));
        endswitch
        varargin([1:2]) = [];
      endwhile
      pd = MultinomialDistribution (probs);

    case 'nakagami'
      mu = 1;
      omega = 1;
      while (numel (varargin) > 0)
        switch (tolower (varargin{1}))
          case 'mu'
            mu = varargin{2};
          case 'omega'
            omega = varargin{2};
          otherwise
            error (strcat ("makedist: unknown parameter for", ...
                           " 'Nakagami' distribution."));
        endswitch
        varargin([1:2]) = [];
      endwhile
      pd = NakagamiDistribution (mu, omega);

    case {'negativebinomial', 'nbin'}
      R = 1;
      P = 0.5;
      while (numel (varargin) > 0)
        switch (tolower (varargin{1}))
          case 'r'
            R = varargin{2};
          case {'ps', 'p'}
            P = varargin{2};
          otherwise
            error (strcat ("makedist: unknown parameter for", ...
                           " 'NegativeBinomial' distribution."));
        endswitch
        varargin([1:2]) = [];
      endwhile
      pd = NegativeBinomialDistribution (R, P);

    case 'normal'
      mu = 0;
      sigma = 1;
      while (numel (varargin) > 0)
        switch (tolower (varargin{1}))
          case 'mu'
            mu = varargin{2};
          case 'sigma'
            sigma = varargin{2};
          otherwise
            error ("makedist: unknown parameter for 'Normal' distribution.");
        endswitch
        varargin([1:2]) = [];
      endwhile
      pd = NormalDistribution (mu, sigma);

    case 'piecewiselinear'
      x = [0, 1];
      Fx = [0, 1];
      while (numel (varargin) > 0)
        switch (tolower (varargin{1}))
          case 'x'
            x = varargin{2};
          case 'fx'
            Fx = varargin{2};
          otherwise
            error (strcat ("makedist: unknown parameter for", ...
                           " 'PiecewiseLinear' distribution."));
        endswitch
        varargin([1:2]) = [];
      endwhile
      pd = PiecewiseLinearDistribution (x, Fx);

    case 'poisson'
      lambda = 1;
      while (numel (varargin) > 0)
        switch (tolower (varargin{1}))
          case 'lambda'
            lambda = varargin{2};
          otherwise
            error ("makedist: unknown parameter for 'Poisson' distribution.");
        endswitch
        varargin([1:2]) = [];
      endwhile
      pd = PoissonDistribution (lambda);

    case 'rayleigh'
      sigma = 1;
      while (numel (varargin) > 0)
        switch (tolower (varargin{1}))
          case {'sigma', 'b'}
            sigma = varargin{2};
          otherwise
            error (strcat ("makedist: unknown parameter for", ...
                           " 'Rayleigh' distribution."));
        endswitch
        varargin([1:2]) = [];
      endwhile
      pd = RayleighDistribution (sigma);

    case 'rician'
      s = 1;
      sigma = 1;
      while (numel (varargin) > 0)
        switch (tolower (varargin{1}))
          case 's'
            s = varargin{2};
          case 'sigma'
            sigma = varargin{2};
          otherwise
            error ("makedist: unknown parameter for 'Rician' distribution.");
        endswitch
        varargin([1:2]) = [];
      endwhile
      pd = RicianDistribution (s, sigma);

    case 'stable'
      alpha = 2;
      beta = 0;
      gam = 1;
      delta = 0;
      while (numel (varargin) > 0)
        switch (tolower (varargin{1}))
          case {'alpha', 's'}
            alpha = varargin{2};
          case 'beta'
            beta = varargin{2};
          case 'gam'
            gam = varargin{2};
          case 'delta'
            delta = varargin{2};
          otherwise
            error ("makedist: unknown parameter for 'Stable' distribution.");
        endswitch
        varargin([1:2]) = [];
      endwhile
      warning ("makedist: 'Stable' distribution not supported yet.");
      pd = [];

    case {'tlocationscale', 'tls'}
      mu = 0;
      sigma = 1;
      df = 5;
      while (numel (varargin) > 0)
        switch (tolower (varargin{1}))
          case {'mu', 's'}
            mu = varargin{2};
          case 'sigma'
            sigma = varargin{2};
          case {'df', 'nu'}
            df = varargin{2};
          otherwise
            error (strcat ("makedist: unknown parameter for", ...
                           " 'tLocationScale' distribution."));
        endswitch
        varargin([1:2]) = [];
      endwhile
      pd = tLocationScaleDistribution (mu, sigma, df);

    case 'triangular'
      A = 0;
      B = 0.5;
      C = 1;
      while (numel (varargin) > 0)
        switch (tolower (varargin{1}))
          case 'a'
            A = varargin{2};
          case 'b'
            B = varargin{2};
          case 'c'
            C = varargin{2};
          otherwise
            error (strcat ("makedist: unknown parameter for", ...
                           " 'Triangular' distribution."));
        endswitch
        varargin([1:2]) = [];
      endwhile
      pd = TriangularDistribution (A, B, C);

    case 'uniform'
      Lower = 0;
      Upper = 1;
      while (numel (varargin) > 0)
        switch (tolower (varargin{1}))
          case 'lower'
            Lower = varargin{2};
          case 'upper'
            Upper = varargin{2};
          otherwise
            error ("makedist: unknown parameter for 'Uniform' distribution.");
        endswitch
        varargin([1:2]) = [];
      endwhile
      pd = UniformDistribution (Lower, Upper);

    case 'weibull'
      lambda = 1;
      k = 1;
      while (numel (varargin) > 0)
        switch (tolower (varargin{1}))
          case {'lambda', 'a'}
            lambda = varargin{2};
          case {'k', 'b'}
            k = varargin{2};
          otherwise
            error ("makedist: unknown parameter for 'Weibull' distribution.");
        endswitch
        varargin([1:2]) = [];
      endwhile
      pd = WeibullDistribution (lambda, k);

  endswitch

endfunction

## Test output
%!test
%! pd = makedist ('beta');
%! assert_equal (class (pd), "BetaDistribution");
%! assert_equal (pd.a, 1);
%! assert_equal (pd.b, 1);
%!test
%! pd = makedist ('beta', 'a', 5);
%! assert_equal (pd.a, 5);
%! assert_equal (pd.b, 1);
%!test
%! pd = makedist ('beta', 'b', 5);
%! assert_equal (pd.a, 1);
%! assert_equal (pd.b, 5);
%!test
%! pd = makedist ('beta', 'a', 3, 'b', 5);
%! assert_equal (pd.a, 3);
%! assert_equal (pd.b, 5);
%!test
%! pd = makedist ('binomial');
%! assert_equal (class (pd), "BinomialDistribution");
%! assert_equal (pd.N, 1);
%! assert_equal (pd.p, 0.5);
%!test
%! pd = makedist ('binomial', 'N', 5);
%! assert_equal (pd.N, 5);
%! assert_equal (pd.p, 0.5);
%!test
%! pd = makedist ('binomial', 'p', 0.2);
%! assert_equal (pd.N, 1);
%! assert_equal (pd.p, 0.2);
%!test
%! pd = makedist ('binomial', 'N', 3, 'p', 0.3);
%! assert_equal (pd.N, 3);
%! assert_equal (pd.p, 0.3);
%!test
%! pd = makedist ('birnbaumsaunders');
%! assert_equal (class (pd), "BirnbaumSaundersDistribution");
%! assert_equal (pd.beta, 1);
%! assert_equal (pd.gamma, 1);
%!test
%! pd = makedist ('birnbaumsaunders', 'beta', 5);
%! assert_equal (pd.beta, 5);
%! assert_equal (pd.gamma, 1);
%!test
%! pd = makedist ('birnbaumsaunders', 'gamma', 5);
%! assert_equal (pd.beta, 1);
%! assert_equal (pd.gamma, 5);
%!test
%! pd = makedist ('birnbaumsaunders', 'beta', 3, 'gamma', 5);
%! assert_equal (pd.beta, 3);
%! assert_equal (pd.gamma, 5);
%!test
%! pd = makedist ('burr');
%! assert_equal (class (pd), "BurrDistribution");
%! assert_equal (pd.alpha, 1);
%! assert_equal (pd.c, 1);
%! assert_equal (pd.k, 1);
%!test
%! pd = makedist ('burr', 'k', 5);
%! assert_equal (pd.alpha, 1);
%! assert_equal (pd.c, 1);
%! assert_equal (pd.k, 5);
%!test
%! pd = makedist ('burr', 'c', 5);
%! assert_equal (pd.alpha, 1);
%! assert_equal (pd.c, 5);
%! assert_equal (pd.k, 1);
%!test
%! pd = makedist ('burr', 'alpha', 3, 'c', 5);
%! assert_equal (pd.alpha, 3);
%! assert_equal (pd.c, 5);
%! assert_equal (pd.k, 1);
%!test
%! pd = makedist ('burr', 'k', 3, 'c', 5);
%! assert_equal (pd.alpha, 1);
%! assert_equal (pd.c, 5);
%! assert_equal (pd.k, 3);
%!test
%! pd = makedist ('exponential');
%! assert_equal (class (pd), "ExponentialDistribution");
%! assert_equal (pd.mu, 1);
%!test
%! pd = makedist ('exponential', 'mu', 5);
%! assert_equal (pd.mu, 5);
%!test
%! pd = makedist ('extremevalue');
%! assert_equal (class (pd), "ExtremeValueDistribution");
%! assert_equal (pd.mu, 0);
%! assert_equal (pd.sigma, 1);
%!test
%! pd = makedist ('extremevalue', 'mu', 5);
%! assert_equal (class (pd), "ExtremeValueDistribution");
%! assert_equal (pd.mu, 5);
%! assert_equal (pd.sigma, 1);
%!test
%! pd = makedist ('ev', 'sigma', 5);
%! assert_equal (class (pd), "ExtremeValueDistribution");
%! assert_equal (pd.mu, 0);
%! assert_equal (pd.sigma, 5);
%!test
%! pd = makedist ('ev', 'mu', -3, 'sigma', 5);
%! assert_equal (class (pd), "ExtremeValueDistribution");
%! assert_equal (pd.mu, -3);
%! assert_equal (pd.sigma, 5);
%!test
%! pd = makedist ('gamma');
%! assert_equal (class (pd), "GammaDistribution");
%! assert_equal (pd.a, 1);
%! assert_equal (pd.b, 1);
%!test
%! pd = makedist ('gamma', 'a', 5);
%! assert_equal (pd.a, 5);
%! assert_equal (pd.b, 1);
%!test
%! pd = makedist ('gamma', 'b', 5);
%! assert_equal (pd.a, 1);
%! assert_equal (pd.b, 5);
%!test
%! pd = makedist ('gamma', 'a', 3, 'b', 5);
%! assert_equal (pd.a, 3);
%! assert_equal (pd.b, 5);
%!test
%! pd = makedist ('GeneralizedExtremeValue');
%! assert_equal (class (pd), "GeneralizedExtremeValueDistribution");
%! assert_equal (pd.k, 0);
%! assert_equal (pd.sigma, 1);
%! assert_equal (pd.mu, 0);
%!test
%! pd = makedist ('GeneralizedExtremeValue', 'k', 5);
%! assert_equal (pd.k, 5);
%! assert_equal (pd.sigma, 1);
%! assert_equal (pd.mu, 0);
%!test
%! pd = makedist ('GeneralizedExtremeValue', 'sigma', 5);
%! assert_equal (pd.k, 0);
%! assert_equal (pd.sigma, 5);
%! assert_equal (pd.mu, 0);
%!test
%! pd = makedist ('GeneralizedExtremeValue', 'k', 3, 'sigma', 5);
%! assert_equal (pd.k, 3);
%! assert_equal (pd.sigma, 5);
%! assert_equal (pd.mu, 0);
%!test
%! pd = makedist ('GeneralizedExtremeValue', 'mu', 3, 'sigma', 5);
%! assert_equal (pd.k, 0);
%! assert_equal (pd.sigma, 5);
%! assert_equal (pd.mu, 3);
%!test
%! pd = makedist ('GeneralizedPareto');
%! assert_equal (class (pd), "GeneralizedParetoDistribution");
%! assert_equal (pd.k, 1);
%! assert_equal (pd.sigma, 1);
%! assert_equal (pd.theta, 1);
%!test
%! pd = makedist ('GeneralizedPareto', 'k', 5);
%! assert_equal (pd.k, 5);
%! assert_equal (pd.sigma, 1);
%! assert_equal (pd.theta, 1);
%!test
%! pd = makedist ('GeneralizedPareto', 'sigma', 5);
%! assert_equal (pd.k, 1);
%! assert_equal (pd.sigma, 5);
%! assert_equal (pd.theta, 1);
%!test
%! pd = makedist ('GeneralizedPareto', 'k', 3, 'sigma', 5);
%! assert_equal (pd.k, 3);
%! assert_equal (pd.sigma, 5);
%! assert_equal (pd.theta, 1);
%!test
%! pd = makedist ('GeneralizedPareto', 'theta', 3, 'sigma', 5);
%! assert_equal (pd.k, 1);
%! assert_equal (pd.sigma, 5);
%! assert_equal (pd.theta, 3);
%!test
%! pd = makedist ('HalfNormal');
%! assert_equal (class (pd), "HalfNormalDistribution");
%! assert_equal (pd.mu, 0);
%! assert_equal (pd.sigma, 1);
%!test
%! pd = makedist ('HalfNormal', 'mu', 5);
%! assert_equal (pd.mu, 5);
%! assert_equal (pd.sigma, 1);
%!test
%! pd = makedist ('HalfNormal', 'sigma', 5);
%! assert_equal (pd.mu, 0);
%! assert_equal (pd.sigma, 5);
%!test
%! pd = makedist ('HalfNormal', 'mu', 3, 'sigma', 5);
%! assert_equal (pd.mu, 3);
%! assert_equal (pd.sigma, 5);
%!test
%! pd = makedist ('InverseGaussian');
%! assert_equal (class (pd), "InverseGaussianDistribution");
%! assert_equal (pd.mu, 1);
%! assert_equal (pd.lambda, 1);
%!test
%! pd = makedist ('InverseGaussian', 'mu', 5);
%! assert_equal (pd.mu, 5);
%! assert_equal (pd.lambda, 1);
%!test
%! pd = makedist ('InverseGaussian', 'lambda', 5);
%! assert_equal (pd.mu, 1);
%! assert_equal (pd.lambda, 5);
%!test
%! pd = makedist ('InverseGaussian', 'mu', 3, 'lambda', 5);
%! assert_equal (pd.mu, 3);
%! assert_equal (pd.lambda, 5);
%!test
%! pd = makedist ('logistic');
%! assert_equal (class (pd), "LogisticDistribution");
%! assert_equal (pd.mu, 0);
%! assert_equal (pd.sigma, 1);
%!test
%! pd = makedist ('logistic', 'mu', 5);
%! assert_equal (pd.mu, 5);
%! assert_equal (pd.sigma, 1);
%!test
%! pd = makedist ('logistic', 'sigma', 5);
%! assert_equal (pd.mu, 0);
%! assert_equal (pd.sigma, 5);
%!test
%! pd = makedist ('logistic', 'mu', 3, 'sigma', 5);
%! assert_equal (pd.mu, 3);
%! assert_equal (pd.sigma, 5);
%!test
%! pd = makedist ('loglogistic');
%! assert_equal (class (pd), "LoglogisticDistribution");
%! assert_equal (pd.mu, 0);
%! assert_equal (pd.sigma, 1);
%!test
%! pd = makedist ('loglogistic', 'mu', 5);
%! assert_equal (pd.mu, 5);
%! assert_equal (pd.sigma, 1);
%!test
%! pd = makedist ('loglogistic', 'sigma', 5);
%! assert_equal (pd.mu, 0);
%! assert_equal (pd.sigma, 5);
%!test
%! pd = makedist ('loglogistic', 'mu', 3, 'sigma', 5);
%! assert_equal (pd.mu, 3);
%! assert_equal (pd.sigma, 5);
%!test
%! pd = makedist ('Lognormal');
%! assert_equal (class (pd), "LognormalDistribution");
%! assert_equal (pd.mu, 0);
%! assert_equal (pd.sigma, 1);
%!test
%! pd = makedist ('Lognormal', 'mu', 5);
%! assert_equal (pd.mu, 5);
%! assert_equal (pd.sigma, 1);
%!test
%! pd = makedist ('Lognormal', 'sigma', 5);
%! assert_equal (pd.mu, 0);
%! assert_equal (pd.sigma, 5);
%!test
%! pd = makedist ('Lognormal', 'mu', -3, 'sigma', 5);
%! assert_equal (pd.mu, -3);
%! assert_equal (pd.sigma, 5);
%!test
%! pd = makedist ('Loguniform');
%! assert_equal (class (pd), "LoguniformDistribution");
%! assert_equal (pd.Lower, 1);
%! assert_equal (pd.Upper, 4);
%!test
%! pd = makedist ('Loguniform', 'Lower', 2);
%! assert_equal (pd.Lower, 2);
%! assert_equal (pd.Upper, 4);
%!test
%! pd = makedist ('Loguniform', 'Lower', 1, 'Upper', 3);
%! assert_equal (pd.Lower, 1);
%! assert_equal (pd.Upper, 3);
%!test
%! pd = makedist ('Multinomial');
%! assert_equal (class (pd), "MultinomialDistribution");
%! assert_equal (pd.Probabilities, [0.5, 0.5]);
%!test
%! pd = makedist ('Multinomial', 'Probabilities', [0.2, 0.3, 0.1, 0.4]);
%! assert_equal (class (pd), "MultinomialDistribution");
%! assert_equal (pd.Probabilities, [0.2, 0.3, 0.1, 0.4]);
%!test
%! pd = makedist ('Nakagami');
%! assert_equal (class (pd), "NakagamiDistribution");
%! assert_equal (pd.mu, 1);
%! assert_equal (pd.omega, 1);
%!test
%! pd = makedist ('Nakagami', 'mu', 5);
%! assert_equal (class (pd), "NakagamiDistribution");
%! assert_equal (pd.mu, 5);
%! assert_equal (pd.omega, 1);
%!test
%! pd = makedist ('Nakagami', 'omega', 0.3);
%! assert_equal (class (pd), "NakagamiDistribution");
%! assert_equal (pd.mu, 1);
%! assert_equal (pd.omega, 0.3);
%!test
%! pd = makedist ('NegativeBinomial');
%! assert_equal (class (pd), "NegativeBinomialDistribution");
%! assert_equal (pd.R, 1);
%! assert_equal (pd.P, 0.5);
%!test
%! pd = makedist ('NegativeBinomial', 'R', 5);
%! assert_equal (class (pd), "NegativeBinomialDistribution");
%! assert_equal (pd.R, 5);
%! assert_equal (pd.P, 0.5);
%!test
%! pd = makedist ('NegativeBinomial', 'p', 0.3);
%! assert_equal (class (pd), "NegativeBinomialDistribution");
%! assert_equal (pd.R, 1);
%! assert_equal (pd.P, 0.3);
%!test
%! pd = makedist ('Normal');
%! assert_equal (class (pd), "NormalDistribution");
%! assert_equal (pd.mu, 0);
%! assert_equal (pd.sigma, 1);
%!test
%! pd = makedist ('Normal', 'mu', 5);
%! assert_equal (class (pd), "NormalDistribution");
%! assert_equal (pd.mu, 5);
%! assert_equal (pd.sigma, 1);
%!test
%! pd = makedist ('Normal', 'sigma', 5);
%! assert_equal (class (pd), "NormalDistribution");
%! assert_equal (pd.mu, 0);
%! assert_equal (pd.sigma, 5);
%!test
%! pd = makedist ('Normal', 'mu', -3, 'sigma', 5);
%! assert_equal (class (pd), "NormalDistribution");
%! assert_equal (pd.mu, -3);
%! assert_equal (pd.sigma, 5);
%!test
%! pd = makedist ('PiecewiseLinear');
%! assert_equal (class (pd), "PiecewiseLinearDistribution");
%! assert_equal (pd.x, [0; 1]);
%! assert_equal (pd.Fx, [0; 1]);
%!test
%! pd = makedist ('PiecewiseLinear', 'x', [0, 1, 2], 'Fx', [0, 0.5, 1]);
%! assert_equal (pd.x, [0; 1; 2]);
%! assert_equal (pd.Fx, [0; 0.5; 1]);
%!test
%! pd = makedist ('Poisson');
%! assert_equal (class (pd), "PoissonDistribution");
%! assert_equal (pd.lambda, 1);
%!test
%! pd = makedist ('Poisson', 'lambda', 5);
%! assert_equal (pd.lambda, 5);
%!test
%! pd = makedist ('Rayleigh');
%! assert_equal (class (pd), "RayleighDistribution");
%! assert_equal (pd.sigma, 1);
%!test
%! pd = makedist ('Rayleigh', 'sigma', 5);
%! assert_equal (pd.sigma, 5);
%!test
%! pd = makedist ('Rician');
%! assert_equal (class (pd), "RicianDistribution");
%! assert_equal (pd.s, 1);
%! assert_equal (pd.sigma, 1);
%!test
%! pd = makedist ('Rician', 's', 3);
%! assert_equal (pd.s, 3);
%! assert_equal (pd.sigma, 1);
%!test
%! pd = makedist ('Rician', 'sigma', 3);
%! assert_equal (pd.s, 1);
%! assert_equal (pd.sigma, 3);
%!test
%! pd = makedist ('Rician', 's', 2, 'sigma', 3);
%! assert_equal (pd.s, 2);
%! assert_equal (pd.sigma, 3);
%!warning
%! pd = makedist ('stable');
%! assert_equal (class (pd), "double");
%! assert_equal (isempty (pd), true);
%!test
%! pd = makedist ('tlocationscale');
%! assert_equal (class (pd), "tLocationScaleDistribution");
%! assert_equal (pd.mu, 0);
%! assert_equal (pd.sigma, 1);
%! assert_equal (pd.nu, 5);
%!test
%! pd = makedist ('tlocationscale', 'mu', 5);
%! assert_equal (pd.mu, 5);
%! assert_equal (pd.sigma, 1);
%! assert_equal (pd.nu, 5);
%!test
%! pd = makedist ('tlocationscale', 'sigma', 2);
%! assert_equal (pd.mu, 0);
%! assert_equal (pd.sigma, 2);
%! assert_equal (pd.nu, 5);
%!test
%! pd = makedist ('tlocationscale', 'mu', 5, 'sigma', 2);
%! assert_equal (pd.mu, 5);
%! assert_equal (pd.sigma, 2);
%! assert_equal (pd.nu, 5);
%!test
%! pd = makedist ('tlocationscale', 'nu', 1, 'sigma', 2);
%! assert_equal (pd.mu, 0);
%! assert_equal (pd.sigma, 2);
%! assert_equal (pd.nu, 1);
%!test
%! pd = makedist ('tlocationscale', 'mu', -2, 'sigma', 3, 'nu', 1);
%! assert_equal (pd.mu, -2);
%! assert_equal (pd.sigma, 3);
%! assert_equal (pd.nu, 1);
%!test
%! pd = makedist ('Triangular');
%! assert_equal (class (pd), "TriangularDistribution");
%! assert_equal (pd.A, 0);
%! assert_equal (pd.B, 0.5);
%! assert_equal (pd.C, 1);
%!test
%! pd = makedist ('Triangular', 'A', -2);
%! assert_equal (pd.A, -2);
%! assert_equal (pd.B, 0.5);
%! assert_equal (pd.C, 1);
%!test
%! pd = makedist ('Triangular', 'A', 0.5, 'B', 0.9);
%! assert_equal (pd.A, 0.5);
%! assert_equal (pd.B, 0.9);
%! assert_equal (pd.C, 1);
%!test
%! pd = makedist ('Triangular', 'A', 1, 'B', 2, 'C', 5);
%! assert_equal (pd.A, 1);
%! assert_equal (pd.B, 2);
%! assert_equal (pd.C, 5);
%!test
%! pd = makedist ('Uniform');
%! assert_equal (class (pd), "UniformDistribution");
%! assert_equal (pd.Lower, 0);
%! assert_equal (pd.Upper, 1);
%!test
%! pd = makedist ('Uniform', 'Lower', -2);
%! assert_equal (pd.Lower, -2);
%! assert_equal (pd.Upper, 1);
%!test
%! pd = makedist ('Uniform', 'Lower', 1, 'Upper', 3);
%! assert_equal (pd.Lower, 1);
%! assert_equal (pd.Upper, 3);
%!test
%! pd = makedist ('Weibull');
%! assert_equal (class (pd), "WeibullDistribution");
%! assert_equal (pd.lambda, 1);
%! assert_equal (pd.k, 1);
%!test
%! pd = makedist ('Weibull', 'lambda', 3);
%! assert_equal (pd.lambda, 3);
%! assert_equal (pd.k, 1);
%!test
%! pd = makedist ('Weibull', 'lambda', 3, 'k', 2);
%! assert_equal (pd.lambda, 3);
%! assert_equal (pd.k, 2);

## Test input validation
%!error <makedist: DISTNAME must be a character vector.> makedist (1)
%!error <makedist: DISTNAME must be a character vector.> makedist (['as';'sd'])
%!error <makedist: unrecognized distribution name.> makedist ('some')
%!error <makedist: optional arguments must be in NAME-VALUE pairs.> ...
%! makedist ('Beta', 'a')
%!error <makedist: unknown parameter for 'Beta' distribution.> ...
%! makedist ('Beta', 'a', 1, 'Q', 23)
%!error <makedist: unknown parameter for 'Binomial' distribution.> ...
%! makedist ('Binomial', 'N', 1, 'Q', 23)
%!error <makedist: unknown parameter for 'BirnbaumSaunders' distribution.> ...
%! makedist ('BirnbaumSaunders', 'N', 1)
%!error <makedist: unknown parameter for 'Burr' distribution.> ...
%! makedist ('Burr', 'lambda', 1, 'sdfs', 34)
%!error <makedist: unknown parameter for 'ExtremeValue' distribution.> ...
%! makedist ('extremevalue', 'mu', 1, 'sdfs', 34)
%!error <makedist: unknown parameter for 'Exponential' distribution.> ...
%! makedist ('exponential', 'mu', 1, 'sdfs', 34)
%!error <makedist: unknown parameter for 'Gamma' distribution.> ...
%! makedist ('Gamma', 'k', 1, 'sdfs', 34)
%!error <makedist: unknown parameter for 'GeneralizedExtremeValue' distribution.> ...
%! makedist ('GeneralizedExtremeValue', 'k', 1, 'sdfs', 34)
%!error <makedist: unknown parameter for 'GeneralizedPareto' distribution.> ...
%! makedist ('GeneralizedPareto', 'k', 1, 'sdfs', 34)
%!error <makedist: unknown parameter for 'HalfNormal' distribution.> ...
%! makedist ('HalfNormal', 'k', 1, 'sdfs', 34)
%!error <makedist: unknown parameter for 'InverseGaussian' distribution.> ...
%! makedist ('InverseGaussian', 'k', 1, 'sdfs', 34)
%!error <makedist: unknown parameter for 'Logistic' distribution.> ...
%! makedist ('Logistic', 'k', 1, 'sdfs', 34)
%!error <makedist: unknown parameter for 'Loglogistic' distribution.> ...
%! makedist ('Loglogistic', 'k', 1, 'sdfs', 34)
%!error <makedist: unknown parameter for 'Lognormal' distribution.> ...
%! makedist ('Lognormal', 'k', 1, 'sdfs', 34)
%!error <makedist: unknown parameter for 'Loguniform' distribution.> ...
%! makedist ('Loguniform', 'k', 1, 'sdfs', 34)
%!error <makedist: unknown parameter for 'Multinomial' distribution.> ...
%! makedist ('Multinomial', 'k', 1, 'sdfs', 34)
%!error <makedist: unknown parameter for 'Nakagami' distribution.> ...
%! makedist ('Nakagami', 'mu', 1, 'sdfs', 34)
%!error <makedist: unknown parameter for 'NegativeBinomial' distribution.> ...
%! makedist ('NegativeBinomial', 'mu', 1, 'sdfs', 34)
%!error <makedist: unknown parameter for 'Normal' distribution.> ...
%! makedist ('Normal', 'mu', 1, 'sdfs', 34)
%!error <makedist: unknown parameter for 'PiecewiseLinear' distribution.> ...
%! makedist ('PiecewiseLinear', 'mu', 1, 'sdfs', 34)
%!error <makedist: unknown parameter for 'Poisson' distribution.> ...
%! makedist ('Poisson', 'mu', 1, 'sdfs', 34)
%!error <makedist: unknown parameter for 'Rayleigh' distribution.> ...
%! makedist ('Rayleigh', 'mu', 1, 'sdfs', 34)
%!error <makedist: unknown parameter for 'Rician' distribution.> ...
%! makedist ('Rician', 'mu', 1, 'sdfs', 34)
%!error <makedist: unknown parameter for 'Stable' distribution.> ...
%! makedist ('Stable', 'mu', 1, 'sdfs', 34)
%!error <makedist: unknown parameter for 'tLocationScale' distribution.> ...
%! makedist ('tLocationScale', 'mu', 1, 'sdfs', 34)
%!error <makedist: unknown parameter for 'Triangular' distribution.> ...
%! makedist ('Triangular', 'mu', 1, 'sdfs', 34)
%!error <makedist: unknown parameter for 'Uniform' distribution.> ...
%! makedist ('Uniform', 'mu', 1, 'sdfs', 34)
%!error <makedist: unknown parameter for 'Weibull' distribution.> ...
%! makedist ('Weibull', 'mu', 1, 'sdfs', 34)
