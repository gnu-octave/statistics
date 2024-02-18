## Copyright (C) 2024 Andreas Bertsatos <abertsatos@biol.uoa.gr>
##
## This file is part of the statistics package for GNU Octave.
##
## This program is free software: you can redistribute it and/or
## modify it under the terms of the GNU General Public License as
## published by the Free Software Foundation, either version 3 of the
## License, or (at your option) any later version.
##
## This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; see the file COPYING.  If not, see
## <http://www.gnu.org/licenses/>.

function pd = makedist (varargin)

  ## Add list of supported probability distribution objects
  PDO = {'Beta'; 'Binomial'; 'BirnbaumSaunders'; 'Burr'; 'ExtremeValue'; ...
         'Exponential'; 'Gamma'; 'GeneralizedExtremeValue'; ...
         'GeneralizedPareto'; 'HalfNormal'; 'InverseGaussian'; ...
         'Logistic'; 'Loglogistic'; 'Lognormal'; 'Loguniform'; ...
         'Multinomial'; 'Nakagami'; 'NegativeBinomial'; 'Normal'; ...
         'PiecewiseLinear'; 'Poisson'; 'Rayleigh'; 'Rician'; ...
         'Stable'; 'tLocationScale'; 'Triangular'; 'Uniform'; 'Weibull'};

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
  elseif (! any (strcmpi (distname, PDO)))
    error ("makedist: unrecognized distribution name.");
  endif

  ## Check for additional arguments being in pairs
  if (mod (numel (varargin), 2) != 0)
    error ("makedist: optional arguments must be in NAME-VALUE pairs.");
  endif

  ## Switch to selected distribution
  switch (tolower (distname))

    case "beta"
      ## Add default parameters
      a = 1;
      b = 1;
      while (numel (varargin) > 0)
        switch (tolower (varargin{1}))
          case "a"
            a = varargin{2};
          case "b"
            b = varargin{2};
          otherwise
            error ("makedist: unknown parameter for 'Beta' distribution.");
        endswitch
        varargin([1:2]) = [];
      endwhile
      pd = [];

    case "binomial"
      N = 1;
      ps = 0.5;
      while (numel (varargin) > 0)
        switch (tolower (varargin{1}))
          case "N"
            N = varargin{2};
          case {"ps", "p"}
            ps = varargin{2};
          otherwise
            error ("makedist: unknown parameter for 'Binomial' distribution.");
        endswitch
        varargin([1:2]) = [];
      endwhile
      pd = [];

    case "birnbaumsaunders"
      beta = 1;
      gamma = 1;
      while (numel (varargin) > 0)
        switch (tolower (varargin{1}))
          case "beta"
            beta = varargin{2};
          case "gamma"
            gamma = varargin{2};
          otherwise
            error (strcat (["makedist: unknown parameter for"], ...
                           [" 'BirnbaumSaunders' distribution."]));
        endswitch
        varargin([1:2]) = [];
      endwhile
      pd = [];

    case "burr"
      lambda = 1;
      c = 1;
      k = 1;
      while (numel (varargin) > 0)
        switch (tolower (varargin{1}))
          case {"lambda", "alpha"}
            lambda = varargin{2};
          case "c"
            c = varargin{2};
          case "k"
            k = varargin{2};
          otherwise
            error ("makedist: unknown parameter for 'Burr' distribution.");
        endswitch
        varargin([1:2]) = [];
      endwhile
      pd = [];

    case "gamma"
      k = 1;
      theta = 1;
      while (numel (varargin) > 0)
        switch (tolower (varargin{1}))
          case {"k", "a"}
            k = varargin{2};
          case {"theta", "b"}
            theta = varargin{2};
          otherwise
            error ("makedist: unknown parameter for 'Gamma' distribution.");
        endswitch
        varargin([1:2]) = [];
      endwhile
      pd = [];

    case "extremevalue"
      mu = 0;
      while (numel (varargin) > 0)
        switch (tolower (varargin{1}))
          case "mu"
            mu = varargin{2};
          otherwise
            error (strcat (["makedist: unknown parameter for"], ...
                           [" 'ExtremeValue' distribution."]));
        endswitch
        varargin([1:2]) = [];
      endwhile
      pd = [];

    case "exponential"
      mu = 0;
      while (numel (varargin) > 0)
        switch (tolower (varargin{1}))
          case "mu"
            mu = varargin{2};
          otherwise
            error (strcat (["makedist: unknown parameter for"], ...
                           [" 'Exponential' distribution."]));
        endswitch
        varargin([1:2]) = [];
      endwhile
      pd = [];

    case "generalizedextremevalue"
      k = 0;
      sigma = 1;
      mu = 0;
      while (numel (varargin) > 0)
        switch (tolower (varargin{1}))
          case "k"
            k = varargin{2};
          case "sigma"
            sigma = varargin{2};
          case "mu"
            mu = varargin{2};
          otherwise
            error (strcat (["makedist: unknown parameter for"], ...
                           [" 'GeneralizedExtremeValue' distribution."]));
        endswitch
        varargin([1:2]) = [];
      endwhile
      pd = [];

    case "generalizedpareto"
    k = 1;
    sigma = 1;
    theta = 1;
      while (numel (varargin) > 0)
        switch (tolower (varargin{1}))
          case "k"
            k = varargin{2};
          case "sigma"
            sigma = varargin{2};
          case "theta"
            theta = varargin{2};
          otherwise
            error (strcat (["makedist: unknown parameter for"], ...
                           [" 'GeneralizedPareto' distribution."]));
        endswitch
        varargin([1:2]) = [];
      endwhile
      pd = [];

    case "halfnormal"
      mu = 0;
      sigma = 1;
      while (numel (varargin) > 0)
        switch (tolower (varargin{1}))
          case "mu"
            mu = varargin{2};
          case "sigma"
            sigma = varargin{2};
          otherwise
            error (strcat (["makedist: unknown parameter for"], ...
                           [" 'HalfNormal' distribution."]));
        endswitch
        varargin([1:2]) = [];
      endwhile
      pd = [];

    case "inversegaussian"
      mu = 1;
      lambda = 1;
      while (numel (varargin) > 0)
        switch (tolower (varargin{1}))
          case "mu"
            mu = varargin{2};
          case "lambda"
            lambda = varargin{2};
          otherwise
            error (strcat (["makedist: unknown parameter for"], ...
                           [" 'InverseGaussian' distribution."]));
        endswitch
        varargin([1:2]) = [];
      endwhile
      pd = [];

    case "logistic"
      mu = 0;
      sigma = 1;
      while (numel (varargin) > 0)
        switch (tolower (varargin{1}))
          case "mu"
            mu = varargin{2};
          case "sigma"
            sigma = varargin{2};
          otherwise
            error (strcat (["makedist: unknown parameter for"], ...
                           [" 'Logistic' distribution."]));
        endswitch
        varargin([1:2]) = [];
      endwhile
      pd = [];

    case "loglogistic"
      mu = 0;
      sigma = 1;
      while (numel (varargin) > 0)
        switch (tolower (varargin{1}))
          case "mu"
            mu = varargin{2};
          case "sigma"
            sigma = varargin{2};
          otherwise
            error (strcat (["makedist: unknown parameter for"], ...
                           [" 'Loglogistic' distribution."]));
        endswitch
        varargin([1:2]) = [];
      endwhile
      pd = [];

    case "lognormal"
      mu = 0;
      sigma = 1;
      while (numel (varargin) > 0)
        switch (tolower (varargin{1}))
          case "mu"
            mu = varargin{2};
          case "sigma"
            sigma = varargin{2};
          otherwise
            error (strcat (["makedist: unknown parameter for"], ...
                           [" 'Lognormal' distribution."]));
        endswitch
        varargin([1:2]) = [];
      endwhile
      pd = [];

    case "loguniform"
      lower = 1;
      upper = 4;
      while (numel (varargin) > 0)
        switch (tolower (varargin{1}))
          case "lower"
            lower = varargin{2};
          case "upper"
            upper = varargin{2};
          otherwise
            error (strcat (["makedist: unknown parameter for"], ...
                           [" 'Loguniform' distribution."]));
        endswitch
        varargin([1:2]) = [];
      endwhile
      pd = [];

    case "multinomial"
      probs = [0.5, 0.5];
      while (numel (varargin) > 0)
        switch (tolower (varargin{1}))
          case "probabilities"
            probs = varargin{2};
          otherwise
            error (strcat (["makedist: unknown parameter for"], ...
                           [" 'Multinomial' distribution."]));
        endswitch
        varargin([1:2]) = [];
      endwhile
      pd = [];

    case "nakagami"
      mu = 0;
      omega = 1;
      while (numel (varargin) > 0)
        switch (tolower (varargin{1}))
          case "mu"
            mu = varargin{2};
          case "omega"
            omega = varargin{2};
          otherwise
            error (strcat (["makedist: unknown parameter for"], ...
                           [" 'Nakagami' distribution."]));
        endswitch
        varargin([1:2]) = [];
      endwhile
      pd = [];

    case "negativebinomial"
      r = 1;
      ps = 0;
      while (numel (varargin) > 0)
        switch (tolower (varargin{1}))
          case "r"
            r = varargin{2};
          case {"ps", "p"}
            ps = varargin{2};
          otherwise
            error (strcat (["makedist: unknown parameter for"], ...
                           [" 'NegativeBinomial' distribution."]));
        endswitch
        varargin([1:2]) = [];
      endwhile
      pd = [];

    case "normal"
      mu = 0;
      sigma = 1;
      while (numel (varargin) > 0)
        switch (tolower (varargin{1}))
          case "mu"
            mu = varargin{2};
          case "sigma"
            sigma = varargin{2};
          otherwise
            error ("makedist: unknown parameter for 'Normal' distribution.");
        endswitch
        varargin([1:2]) = [];
      endwhile
      pd = [];

    case "piecewiselinear"
      x = [0, 1];
      Fx = [0, 1];
      while (numel (varargin) > 0)
        switch (tolower (varargin{1}))
          case "x"
            x = varargin{2};
          case "fx"
            Fx = varargin{2};
          otherwise
            error (strcat (["makedist: unknown parameter for"], ...
                           [" 'PiecewiseLinear' distribution."]));
        endswitch
        varargin([1:2]) = [];
      endwhile
      pd = PiecewiseLinearDistribution (x, Fx);

    case "poisson"
      lambda = 1;
      while (numel (varargin) > 0)
        switch (tolower (varargin{1}))
          case "lambda"
            lambda = varargin{2};
          otherwise
            error ("makedist: unknown parameter for 'Poisson' distribution.");
        endswitch
        varargin([1:2]) = [];
      endwhile
      pd = PoissonDistribution (lambda);

    case "rayleigh"
      sigma = 1;
      while (numel (varargin) > 0)
        switch (tolower (varargin{1}))
          case {"sigma", "b"}
            sigma = varargin{2};
          otherwise
            error (strcat (["makedist: unknown parameter for"], ...
                           [" 'Rayleigh' distribution."]));
        endswitch
        varargin([1:2]) = [];
      endwhile
      pd = RayleighDistribution (sigma);

    case "rician"
      nu = 1;
      sigma = 1;
      while (numel (varargin) > 0)
        switch (tolower (varargin{1}))
          case {"nu", "s"}
            nu = varargin{2};
          case "sigma"
            sigma = varargin{2};
          otherwise
            error ("makedist: unknown parameter for 'Rician' distribution.");
        endswitch
        varargin([1:2]) = [];
      endwhile
      pd = RicianDistribution (nu, sigma);

    case "stable"
      alpha = 2;
      beta = 0;
      gam = 1;
      delta = 0;
      while (numel (varargin) > 0)
        switch (tolower (varargin{1}))
          case {"alpha", "s"}
            alpha = varargin{2};
          case "beta"
            beta = varargin{2};
          case "gam"
            gam = varargin{2};
          case "delta"
            delta = varargin{2};
          otherwise
            error ("makedist: unknown parameter for 'Stable' distribution.");
        endswitch
        varargin([1:2]) = [];
      endwhile
      warning ("makedist: 'Stable' distribution not supported yet.");
      pd = [];

    case "tlocationscale"
      mu = 0;
      sigma = 1;
      df = 5;
      while (numel (varargin) > 0)
        switch (tolower (varargin{1}))
          case {"mu", "s"}
            mu = varargin{2};
          case "sigma"
            sigma = varargin{2};
          case {"df", "nu"}
            df = varargin{2};
          otherwise
            error (strcat (["makedist: unknown parameter for"], ...
                           [" 'tLocationScale' distribution."]));
        endswitch
        varargin([1:2]) = [];
      endwhile
      pd = tLocationScaleDistribution (mu, sigma, df);

    case "triangular"
      A = 0;
      B = 0.5;
      C = 1;
      while (numel (varargin) > 0)
        switch (tolower (varargin{1}))
          case "a"
            A = varargin{2};
          case "b"
            B = varargin{2};
          case "c"
            C = varargin{2};
          otherwise
            error (strcat (["makedist: unknown parameter for"], ...
                           [" 'Triangular' distribution."]));
        endswitch
        varargin([1:2]) = [];
      endwhile
      pd = TriangularDistribution (A, B, C);

    case "uniform"
      Lower = 0;
      Upper = 1;
      while (numel (varargin) > 0)
        switch (tolower (varargin{1}))
          case "lower"
            Lower = varargin{2};
          case "upper"
            Upper = varargin{2};
          otherwise
            error ("makedist: unknown parameter for 'Uniform' distribution.");
        endswitch
        varargin([1:2]) = [];
      endwhile
      pd = UniformDistribution (Lower, Upper);

    case "weibull"
      lambda = 1;
      k = 1;
      while (numel (varargin) > 0)
        switch (tolower (varargin{1}))
          case {"lambda", "a"}
            lambda = varargin{2};
          case {"k", "b"}
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
%! pd = makedist ("PiecewiseLinear");
%! assert (class (pd), "PiecewiseLinearDistribution");
%! assert (pd.x, [0; 1]);
%! assert (pd.Fx, [0; 1]);
%!test
%! pd = makedist ("PiecewiseLinear", "x", [0, 1, 2], "Fx", [0, 0.5, 1]);
%! assert (pd.x, [0; 1; 2]);
%! assert (pd.Fx, [0; 0.5; 1]);
%!test
%! pd = makedist ("Poisson");
%! assert (class (pd), "PoissonDistribution");
%! assert (pd.lambda, 1);
%!test
%! pd = makedist ("Poisson", "lambda", 5);
%! assert (pd.lambda, 5);
%!test
%!test
%! pd = makedist ("Rayleigh");
%! assert (class (pd), "RayleighDistribution");
%! assert (pd.sigma, 1);
%!test
%! pd = makedist ("Rayleigh", "sigma", 5);
%! assert (pd.sigma, 5);
%!test
%! pd = makedist ("Rician");
%! assert (class (pd), "RicianDistribution");
%! assert (pd.nu, 1);
%! assert (pd.sigma, 1);
%!test
%! pd = makedist ("Rician", "nu", 3);
%! assert (pd.nu, 3);
%! assert (pd.sigma, 1);
%!test
%! pd = makedist ("Rician", "sigma", 3);
%! assert (pd.nu, 1);
%! assert (pd.sigma, 3);
%!test
%! pd = makedist ("Rician", "nu", 2, "sigma", 3);
%! assert (pd.nu, 2);
%! assert (pd.sigma, 3);
%!warning
%! pd = makedist ("stable");
%! assert (class (pd), "double");
%! assert (isempty (pd), true);
%!test
%! pd = makedist ("tlocationscale");
%! assert (class (pd), "tLocationScaleDistribution");
%! assert (pd.mu, 0);
%! assert (pd.sigma, 1);
%! assert (pd.nu, 5);
%!test
%! pd = makedist ("tlocationscale", "mu", 5);
%! assert (pd.mu, 5);
%! assert (pd.sigma, 1);
%! assert (pd.nu, 5);
%!test
%! pd = makedist ("tlocationscale", "sigma", 2);
%! assert (pd.mu, 0);
%! assert (pd.sigma, 2);
%! assert (pd.nu, 5);
%!test
%! pd = makedist ("tlocationscale", "mu", 5, "sigma", 2);
%! assert (pd.mu, 5);
%! assert (pd.sigma, 2);
%! assert (pd.nu, 5);
%!test
%! pd = makedist ("tlocationscale", "nu", 1, "sigma", 2);
%! assert (pd.mu, 0);
%! assert (pd.sigma, 2);
%! assert (pd.nu, 1);
%!test
%! pd = makedist ("tlocationscale", "mu", -2, "sigma", 3, "nu", 1);
%! assert (pd.mu, -2);
%! assert (pd.sigma, 3);
%! assert (pd.nu, 1);
%!test
%! pd = makedist ("Triangular");
%! assert (class (pd), "TriangularDistribution");
%! assert (pd.A, 0);
%! assert (pd.B, 0.5);
%! assert (pd.C, 1);
%!test
%! pd = makedist ("Triangular", "A", -2);
%! assert (pd.A, -2);
%! assert (pd.B, 0.5);
%! assert (pd.C, 1);
%!test
%! pd = makedist ("Triangular", "A", 0.5, "B", 0.9);
%! assert (pd.A, 0.5);
%! assert (pd.B, 0.9);
%! assert (pd.C, 1);
%!test
%! pd = makedist ("Triangular", "A", 1, "B", 2, "C", 5);
%! assert (pd.A, 1);
%! assert (pd.B, 2);
%! assert (pd.C, 5);
%!test
%! pd = makedist ("Uniform");
%! assert (class (pd), "UniformDistribution");
%! assert (pd.Lower, 0);
%! assert (pd.Upper, 1);
%!test
%! pd = makedist ("Uniform", "Lower", -2);
%! assert (pd.Lower, -2);
%! assert (pd.Upper, 1);
%!test
%! pd = makedist ("Uniform", "Lower", 1, "Upper", 3);
%! assert (pd.Lower, 1);
%! assert (pd.Upper, 3);
%!test
%! pd = makedist ("Weibull");
%! assert (class (pd), "WeibullDistribution");
%! assert (pd.lambda, 1);
%! assert (pd.k, 1);
%!test
%! pd = makedist ("Weibull", "lambda", 3);
%! assert (pd.lambda, 3);
%! assert (pd.k, 1);
%!test
%! pd = makedist ("Weibull", "lambda", 3, "k", 2);
%! assert (pd.lambda, 3);
%! assert (pd.k, 2);

## Test input validation
%!error <makedist: DISTNAME must be a character vector.> makedist (1)
%!error <makedist: DISTNAME must be a character vector.> makedist (["as";"sd"])
%!error <makedist: unrecognized distribution name.> makedist ("some")
%!error <makedist: optional arguments must be in NAME-VALUE pairs.> ...
%! makedist ("Beta", "a")
%!error <makedist: unknown parameter for 'Beta' distribution.> ...
%! makedist ("Beta", "a", 1, "Q", 23)
%!error <makedist: unknown parameter for 'Binomial' distribution.> ...
%! makedist ("Binomial", "N", 1, "Q", 23)
%!error <makedist: unknown parameter for 'BirnbaumSaunders' distribution.> ...
%! makedist ("BirnbaumSaunders", "N", 1)
%!error <makedist: unknown parameter for 'Burr' distribution.> ...
%! makedist ("Burr", "lambda", 1, "sdfs", 34)
%!error <makedist: unknown parameter for 'ExtremeValue' distribution.> ...
%! makedist ("extremevalue", "mu", 1, "sdfs", 34)
%!error <makedist: unknown parameter for 'Exponential' distribution.> ...
%! makedist ("exponential", "mu", 1, "sdfs", 34)
%!error <makedist: unknown parameter for 'Gamma' distribution.> ...
%! makedist ("Gamma", "k", 1, "sdfs", 34)
%!error <makedist: unknown parameter for 'GeneralizedExtremeValue' distribution.> ...
%! makedist ("GeneralizedExtremeValue", "k", 1, "sdfs", 34)
%!error <makedist: unknown parameter for 'GeneralizedPareto' distribution.> ...
%! makedist ("GeneralizedPareto", "k", 1, "sdfs", 34)
%!error <makedist: unknown parameter for 'HalfNormal' distribution.> ...
%! makedist ("HalfNormal", "k", 1, "sdfs", 34)
%!error <makedist: unknown parameter for 'InverseGaussian' distribution.> ...
%! makedist ("InverseGaussian", "k", 1, "sdfs", 34)
%!error <makedist: unknown parameter for 'Logistic' distribution.> ...
%! makedist ("Logistic", "k", 1, "sdfs", 34)
%!error <makedist: unknown parameter for 'Loglogistic' distribution.> ...
%! makedist ("Loglogistic", "k", 1, "sdfs", 34)
%!error <makedist: unknown parameter for 'Lognormal' distribution.> ...
%! makedist ("Lognormal", "k", 1, "sdfs", 34)
%!error <makedist: unknown parameter for 'Loguniform' distribution.> ...
%! makedist ("Loguniform", "k", 1, "sdfs", 34)
%!error <makedist: unknown parameter for 'Multinomial' distribution.> ...
%! makedist ("Multinomial", "k", 1, "sdfs", 34)
%!error <makedist: unknown parameter for 'Nakagami' distribution.> ...
%! makedist ("Nakagami", "mu", 1, "sdfs", 34)
%!error <makedist: unknown parameter for 'NegativeBinomial' distribution.> ...
%! makedist ("NegativeBinomial", "mu", 1, "sdfs", 34)
%!error <makedist: unknown parameter for 'Normal' distribution.> ...
%! makedist ("Normal", "mu", 1, "sdfs", 34)
%!error <makedist: unknown parameter for 'PiecewiseLinear' distribution.> ...
%! makedist ("PiecewiseLinear", "mu", 1, "sdfs", 34)
%!error <makedist: unknown parameter for 'Poisson' distribution.> ...
%! makedist ("Poisson", "mu", 1, "sdfs", 34)
%!error <makedist: unknown parameter for 'Rayleigh' distribution.> ...
%! makedist ("Rayleigh", "mu", 1, "sdfs", 34)
%!error <makedist: unknown parameter for 'Rician' distribution.> ...
%! makedist ("Rician", "mu", 1, "sdfs", 34)
%!error <makedist: unknown parameter for 'Stable' distribution.> ...
%! makedist ("Stable", "mu", 1, "sdfs", 34)
%!error <makedist: unknown parameter for 'tLocationScale' distribution.> ...
%! makedist ("tLocationScale", "mu", 1, "sdfs", 34)
%!error <makedist: unknown parameter for 'Triangular' distribution.> ...
%! makedist ("Triangular", "mu", 1, "sdfs", 34)
%!error <makedist: unknown parameter for 'Uniform' distribution.> ...
%! makedist ("Uniform", "mu", 1, "sdfs", 34)
%!error <makedist: unknown parameter for 'Weibull' distribution.> ...
%! makedist ("Weibull", "mu", 1, "sdfs", 34)
