## Copyright (C) 2013 Pantxo Diribarne
## Copyright (C) 2022 Andreas Bertsatos <abertsatos@biol.uoa.gr>
##
## This file is part of the statistics package for GNU Octave.
##
## This program is free software; you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn {Function File} {@var{retval} =} cdf (@var{name}, @var{X}, @dots{})
## Return cumulative density function of @var{name} function for value
## @var{x}.
## This is a wrapper around various @var{name}cdf and @var{name}_cdf
## functions. See the individual functions help to learn the signification of
## the arguments after @var{x}. Supported functions and corresponding number of
## additional arguments are:
##
## @multitable @columnfractions 0.02 0.3 0.45 0.2
## @headitem @tab function @tab alternative @tab args
## @item @tab "bbs" @tab "Birnbaum-Saunders" @tab 3
## @item @tab "beta" @tab @tab 2
## @item @tab "bino" @tab "binomial" @tab 2
## @item @tab "bino" @tab "binomial" @tab 3 include 'upper'
## @item @tab "burr" @tab "Burr" @tab 2
## @item @tab "cauchy" @tab "Cauchy" @tab 2
## @item @tab "chi2" @tab "chi-square" @tab 1
## @item @tab "copula" @tab "Copula family" @tab 2
## @item @tab "copula" @tab "Copula family" @tab 3 include nu
## @item @tab "discrete" @tab "univariate discrete" @tab 2
## @item @tab "empirical" @tab "univariate empirical" @tab 1
## @item @tab "exp" @tab "exponential" @tab 1
## @item @tab "f" @tab @tab 2
## @item @tab "gam" @tab "gamma" @tab  2
## @item @tab "geo" @tab "geometric" @tab 1
## @item @tab "gev" @tab "generalized extreme value" @tab  3
## @item @tab "gp" @tab "generalized Pareto" @tab 3
## @item @tab "hyge" @tab "hypergeometric" @tab 3
## @item @tab "jsu" @tab "Johnson SU" @tab 2
## @item @tab "ks" @tab "Kolmogorov-Smirnov" @tab 1
## @item @tab "laplace" @tab "Laplace" @tab 0
## @item @tab "logistic" @tab @tab 0
## @item @tab "logn" @tab "lognormal" @tab 0  defaults: mu=0,sigma=1
## @item @tab "logn" @tab "lognormal" @tab 2
## @item @tab "mvn" @tab "multivariate normal" @tab 2
## @item @tab "mvn" @tab "multivariate normal" @tab 3 include low limit a
## @item @tab "mvt" @tab "multivariate Student" @tab 2
## @item @tab "mvt" @tab "multivariate Student" @tab 3 include low limit a
## @item @tab "naka" @tab "Nakagami" @tab 2
## @item @tab "nbin" @tab "negative binomial" @tab 2
## @item @tab "norm" @tab "normal" @tab 0  defaults: mu=0,sigma=1
## @item @tab "norm" @tab "normal" @tab 2
## @item @tab "poiss" @tab "Poisson" @tab 1
## @item @tab "rayl" @tab "Rayleigh" @tab 1
## @item @tab "stdnormal" @tab "standard normal" @tab 0
## @item @tab "t" @tab @tab 1
## @item @tab "tri" @tab "triangular" @tab 3
## @item @tab "unid" @tab "uniform discrete" @tab 1
## @item @tab "unif" @tab "uniform" @tab 0  defaults: a=0,b=1
## @item @tab "unif" @tab "uniform" @tab 2
## @item @tab "wbl" @tab "Weibull" @tab 2
## @end multitable
##
## @seealso{pdf, rnd}
## @end deftypefn

function [retval] = cdf (varargin)
  ## implemented functions
  persistent allcdf = { ...
            {"bbs", "Birnbaum-Saunders"}, @bbscdf, 3, ...
            {"beta"}, @betacdf, 2, ...
            {"bino", "binomial"}, @binocdf, 2, ...
            {"bino", "binomial"}, @binocdf, 3, ... ## include 'upper'
            {"burr", "Burr"}, @burrcdf, 2, ...
            {"cauchy", "Cauchy"}, @cauchy_cdf, 2, ...
            {"chi2", "chi-square"}, @chi2cdf, 1, ...
            {"copula", "Copula family"}, @copulacdf, 2, ...
            {"copula", "Copula family"}, @copulacdf, 3, ... ## include nu
            {"discrete", "univariate discrete"}, @discrete_cdf, 2, ...
            {"empirical", "univariate empirical"}, @empirical_cdf, 1, ...
            {"exp", "exponential"}, @expcdf, 1, ...
            {"f"}, @fcdf, 2, ...
            {"gam", "gamma"}, @gamcdf, 2, ...
            {"geo", "geometric"}, @geocdf, 1, ...
            {"gev", "generalized extreme value"}, @gevcdf, 3, ...
            {"gp", "generalized Pareto"}, @gpcdf, 3, ...
            {"hyge", "hypergeometric"}, @hygecdf, 3, ...
            {"jsu", "Johnson SU"}, @jsucdf, 2, ...
            {"ks", "kolmogorov-smirnov"}, @kolmogorov_smirnov_cdf, 1, ...
            {"laplace", "Laplace"}, @laplace_cdf, 0, ...
            {"logistic"}, @logistic_cdf, 0, ...
            {"logn", "lognormal"}, @logncdf, 0, ... ## mu = 0, sigma = 1
            {"logn", "lognormal"}, @logncdf, 2, ...
            {"mvn", "multivariate normal"}, @mvncdf, 2, ...
            {"mvn", "multivariate normal"}, @mvncdf, 3, ... ## include alpha
            {"mvt", "multivariate Student"}, @mvncdf, 2, ...
            {"mvt", "multivariate Student"}, @mvncdf, 3, ... ## include alpha
            {"naka", "Nakagami"}, @nakacdf, 2, ...
            {"nbin", "negative binomial"}, @nbincdf, 2, ...
            {"norm", "normal"}, @normcdf, 0, ... ## mu = 0, sigma = 1
            {"norm", "normal"}, @normcdf, 2, ...
            {"poiss", "Poisson"}, @poisscdf, 1, ...
            {"rayl", "Rayleigh"}, @raylcdf, 1, ...
            {"stdnormal",  "standard normal"}, @stdnormal_cdf, 0, ...
            {"t"}, @tcdf, 1, ...
            {"tri", "triangular"}, @tricdf, 3, ...
            {"unit", "uniform discrete"}, @unidcdf, 1, ...
            {"unif", "uniform"}, @unifcdf, 0, ...
            {"unif", "uniform"}, @unifcdf, 2, ...
            {"wbl", "Weibull"}, @wblcdf, 2};

  if (numel (varargin) < 2 || ! ischar (varargin{1}))
    print_usage ();
  endif

  name = varargin{1};
  x = varargin{2};

  varargin(1:2) = [];
  nargs = numel (varargin);

  cdfnames = allcdf(1:3:end);
  cdfhdl = allcdf(2:3:end);
  cdfargs = allcdf(3:3:end);

  idx = cellfun (@(x, y)any(strcmpi (name, x) & nargs == y), cdfnames, cdfargs);
  ## Add special list
  special = {"copula", "Copula family", "mvn", "multivariate normal", ...
             "mvt", "multivariate Student"};

  if (any (idx))
    if (nargs == cdfargs{idx} && ! any (strcmpi (name, special)))
      retval = feval (cdfhdl{idx}, x, varargin{:});
    elseif (nargs == cdfargs{idx} && any (strcmpi (name, special)))
      if ((any (strcmpi (name, special(3:6))) && nargs == 3) ||
          any (strcmpi (name, special(1:2))))
        retval = feval (cdfhdl{idx}, varargin{1}, x, varargin{2:end});
      else
        retval = feval (cdfhdl{idx}, x, varargin{:});
      endif
    else
      error ("cdf: %s requires %d arguments", name, cdfargs{idx})
    endif
  else
    error ("cdf: %s not implemented", name);
  endif

endfunction

%!test
%! assert (cdf ("norm", 1, 0, 1), normcdf (1, 0, 1))
%!test
%! x = [0.2:0.2:0.6; 0.2:0.2:0.6];
%! theta = [1; 2];
%! assert (cdf ("copula", x, "Clayton", theta), copulacdf ("Clayton", x, theta))
%!test
%! x = [-1, 0, 1, 2, Inf];
%! assert (cdf ("bbs", x, ones (1,5), ones (1,5), zeros (1,5)), ...
%!         bbscdf (x, ones (1,5), ones (1,5), zeros (1,5)))
%!test
%! x = [1 2];
%! mu = [0.5 1.5];
%! sigma = [1.0 0.5; 0.5 1.0];
%! assert (cdf ("multivariate normal", x, mu, sigma), ...
%!         mvncdf (x, mu, sigma), 0.01)
%! a = [-inf 0];
%! assert (cdf ("mvn", x, a, mu, sigma), ...
%!         mvncdf (a, x, mu, sigma), 0.01)
