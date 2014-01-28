## Copyright (C) 2013 Pantxo Diribarne
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
## @item @tab "beta" @tab "beta" @tab 2
## @item @tab "bino" @tab "binomial" @tab 2
## @item @tab "cauchy" @tab @tab 2
## @item @tab "chi2" @tab "chisquare" @tab 1
## @item @tab "discrete" @tab @tab 2
## @item @tab "exp" @tab "exponential" @tab 1
## @item @tab "f" @tab @tab 2
## @item @tab "gam" @tab "gamma" @tab  2
## @item @tab "geo" @tab "geometric" @tab 1
## @item @tab "gev" @tab "generalized extreme value" @tab  3
## @item @tab "hyge" @tab "hypergeometric" @tab 3
## @item @tab "kolmogorov_smirnov" @tab @tab 1
## @item @tab "laplace" @tab @tab 2
## @item @tab "logistic" @tab  @tab 0
## @item @tab "logn" @tab "lognormal" @tab 2
## @item @tab "norm" @tab "normal" @tab 2
## @item @tab "poiss" @tab "poisson" @tab 2
## @item @tab "rayl" @tab "rayleigh" @tab 1
## @item @tab "t" @tab @tab 1
## @item @tab "unif" @tab "uniform" @tab 2
## @item @tab "wbl" @tab "weibull" @tab 2
## @end multitable
## 
## @seealso{betacdf, binocdf, cauchy_cdf, chi2cdf, discrete_cdf,
## expcdf, fcdf, gamcdf, geocdf, gevcdf, hygecdf,
## kolmogorov_smirnov_cdf, laplace_cdf, logistic_cdf, logncdf,
## normcdf, poisscdf, raylcdf, tcdf, unifcdf, wblcdf}
## @end deftypefn

function [retval] = cdf (varargin)
  ## implemented functions
  persistent allcdf = {{"beta", "beta"}, @betacdf, 2, ...
            {"bino", "binomial"}, @binocdf, 2, ...
            {"cauchy"}, @cauchy_cdf, 2, ...
            {"chi2", "chisquare"}, @chi2cdf, 1, ...
            {"discrete"}, @discrete_cdf, 2, ...
            {"exp", "exponential"}, @expcdf, 1, ...
            {"f"}, @fcdf, 2, ...
            {"gam", "gamma"}, @gamcdf, 2, ...
            {"geo", "geometric"}, @geocdf, 1, ...
            {"gev", "generalized extreme value"}, @gevcdf, 3, ...
            {"hyge", "hypergeometric"}, @hygecdf, 3, ...
            {"kolmogorov_smirnov"}, @kolmogorov_smirnov_cdf, 1, ...
            {"laplace"}, @laplace_cdf, 2, ...
            {"logistic"}, @logistic_cdf, 0, ... # ML has 2 args here
            {"logn", "lognormal"}, @logncdf, 2, ...
            {"norm", "normal"}, @normcdf, 2, ...
            {"poiss", "poisson"}, @poisscdf, 2, ...
            {"rayl", "rayleigh"}, @raylcdf, 1, ...
            {"t"}, @tcdf, 1, ...
            {"unif", "uniform"}, @unifcdf, 2, ...
            {"wbl", "weibull"}, @wblcdf, 2};

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

  idx = cellfun (@(x) any (strcmpi (name, x)), cdfnames);
  
  if (any (idx))
    if (nargs == cdfargs{idx})
      retval = feval (cdfhdl{idx}, x, varargin{:});
    else
      error ("cdf: %s requires %d arguments", name, cdfargs{idx})
    endif
  else
    error ("cdf: %s not implemented", name);
  endif
  
endfunction
