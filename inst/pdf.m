## Copyright (C) 2016 Andreas Stahel
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
## @deftypefn  {statistics} @var{retval} = pdf (@var{name}, @var{X}, @dots{})
##
## Return probability density function of @var{name} function for value @var{x}.
##
## This is a wrapper around various @var{name}pdf and @var{name}_pdf
## functions. See the individual functions help to learn the signification of
## the arguments after @var{x}. Supported functions and corresponding number of
## additional arguments are:
##
## @multitable @columnfractions 0.02 0.3 0.45 0.2
## @headitem @tab function @tab alternative @tab args
## @item @tab "bbs" @tab "Birnbaum-Saunders" @tab 3
## @item @tab "beta" @tab @tab 2
## @item @tab "bino" @tab "binomial" @tab 2
## @item @tab "burr" @tab "Burr" @tab 3
## @item @tab "cauchy" @tab "Cauchy" @tab 0  defaults: loc=0,scale=1
## @item @tab "cauchy" @tab "Cauchy" @tab 2
## @item @tab "chi2" @tab "chi-square" @tab 1
## @item @tab "copula" @tab "Copula family" @tab 2
## @item @tab "discrete" @tab "univariate discrete" @tab 2
## @item @tab "empirical" @tab "univariate empirical" @tab 2
## @item @tab "exp" @tab "exponential" @tab 1
## @item @tab "f" @tab @tab 2
## @item @tab "gam" @tab "gamma" @tab  2
## @item @tab "geo" @tab "geometric" @tab 1
## @item @tab "gev" @tab "generalized extreme value" @tab  3
## @item @tab "gp" @tab "generalized Pareto" @tab 3
## @item @tab "hyge" @tab "hypergeometric" @tab 3
## @item @tab "iwish" @tab "inverse Wishart" @tab 2
## @item @tab "iwish" @tab "inverse Wishart" @tab 3 set log_y=true
## @item @tab "jsu" @tab "Johnson SU" @tab 2
## @item @tab "laplace" @tab "Laplace" @tab 0
## @item @tab "logistic" @tab @tab 0
## @item @tab "logn" @tab "lognormal" @tab 0  defaults: mu=0,sigma=1
## @item @tab "logn" @tab "lognormal" @tab 2
## @item @tab "mn" @tab "multinomial" @tab 1
## @item @tab "mvn" @tab "multivariate normal" @tab 0  defaults: mu=0,sigma=1
## @item @tab "mvn" @tab "multivariate normal" @tab 1  defaults: sigma=1
## @item @tab "mvn" @tab "multivariate normal" @tab 2
## @item @tab "mvt" @tab "multivariate Student" @tab 2
## @item @tab "naka" @tab "Nakagami" @tab 2
## @item @tab "nbin" @tab "negative binomial" @tab 2
## @item @tab "norm" @tab "normal" @tab 2
## @item @tab "poiss" @tab "Poisson" @tab 1
## @item @tab "rayl" @tab "Rayleigh" @tab 1
## @item @tab "stdnormal" @tab "standard normal" @tab 0
## @item @tab "t" @tab @tab 1
## @item @tab "tri" @tab "triangular" @tab 3
## @item @tab "unid" @tab "uniform discrete" @tab 1
## @item @tab "unif" @tab "uniform" @tab 0  defaults: a=0,b=1
## @item @tab "unif" @tab "uniform" @tab 2
## @item @tab "vm" @tab "Von Mises" @tab 2
## @item @tab "wbl" @tab "Weibull" @tab 0  defaults: scale=0,shape=1
## @item @tab "wbl" @tab "Weibull" @tab 1  defaults: shape=1
## @item @tab "wbl" @tab "Weibull" @tab 2
## @item @tab "wish" @tab "Wishart" @tab 2
## @item @tab "wish" @tab "Wishart" @tab 3 set log_y=true
## @end multitable
##
## @seealso{cdf, rnd}
## @end deftypefn

function [retval] = pdf (varargin)
  ## implemented functions
  persistent allpdf = { ...
            {"bbs", "Birnbaum-Saunders"}, @bbspdf, 3, ...
            {"beta"}, @betapdf, 2, ...
            {"bino", "binomial"}, @binopdf, 2, ...
            {"burr", "Burr"}, @burrpdf, 3, ...
            {"cauchy", "Cauchy"}, @cauchy_pdf, 0, ... ## loc = 0, scale = 1
            {"cauchy", "Cauchy"}, @cauchy_pdf, 2, ...
            {"chi2", "chi-square"}, @chi2pdf, 1, ...
            {"copula", "Copula family"}, @copulapdf, 2, ...
            {"discrete", "univariate discrete"}, @discrete_pdf, 2, ...
            {"empirical", "univariate empirical"}, @empirical_pdf, 1, ...
            {"exp", "exponential"}, @exppdf, 1, ...
            {"f"}, @fpdf, 2, ...
            {"gam", "gamma"}, @gampdf, 2, ...
            {"geo", "geometric"}, @geopdf, 1, ...
            {"gev", "generalized extreme value"}, @gevpdf, 3, ...
            {"gp", "generalized Pareto"}, @gppdf, 3, ...
            {"hyge", "hypergeometric"}, @hygepdf, 3, ...
            {"iwish", "inverse Wishart"}, @iwishpdf, 2, ...
            {"iwish", "inverse Wishart"}, @iwishpdf, 3, ... ## include log_y
            {"jsu", "Johnson SU"}, @jsupdf, 2, ...
            {"laplace", "Laplace"}, @laplace_pdf, 0, ...
            {"logistic"}, @logistic_pdf, 0, ...
            {"logn", "lognormal"}, @lognpdf, 0, ... ## mu = 0, sigma = 1
            {"logn", "lognormal"}, @lognpdf, 2, ...
            {"mn", "multinomial"}, @mnpdf, 1, ...
            {"mvn", "multivariate normal"}, @mvnpdf, 0, ...
            {"mvn", "multivariate normal"}, @mvnpdf, 1, ... ## include mu
            {"mvn", "multivariate normal"}, @mvnpdf, 2, ... ## include mu, sigma
            {"mvt", "multivariate Student"}, @mvnpdf, 2, ...
            {"naka", "Nakagami"}, @nakapdf, 2, ...
            {"nbin", "negative binomial"}, @nbinpdf, 2, ...
            {"norm", "normal"}, @normpdf, 2, ...
            {"poiss", "Poisson"}, @poisspdf, 1, ...
            {"rayl", "Rayleigh"}, @raylpdf, 1, ...
            {"stdnormal",  "standard normal"}, @stdnormal_pdf, 0, ...
            {"t"}, @tpdf, 1, ...
            {"tri", "triangular"}, @tripdf, 3, ...
            {"unit", "uniform discrete"}, @unidpdf, 1, ...
            {"unif", "uniform"}, @unifpdf, 0, ...
            {"unif", "uniform"}, @unifpdf, 2, ...
            {"vm", "Von Mises"}, @vmpdf, 2, ...
            {"wbl", "Weibull"}, @wblcdf, 0, ...
            {"wbl", "Weibull"}, @wblpdf, 1, ... ## include scale
            {"wbl", "Weibull"}, @wblpdf, 2, ... ## include scale, shape
            {"wish", "Wishart"}, @iwishpdf, 2, ...
            {"wish", "Wishart"}, @iwishpdf, 3}; ## include log_y

  if (numel (varargin) < 2 || ! ischar (varargin{1}))
    print_usage ();
  endif

  name = varargin{1};
  x = varargin{2};

  varargin(1:2) = [];
  nargs = numel (varargin);

  pdfnames = allpdf(1:3:end);
  pdfhdl = allpdf(2:3:end);
  pdfargs = allpdf(3:3:end);

  idx = cellfun (@(x, y)any(strcmpi (name, x) & nargs == y), pdfnames, pdfargs);
  ## Add special list
  special = {"copula", "Copula family"};

  if (any (idx))
    if (nargs == pdfargs{idx} && ! any (strcmpi (name, special)))
      retval = feval (pdfhdl{idx}, x, varargin{:});
    elseif (nargs == pdfargs{idx} && any (strcmpi (name, special)))
      retval = feval (pdfhdl{idx}, varargin{1}, x, varargin{2:end});
    else
      error ("pdf: %s requires %d arguments", name, pdfargs{idx})
    endif
  else
    error ("pdf: %s not implemented", name);
  endif

endfunction

%!test
%! assert(pdf ('norm', 1, 0, 1), normpdf (1, 0, 1))
%!test
%! x = [0.2:0.2:0.6; 0.2:0.2:0.6];
%! theta = [1; 2];
%! assert (pdf ("copula", x, "Clayton", theta), copulapdf ("Clayton", x, theta))
%!test
%! x = [-1, 0, 1, 2, Inf];
%! assert (pdf ("bbs", x, ones (1,5), ones (1,5), zeros (1,5)), ...
%!         bbspdf (x, ones (1,5), ones (1,5), zeros (1,5)))
%!test
%! x = [1 2];
%! mu = [0.5 1.5];
%! sigma = [1.0 0.5; 0.5 1.0];
%! assert (pdf ("multivariate normal", x, mu, sigma), ...
%!         mvnpdf (x, mu, sigma), 0.001)
