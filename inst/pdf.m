## Copyright (C) 2016 Andreas Stahel
## strongly based on cdf.m by  2013 Pantxo Diribarne
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
## @deftypefn {Function File} {@var{retval} =} pdf (@var{name}, @var{X}, @dots{})
## Return probability density function of @var{name} function for value
## @var{x}.
## This is a wrapper around various @var{name}pdf and @var{name}_pdf
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
## @item @tab "poiss" @tab "poisson" @tab 1
## @item @tab "rayl" @tab "rayleigh" @tab 1
## @item @tab "t" @tab @tab 1
## @item @tab "unif" @tab "uniform" @tab 2
## @item @tab "wbl" @tab "weibull" @tab 2
## @end multitable
## 
## @seealso{betapdf, binopdf, cauchy_pdf, chi2pdf, discrete_pdf,
## exppdf, fpdf, gampdf, geopdf, gevpdf, hygepdf, laplace_pdf,
## logistic_pdf, lognpdf, normpdf, poisspdf, raylpdf, tpdf,
## unifpdf, wblpdf}
## @end deftypefn

function [retval] = pdf (varargin)
  ## implemented functions
  persistent allpdf = {{"beta", "beta"}, @betapdf, 2, ...
            {"bino", "binomial"}, @binopdf, 2, ...
            {"cauchy"}, @cauchy_pdf, 2, ...
            {"chi2", "chisquare"}, @chi2pdf, 1, ...
            {"discrete"}, @discrete_pdf, 2, ...
            {"exp", "exponential"}, @exppdf, 1, ...
            {"f"}, @fpdf, 2, ...
            {"gam", "gamma"}, @gampdf, 2, ...
            {"geo", "geometric"}, @geopdf, 1, ...
            {"gev", "generalized extreme value"}, @gevpdf, 3, ...
            {"hyge", "hypergeometric"}, @hygepdf, 3, ...
            {"laplace"}, @laplace_pdf, 1, ...
            {"logistic"}, @logistic_pdf, 0, ... # ML has 2 args here
            {"logn", "lognormal"}, @lognpdf, 2, ...
            {"norm", "normal"}, @normpdf, 2, ...
            {"poiss", "poisson"}, @poisspdf, 1, ...
            {"rayl", "rayleigh"}, @raylpdf, 1, ...
            {"t"}, @tpdf, 1, ...
            {"unif", "uniform"}, @unifpdf, 2, ...
            {"wbl", "weibull"}, @wblpdf, 2};

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

  idx = cellfun (@(x) any (strcmpi (name, x)), pdfnames);
  
  if (any (idx))
    if (nargs == pdfargs{idx})
      retval = feval (pdfhdl{idx}, x, varargin{:});
    else
      error ("pdf: %s requires %d arguments", name, pdfargs{idx})
    endif
  else
    error ("pdf: %s not implemented", name);
  endif
  
endfunction

%!test
%! assert(pdf ('norm', 1, 0, 1), normpdf (1, 0, 1))