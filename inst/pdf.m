## Copyright (C) 2016 Andreas Stahel
## Copyright (C) 2022-2023 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {statistics} {@var{y} =} pdf (@var{name}, @var{x}, @var{A})
## @deftypefnx {statistics} {@var{y} =} pdf (@var{name}, @var{x}, @var{A}, @var{B})
## @deftypefnx {statistics} {@var{y} =} pdf (@var{name}, @var{x}, @var{A}, @var{B}, @var{C})
##
## Return the PDF of a univariate distribution evaluated at @var{x}.
##
## @code{pdf} is a wrapper for the univariate cumulative distribution functions
## available in the statistics package.  See the corresponding functions' help
## to learn the signification of the parameters after @var{x}.
##
## @code{@var{y} = pdf (@var{name}, @var{x}, @var{A})} returns the CDF for the
## one-parameter distribution family specified by @var{name} and the
## distribution parameter @var{A}, evaluated at the values in @var{x}.
##
## @code{@var{y} = pdf (@var{name}, @var{x}, @var{A}, @var{B})} returns the CDF
## for the two-parameter distribution family specified by @var{name} and the
## distribution parameters @var{A} and @var{B}, evaluated at the values in
## @var{x}.
##
## @code{@var{y} = pdf (@var{name}, @var{x}, @var{A}, @var{B}, @var{C})} returns
## the CDF for the three-parameter distribution family specified by @var{name}
## and the distribution parameters @var{A}, @var{B}, and @var{C}, evaluated at
## the values in @var{x}.
##
## @var{name} must be a char string of the name or the abbreviation of the
## desired cumulative distribution function as listed in the followng table.
## The last column shows the number of required parameters that should be parsed
## after @var{x} to the desired PDF.
##
## @multitable @columnfractions 0.4 0.05 0.2 0.05 0.3
## @headitem Distribution Name @tab @tab Abbreviation @tab @tab Input Parameters
## @item @qcode{"Beta"} @tab @tab @qcode{"beta"} @tab @tab 2
## @item @qcode{"Binomial"} @tab @tab @qcode{"bino"} @tab @tab 2
## @item @qcode{"Birnbaum-Saunders"} @tab @tab @qcode{"bisa"} @tab @tab 2
## @item @qcode{"Burr"} @tab  @tab @qcode{"burr"} @tab  @tab 3
## @item @qcode{"Cauchy"} @tab @tab @qcode{"cauchy"} @tab @tab 2
## @item @qcode{"Chi-square"} @tab @tab @qcode{"chi2"} @tab @tab 1
## @item @qcode{"Extreme Value"} @tab @tab @qcode{"ev"} @tab @tab 2
## @item @qcode{"Exponential"} @tab @tab @qcode{"exp"} @tab @tab 1
## @item @qcode{"F-Distribution"} @tab @tab @qcode{"f"} @tab @tab 2
## @item @qcode{"Gamma"} @tab @tab @qcode{"gam"} @tab @tab 2
## @item @qcode{"Geometric"} @tab @tab @qcode{"geo"} @tab @tab 1
## @item @qcode{"Generalized Extreme Value"} @tab @tab @qcode{"gev"} @tab @tab 3
## @item @qcode{"Generalized Pareto"} @tab @tab @qcode{"gp"} @tab @tab 3
## @item @qcode{"Gumbel"} @tab @tab @qcode{"gumbel"} @tab @tab 2
## @item @qcode{"Hypergeometric"} @tab @tab @qcode{"hyge"} @tab @tab 3
## @item @qcode{"Laplace"} @tab @tab @qcode{"laplace"} @tab @tab 2
## @item @qcode{"Logistic"} @tab @tab @qcode{"logi"} @tab @tab 2
## @item @qcode{"Log-Logistic"} @tab @tab @qcode{"logl"} @tab @tab 2
## @item @qcode{"Lognormal"} @tab @tab @qcode{"logn"} @tab @tab 2
## @item @qcode{"Nakagami"} @tab @tab @qcode{"naka"} @tab @tab 2
## @item @qcode{"Negative Binomial"} @tab @tab @qcode{"nbin"} @tab @tab 2
## @item @qcode{"Noncentral F-Distribution"} @tab @tab @qcode{"ncf"} @tab @tab 3
## @item @qcode{"Noncentral Student T"} @tab @tab @qcode{"nct"} @tab @tab 2
## @item @qcode{"Noncentral Chi-Square"} @tab @tab @qcode{"ncx2"} @tab @tab 2
## @item @qcode{"Normal"} @tab @tab @qcode{"norm"} @tab @tab 2
## @item @qcode{"Poisson"} @tab @tab @qcode{"poiss"} @tab @tab 1
## @item @qcode{"Rayleigh"} @tab @tab @qcode{"rayl"} @tab @tab 1
## @item @qcode{"Student T"} @tab @tab @qcode{"t"} @tab @tab 1
## @item @qcode{"Triangular"} @tab @tab @qcode{"tri"} @tab @tab 3
## @item @qcode{"Discrete Uniform"} @tab @tab @qcode{"unid"} @tab @tab 1
## @item @qcode{"Uniform"} @tab @tab @qcode{"unif"} @tab @tab 2
## @item @qcode{"Weibull"} @tab @tab @qcode{"wbl"} @tab @tab 2
## @end multitable
##
## @seealso{cdf, icdf, random, betapdf, binopdf, bisapdf, burrpdf, cauchypdf,
## chi2pdf, evpdf, exppdf, fpdf, gampdf, geopdf, gevpdf, gppdf, hygepdf,
## laplacepdf, logipdf, loglpdf, lognpdf, nakapdf, nbinpdf, ncfpdf, nctpdf,
## ncx2pdf, normpdf, poisspdf, raylpdf, tpdf, tripdf, unidpdf, unifpdf, wblpdf}
## @end deftypefn

function y = pdf (name, x, varargin)

  ## implemented functions
  persistent allDF = { ...
    {"beta"     , "Beta"},                      @betapdf,      2, ...
    {"bino"     , "Binomial"},                  @binopdf,      2, ...
    {"bisa"     , "Birnbaum-Saunders"},         @bisapdf,      2, ...
    {"burr"     , "Burr"},                      @burrpdf,      3, ...
    {"cauchy"   , "Cauchy"},                    @cauchypdf,    2, ...
    {"chi2"     , "Chi-squared"},               @chi2pdf,      1, ...
    {"ev"       , "Extreme Value"},             @evpdf,        2, ...
    {"exp"      , "Exponential"},               @exppdf,       1, ...
    {"f"        , "F-Distribution"},            @fpdf,         2, ...
    {"gam"      , "Gamma"},                     @gampdf,       2, ...
    {"geo"      , "Geometric"},                 @geopdf,       1, ...
    {"gev"      , "Generalized Extreme Value"}, @gevpdf,       3, ...
    {"gp"       , "Generalized Pareto"},        @gppdf,        3, ...
    {"gumbel"   , "Gumbel"},                    @gpcdf,        2, ...
    {"hyge"     , "Hypergeometric"},            @hygepdf,      3, ...
    {"laplace"  , "Laplace"},                   @laplacepdf,   2, ...
    {"logi"     , "Logistic"},                  @logipdf,      2, ...
    {"logl"     , "Log-Logistic"},              @loglpdf,      2, ...
    {"logn"     , "Lognormal"},                 @lognpdf,      2, ...
    {"naka"     , "Nakagami"},                  @nakapdf,      2, ...
    {"nbin"     , "Negative Binomial"},         @nbinpdf,      2, ...
    {"ncf"      , "Noncentral F-Distribution"}, @ncfpdf,       3, ...
    {"nct"      , "Noncentral Student T"},      @nctpdf,       2, ...
    {"ncx2"     , "Noncentral Chi-squared"},    @ncx2pdf,      2, ...
    {"norm"     , "Normal"},                    @normpdf,      2, ...
    {"poiss"    , "Poisson"},                   @poisspdf,     1, ...
    {"rayl"     , "Rayleigh"},                  @raylpdf,      1, ...
    {"t"        , "Student T"},                 @tpdf,         1, ...
    {"tri"      , "Triangular"},                @tripdf,       3, ...
    {"unid"     , "Discrete Uniform"},          @unidpdf,      1, ...
    {"unif"     , "Uniform"},                   @unifpdf,      2, ...
    {"wbl"      , "Weibull"},                   @wblpdf,       2};

  if (! ischar (name))
    error ("pdf: distribution NAME must a char string.");
  endif

  ## Check X being numeric and real
  if (! isnumeric (x))
    error ("pdf: X must be numeric.");
  elseif (! isreal (x))
    error ("pdf: values in X must be real.");
  endif

  ## Get number of arguments
  nargs = numel (varargin);

  ## Get available functions
  pdfnames = allDF(1:3:end);
  pdfhandl = allDF(2:3:end);
  pdf_args = allDF(3:3:end);

  ## Search for PDF function
  idx = cellfun (@(x)any(strcmpi (name, x)), pdfnames);

  if (any (idx))

    if (nargs == pdf_args{idx})
      ## Check that all distribution parameters are numeric
      if (! all (cellfun (@(x)isnumeric(x), (varargin))))
        error ("pdf: distribution parameters must be numeric.");
      endif
      ## Call appropriate iCDF
      y = feval (pdfhandl{idx}, x, varargin{:});

    else
      if (pdf_args{idx} == 1)
        error ("pdf: %s distribution requires 1 parameter.", name);
      else
        error ("pdf: %s distribution requires %d parameters.", ...
               name, pdf_args{idx});
      endif

    endif

  else
    error ("pdf: %s distribution is not implemented in Statistics.", name);
  endif

endfunction

## Test results
%!shared x
%! x = [1:5];
%!assert (pdf ("Beta", x, 5, 2), betapdf (x, 5, 2))
%!assert (pdf ("beta", x, 5, 2), betapdf (x, 5, 2))
%!assert (pdf ("Binomial", x, 5, 2), binopdf (x, 5, 2))
%!assert (pdf ("bino", x, 5, 2), binopdf (x, 5, 2))
%!assert (pdf ("Birnbaum-Saunders", x, 5, 2), bisapdf (x, 5, 2))
%!assert (pdf ("bisa", x, 5, 2), bisapdf (x, 5, 2))
%!assert (pdf ("Burr", x, 5, 2, 2), burrpdf (x, 5, 2, 2))
%!assert (pdf ("burr", x, 5, 2, 2), burrpdf (x, 5, 2, 2))
%!assert (pdf ("Cauchy", x, 5, 2), cauchypdf (x, 5, 2))
%!assert (pdf ("cauchy", x, 5, 2), cauchypdf (x, 5, 2))
%!assert (pdf ("Chi-squared", x, 5), chi2pdf (x, 5))
%!assert (pdf ("chi2", x, 5), chi2pdf (x, 5))
%!assert (pdf ("Extreme Value", x, 5, 2), evpdf (x, 5, 2))
%!assert (pdf ("ev", x, 5, 2), evpdf (x, 5, 2))
%!assert (pdf ("Exponential", x, 5), exppdf (x, 5))
%!assert (pdf ("exp", x, 5), exppdf (x, 5))
%!assert (pdf ("F-Distribution", x, 5, 2), fpdf (x, 5, 2))
%!assert (pdf ("f", x, 5, 2), fpdf (x, 5, 2))
%!assert (pdf ("Gamma", x, 5, 2), gampdf (x, 5, 2))
%!assert (pdf ("gam", x, 5, 2), gampdf (x, 5, 2))
%!assert (pdf ("Geometric", x, 5), geopdf (x, 5))
%!assert (pdf ("geo", x, 5), geopdf (x, 5))
%!assert (pdf ("Generalized Extreme Value", x, 5, 2, 2), gevpdf (x, 5, 2, 2))
%!assert (pdf ("gev", x, 5, 2, 2), gevpdf (x, 5, 2, 2))
%!assert (pdf ("Generalized Pareto", x, 5, 2, 2), gppdf (x, 5, 2, 2))
%!assert (pdf ("gp", x, 5, 2, 2), gppdf (x, 5, 2, 2))
%!assert (pdf ("Gumbel", x, 5, 2), gumbelpdf (x, 5, 2))
%!assert (pdf ("gumbel", x, 5, 2), gumbelpdf (x, 5, 2))
%!assert (pdf ("Hypergeometric", x, 5, 2, 2), hygepdf (x, 5, 2, 2))
%!assert (pdf ("hyge", x, 5, 2, 2), hygepdf (x, 5, 2, 2))
%!assert (pdf ("Laplace", x, 5, 2), laplacepdf (x, 5, 2))
%!assert (pdf ("laplace", x, 5, 2), laplacepdf (x, 5, 2))
%!assert (pdf ("Logistic", x, 5, 2), logipdf (x, 5, 2))
%!assert (pdf ("logi", x, 5, 2), logipdf (x, 5, 2))
%!assert (pdf ("Log-Logistic", x, 5, 2), loglpdf (x, 5, 2))
%!assert (pdf ("logl", x, 5, 2), loglpdf (x, 5, 2))
%!assert (pdf ("Lognormal", x, 5, 2), lognpdf (x, 5, 2))
%!assert (pdf ("logn", x, 5, 2), lognpdf (x, 5, 2))
%!assert (pdf ("Nakagami", x, 5, 2), nakapdf (x, 5, 2))
%!assert (pdf ("naka", x, 5, 2), nakapdf (x, 5, 2))
%!assert (pdf ("Negative Binomial", x, 5, 2), nbinpdf (x, 5, 2))
%!assert (pdf ("nbin", x, 5, 2), nbinpdf (x, 5, 2))
%!assert (pdf ("Noncentral F-Distribution", x, 5, 2, 2), ncfpdf (x, 5, 2, 2))
%!assert (pdf ("ncf", x, 5, 2, 2), ncfpdf (x, 5, 2, 2))
%!assert (pdf ("Noncentral Student T", x, 5, 2), nctpdf (x, 5, 2))
%!assert (pdf ("nct", x, 5, 2), nctpdf (x, 5, 2))
%!assert (pdf ("Noncentral Chi-Squared", x, 5, 2), ncx2pdf (x, 5, 2))
%!assert (pdf ("ncx2", x, 5, 2), ncx2pdf (x, 5, 2))
%!assert (pdf ("Normal", x, 5, 2), normpdf (x, 5, 2))
%!assert (pdf ("norm", x, 5, 2), normpdf (x, 5, 2))
%!assert (pdf ("Poisson", x, 5), poisspdf (x, 5))
%!assert (pdf ("poiss", x, 5), poisspdf (x, 5))
%!assert (pdf ("Rayleigh", x, 5), raylpdf (x, 5))
%!assert (pdf ("rayl", x, 5), raylpdf (x, 5))
%!assert (pdf ("Student T", x, 5), tpdf (x, 5))
%!assert (pdf ("t", x, 5), tpdf (x, 5))
%!assert (pdf ("Triangular", x, 5, 2, 2), tripdf (x, 5, 2, 2))
%!assert (pdf ("tri", x, 5, 2, 2), tripdf (x, 5, 2, 2))
%!assert (pdf ("Discrete Uniform", x, 5), unidpdf (x, 5))
%!assert (pdf ("unid", x, 5), unidpdf (x, 5))
%!assert (pdf ("Uniform", x, 5, 2), unifpdf (x, 5, 2))
%!assert (pdf ("unif", x, 5, 2), unifpdf (x, 5, 2))
%!assert (pdf ("Weibull", x, 5, 2), wblpdf (x, 5, 2))
%!assert (pdf ("wbl", x, 5, 2), wblpdf (x, 5, 2))

## Test input validation
%!error<pdf: distribution NAME must a char string.> pdf (1)
%!error<pdf: distribution NAME must a char string.> pdf ({"beta"})
%!error<pdf: X must be numeric.> pdf ("beta", {[1 2 3 4 5]})
%!error<pdf: X must be numeric.> pdf ("beta", "text")
%!error<pdf: values in X must be real.> pdf ("beta", 1+i)
%!error<pdf: distribution parameters must be numeric.> ...
%! pdf ("Beta", x, "a", 2)
%!error<pdf: distribution parameters must be numeric.> ...
%! pdf ("Beta", x, 5, "")
%!error<pdf: distribution parameters must be numeric.> ...
%! pdf ("Beta", x, 5, {2})
%!error<pdf: chi2 distribution requires 1 parameter.> pdf ("chi2", x)
%!error<pdf: Beta distribution requires 2 parameters.> pdf ("Beta", x, 5)
%!error<pdf: Burr distribution requires 3 parameters.> pdf ("Burr", x, 5)
%!error<pdf: Burr distribution requires 3 parameters.> pdf ("Burr", x, 5, 2)
