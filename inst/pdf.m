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

function y = pdf (name, varargin)

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
    {"hyge"     , "Hypergeometric"},            @hygepdf,      3, ...
    {"laplace"  , "Laplace"},                   @laplacepdf,   2, ...
    {"logi"     , "Logistic"},                  @logipdf,      2, ...
    {"logl"     , "Log-Logistic"},              @logipdf,      2, ...
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

  if (numel (varargin) < 1 || ! ischar (name))
    print_usage ();
  endif

  ## Get data
  x = varargin{1};
  varargin(1) = [];

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
      if (cdf_args{idx} == 1)
        error ("pdf: %s requires 1 parameter.", name);
      else
        error ("pdf: %s requires %d parameters.", name, pdf_args{idx});
      endif

    endif

  else
    error ("pdf: %s distribution is not implemented in Statistics.", name);
  endif

endfunction

## Test results
%!test
%! assert(pdf ("norm", 1, 0, 1), normpdf (1, 0, 1))
