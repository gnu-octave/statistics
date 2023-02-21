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
## @deftypefn  {statistics} {@var{retval} =} pdf (@var{name}, @var{x}, @dots{})
##
## Return the PDF of @var{name} distribution function for value @var{x}.
##
## This is a wrapper around various @qcode{name}pdf and @qcode{name}_pdf
## functions. See the corresponding functions' help to learn the signification
## of the arguments after @var{x}.
##
## @var{name} must be a char string of the name or the abbreviation of the
## desired probability distribution function as listed in the followng table.
## The last column shows the maximum number of extra arguments that can be
## passed to the desired PDF.
##
## @multitable @columnfractions 0.45 0.2 0.35
## @headitem Distribution Name @tab Abbreviation @tab Max. Extra Arguments
## @item @qcode{"Birnbaum-Saunders"} @tab @qcode{"bbs"} @tab 3
## @item @qcode{"Beta"} @tab @qcode{"beta"} @tab 2
## @item @qcode{"Binomial"} @tab @qcode{"bino"} @tab 2
## @item @qcode{"Burr"} @tab @qcode{"burr"} @tab 3
## @item @qcode{"Bivariate Normal"} @tab @qcode{"bvn"} @tab 2
## @item @qcode{"Cauchy"} @tab @qcode{"cauchy"} @tab 2
## @item @qcode{"Chi-square"} @tab @qcode{"chi2"} @tab 1
## @item @qcode{"Copula Family"} @tab @qcode{"copula"} @tab 2
## @item @qcode{"Extreme Value"} @tab @qcode{"ev"} @tab 2
## @item @qcode{"Exponential"} @tab @qcode{"exp"} @tab 1
## @item @qcode{"F-Distribution"} @tab @qcode{"f"} @tab 2
## @item @qcode{"Gamma"} @tab @qcode{"gam"} @tab 2
## @item @qcode{"Geometric"} @tab @qcode{"geo"} @tab 1
## @item @qcode{"Generalized Extreme Value"} @tab @qcode{"gev"} @tab 3
## @item @qcode{"Generalized Pareto"} @tab @qcode{"gp"} @tab 3
## @item @qcode{"Hypergeometric"} @tab @qcode{"hyge"} @tab 4
## @item @qcode{"Inverse Wishart"} @tab @qcode{"iwish"} @tab 3
## @item @qcode{"Johnson SU"} @tab @qcode{"jsu"} @tab 2
## @item @qcode{"Laplace"} @tab @qcode{"laplace"} @tab 2
## @item @qcode{"Logistic"} @tab @qcode{"logistic"} @tab 2
## @item @qcode{"Lognormal"} @tab @qcode{"logn"} @tab 2
## @item @qcode{"Multinomial"} @tab @qcode{"mn"} @tab 1
## @item @qcode{"Multivariate Normal"} @tab @qcode{"mvn"} @tab 2
## @item @qcode{"Multivariate Student T"} @tab @qcode{"mvt"} @tab 2
## @item @qcode{"Nakagami"} @tab @qcode{"naka"} @tab 2
## @item @qcode{"Negative Binomial"} @tab @qcode{"nbin"} @tab 2
## @item @qcode{"Noncentral F-Distribution"} @tab @qcode{"ncf"} @tab 3
## @item @qcode{"Noncentral Student T"} @tab @qcode{"nct"} @tab 2
## @item @qcode{"Noncentral Chi-Square"} @tab @qcode{"ncx2"} @tab 2
## @item @qcode{"Normal"} @tab @qcode{"norm"} @tab 2
## @item @qcode{"Poisson"} @tab @qcode{"poiss"} @tab 1
## @item @qcode{"Rayleigh"} @tab @qcode{"rayl"} @tab 1
## @item @qcode{"Standard Normal"} @tab @qcode{"stdnormal"} @tab 0
## @item @qcode{"Student T"} @tab @qcode{"t"} @tab 1
## @item @qcode{"Triangular"} @tab @qcode{"tri"} @tab 3
## @item @qcode{"Discrete Uniform"} @tab @qcode{"unid"} @tab 1
## @item @qcode{"Uniform"} @tab @qcode{"unif"} @tab 2
## @item @qcode{"Von Mises"} @tab @qcode{"vm"} @tab 2
## @item @qcode{"Weibull"} @tab @qcode{"wbl"} @tab 2
## @item @qcode{"Wishart"} @tab @qcode{"wish"} @tab 3
## @end multitable
##
## @seealso{cdf, icdf, random, bbspdf, betapdf, binopdf, burrpdf, bvnpdf,
## cauchy_pdf, chi2pdf, copulapdf, evpdf, exppdf, fpdf, gampdf, geopdf, gevpdf,
## gppdf, hygepdf, iwishpdf, jsupdf, laplace_pdf, logistic_pdf, lognpdf, mvnpdf,
## mvtpdf, nakapdf, nbinpdf, ncfpdf, nctpdf, ncx2pdf, normpdf, poisspdf,
## raylpdf, stdnormal_pdf, tpdf, tripdf, unidpdf, unifpdf, vmpdf, wblpdf,
## wishpdf}
## @end deftypefn

function [retval] = pdf (name, varargin)
  ## implemented functions
  persistent allpdf = { ...
    {"bbs"      , "Birnbaum-Saunders"},         @bbspdf,       3, ...
    {"beta"     , "Beta"},                      @betapdf,      2, ...
    {"bino"     , "Binomial"},                  @binopdf,      2, ...
    {"burr"     , "Burr"},                      @burrpdf,      3, ...
    {"bvn"      , "Bivariate Normal"},          @bvnpdf,       2, ...
    {"cauchy"   , "Cauchy"},                    @cauchy_pdf,   2, ...
    {"chi2"     , "Chi-square"},                @chi2pdf,      1, ...
    {"copula"   , "Copula Family"},             @copulapdf,    2, ...
    {"ev"       , "Extreme Value"},             @evpdf,        2, ...
    {"exp"      , "Exponential"},               @exppdf,       1, ...
    {"f"        , "F-Distribution"},            @fpdf,         2, ...
    {"gam"      , "Gamma"},                     @gampdf,       2, ...
    {"geo"      , "Geometric"},                 @geopdf,       1, ...
    {"gev"      , "Generalized Extreme Value"}, @gevpdf,       3, ...
    {"gp"       , "Generalized Pareto"},        @gppdf,        3, ...
    {"hyge"     , "Hypergeometric"},            @hygepdf,      4, ...
    {"iwish"    , "Inverse Wishart"},           @iwishpdf,     3, ...
    {"jsu"      , "Johnson SU"},                @jsupdf,       2, ...
    {"laplace"  , "Laplace"},                   @laplace_pdf,  2, ...
    {"logistic" , "Logistic"},                  @logistic_pdf, 2, ...
    {"logn"     , "Lognormal"},                 @lognpdf,      2, ...
    {"mvn"      , "Multivariate Normal"},       @mvnpdf,       2, ...
    {"mvt"      , "Multivariate Student T"},    @mvtpdf,       2, ...
    {"naka"     , "Nakagami"},                  @nakapdf,      2, ...
    {"nbin"     , "Negative Binomial"},         @nbinpdf,      2, ...
    {"ncf"      , "Noncentral F-Distribution"}, @ncfpdf,       3, ...
    {"nct"      , "Noncentral Student T"},      @nctpdf,       2, ...
    {"ncx2"     , "Noncentral Chi-Square"},     @ncx2pdf,      2, ...
    {"norm"     , "Normal"},                    @normpdf,      2, ...
    {"poiss"    , "Poisson"},                   @poisspdf,     1, ...
    {"rayl"     , "Rayleigh"},                  @raylpdf,      1, ...
    {"stdnormal", "Standard Normal"},           @stdnormal_pdf,0, ...
    {"t"        , "Student T"},                 @tpdf,         1, ...
    {"tri"      , "Triangular"},                @tripdf,       3, ...
    {"unid"     , "Discrete Uniform"},          @unidpdf,      1, ...
    {"unif"     , "Uniform"},                   @unifpdf,      2, ...
    {"vm"       , "Von Mises"},                 @vmpdf,        2, ...
    {"wbl"      , "Weibull"},                   @wblpdf,       2, ...
    {"wish"     , "Wishart"},                   @wishpdf,      3};

  if (numel (varargin) < 1 || ! ischar (name))
    print_usage ();
  endif

  x = varargin{1};
  varargin(1) = [];
  nargs = numel (varargin);

  pdfnames = allpdf(1:3:end);
  pdfhandl = allpdf(2:3:end);
  pdf_args = allpdf(3:3:end);

  idx = cellfun (@(x)any(strcmpi (name, x)), pdfnames);
  ## Add special list
  special = {"copula", "Copula family"};

  if (any (idx))

    if (nargs > pdf_args{idx})
      if (pdf_args{idx} == 1)
        error ("pdf: %s takes only 1 extra argument.", name);
      else
        error ("pdf: %s takes up to %d extra arguments.", name, pdf_args{idx});
      endif
    endif

    if (! any (strcmpi (name, special)))
      retval = feval (pdfhandl{idx}, x, varargin{:});

    else
      retval = feval (pdfhandl{idx}, varargin{1}, x, varargin{2:end});

    endif

  else
    error ("pdf: %s distribution is not implemented in Statistics.", name);
  endif

endfunction

%!test
%! assert(pdf ("norm", 1, 0, 1), normpdf (1, 0, 1))
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

%!error pdf ("Birnbaum-Saunders", 1, 1, 2, 3, 4)
%!error pdf ("some", 1, 2)
