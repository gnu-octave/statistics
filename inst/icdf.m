# Copyright (C) 2023 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {statistics} {@var{retval} =} icdf (@var{name}, @var{p}, @dots{})
##
## Return the inverse CDF of @var{name} distribution function for value @var{p}.
##
## This is a wrapper for the @qcode{name}inv and @qcode{name}_inv functions
## available in the statistics package. See the corresponding functions' help to
## learn the signification of the arguments after @var{p}.
##
## @var{name} must be a char string of the name or the abbreviation of the
## desired quantile distribution function as listed in the followng table.
## The last column shows the maximum number of extra arguments that must be
## passed to the desired inverse CDF.
##
## @multitable @columnfractions 0.4 0.05 0.2 0.05 0.3
## @headitem Distribution Name @tab @tab Abbreviation @tab @tab Max Arguments
## @item @qcode{"Birnbaum-Saunders"} @tab @tab @qcode{"bbs"} @tab @tab 3
## @item @qcode{"Beta"} @tab @tab @qcode{"beta"} @tab @tab 2
## @item @qcode{"Binomial"} @tab @tab @qcode{"bino"} @tab @tab 2
## @item @qcode{"Burr"} @tab @tab @qcode{"burr"} @tab @tab 3
## @item @qcode{"Cauchy"} @tab @tab @qcode{"cauchy"} @tab @tab 2
## @item @qcode{"Chi-square"} @tab @tab @qcode{"chi2"} @tab @tab 1
## @item @qcode{"Extreme Value"} @tab @tab @qcode{"ev"} @tab @tab 4
## @item @qcode{"Exponential"} @tab @tab @qcode{"exp"} @tab @tab 1
## @item @qcode{"F-Distribution"} @tab @tab @qcode{"f"} @tab @tab 2
## @item @qcode{"Gamma"} @tab @tab @qcode{"gam"} @tab @tab 2
## @item @qcode{"Geometric"} @tab @tab @qcode{"geo"} @tab @tab 1
## @item @qcode{"Generalized Extreme Value"} @tab @tab @qcode{"gev"} @tab @tab 3
## @item @qcode{"Generalized Pareto"} @tab @tab @qcode{"gp"} @tab @tab 3
## @item @qcode{"Hypergeometric"} @tab @tab @qcode{"hyge"} @tab @tab 3
## @item @qcode{"Laplace"} @tab @tab @qcode{"laplace"} @tab @tab 2
## @item @qcode{"Logistic"} @tab @tab @qcode{"logistic"} @tab @tab 2
## @item @qcode{"Lognormal"} @tab @tab @qcode{"logn"} @tab @tab 2
## @item @qcode{"Nakagami"} @tab @tab @qcode{"naka"} @tab @tab 2
## @item @qcode{"Negative Binomial"} @tab @tab @qcode{"nbin"} @tab @tab 2
## @item @qcode{"Noncentral F-Distribution"} @tab @tab @qcode{"ncf"} @tab @tab 3
## @item @qcode{"Noncentral Student T"} @tab @tab @qcode{"nct"} @tab @tab 2
## @item @qcode{"Noncentral Chi-Square"} @tab @tab @qcode{"ncx2"} @tab @tab 2
## @item @qcode{"Normal"} @tab @tab @qcode{"norm"} @tab @tab 2
## @item @qcode{"Poisson"} @tab @tab @qcode{"poiss"} @tab @tab 1
## @item @qcode{"Rayleigh"} @tab @tab @qcode{"rayl"} @tab @tab 1
## @item @qcode{"Standard Normal"} @tab @tab @qcode{"stdnormal"} @tab @tab 0
## @item @qcode{"Student T"} @tab @tab @qcode{"t"} @tab @tab 1
## @item @qcode{"Triangular"} @tab @tab @qcode{"tri"} @tab @tab 3
## @item @qcode{"Discrete Uniform"} @tab @tab @qcode{"unid"} @tab @tab 1
## @item @qcode{"Uniform"} @tab @tab @qcode{"unif"} @tab @tab 2
## @item @qcode{"Weibull"} @tab @tab @qcode{"wbl"} @tab @tab 2
## @end multitable
##
## @seealso{cdf, pdf, random, bbsinv, betainv, binoinv, burrinv, cauchy_inv,
## chi2inv, evinv, expinv, finv, gaminv, geoinv, gevinv, gpinv, hygeinv,
## laplace_inv, logistic_inv, logninv, nakainv, nbininv, ncfinv, nctinv,
## ncx2inv, norminv, poissinv, raylinv, stdnormal_inv, tinv, triinv, unidinv,
## unifinv, wblinv}
## @end deftypefn

function [retval] = icdf (name, varargin)
  ## implemented functions
  persistent allpdf = { ...
    {"bbs"      , "Birnbaum-Saunders"},         @bbsinv,       3, ...
    {"beta"     , "Beta"},                      @betainv,      2, ...
    {"bino"     , "Binomial"},                  @binoinv,      2, ...
    {"burr"     , "Burr"},                      @burrinv,      3, ...
    {"cauchy"   , "Cauchy"},                    @cauchy_inv,   2, ...
    {"chi2"     , "Chi-square"},                @chi2inv,      1, ...
    {"ev"       , "Extreme Value"},             @evinv,        4, ...
    {"exp"      , "Exponential"},               @expinv,       1, ...
    {"f"        , "F-Distribution"},            @finv,         2, ...
    {"gam"      , "Gamma"},                     @gaminv,       2, ...
    {"geo"      , "Geometric"},                 @geoinv,       1, ...
    {"gev"      , "Generalized Extreme Value"}, @gevinv,       3, ...
    {"gp"       , "Generalized Pareto"},        @gpinv,        3, ...
    {"hyge"     , "Hypergeometric"},            @hygeinv,      3, ...
    {"laplace"  , "Laplace"},                   @laplace_inv,  2, ...
    {"logistic" , "Logistic"},                  @logistic_inv, 2, ...
    {"logn"     , "Lognormal"},                 @logninv,      2, ...
    {"naka"     , "Nakagami"},                  @nakainv,      2, ...
    {"nbin"     , "Negative Binomial"},         @nbininv,      2, ...
    {"ncf"      , "Noncentral F-Distribution"}, @ncfinv,       3, ...
    {"nct"      , "Noncentral Student T"},      @nctinv,       2, ...
    {"ncx2"     , "Noncentral Chi-Square"},     @ncx2inv,      2, ...
    {"norm"     , "Normal"},                    @norminv,      2, ...
    {"poiss"    , "Poisson"},                   @poissinv,     1, ...
    {"rayl"     , "Rayleigh"},                  @raylinv,      1, ...
    {"stdnormal", "Standard Normal"},           @stdnormal_inv,0, ...
    {"t"        , "Student T"},                 @tinv,         1, ...
    {"tri"      , "Triangular"},                @triinv,       3, ...
    {"unid"     , "Discrete Uniform"},          @unidinv,      1, ...
    {"unif"     , "Uniform"},                   @unifinv,      2, ...
    {"wbl"      , "Weibull"},                   @wblinv,       2};

  if (numel (varargin) < 1 || ! ischar (name))
    print_usage ();
  endif

  p = varargin{1};
  varargin(1) = [];
  nargs = numel (varargin);

  icdfnames = allpdf(1:3:end);
  icdfhandl = allpdf(2:3:end);
  icdf_args = allpdf(3:3:end);

  idx = cellfun (@(x)any(strcmpi (name, x)), icdfnames);
  ## Add special list
  special = {"copula", "Copula family"};

  if (any (idx))

    if (nargs > icdf_args{idx})
      if (icdf_args{idx} == 1)
        error ("icdf: %s takes only 1 extra argument.", name);
      else
        error ("icdf: %s takes up to %d extra arguments.", ...
               name, icdf_args{idx});
      endif
    endif

    retval = feval (icdfhandl{idx}, p, varargin{:});

  else
    error ("pdf: %s distribution is not implemented in Statistics.", name);
  endif

endfunction

%!test
%! assert(icdf ("norm", 0.05, 0, 1), norminv (0.05, 0, 1))
%!test
%! p = [0.05, 0.1, 0.2];
%! assert (icdf ("bbs", p, ones (1,3), ones (1,3), zeros (1,3)), ...
%!         bbsinv (p, 1, 1, 0))
%!test
%! p = [0.001 0.05];
%! df = 3;
%! assert (icdf ("Chi-square", p, df), chi2inv (p, df), 1e-10)

%!error icdf ("Birnbaum-Saunders", 1, 1, 2, 3, 4)
%!error icdf ("some", 1, 2)
