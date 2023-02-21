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
## This is a wrapper around various @qcode{name}inv and @qcode{name}_inv
## functions. See the corresponding functions' help to learn the signification
## of the arguments after @var{p}.
##
## @var{name} must be a char string of the name or the abbreviation of the
## desired quantile distribution function as listed in the followng table.
## The last column shows the maximum number of extra arguments that must be
## passed to the desired inverse CDF.
##
## @multitable @columnfractions 0.45 0.2 0.35
## @headitem Distribution Name @tab Abbreviation @tab Max. Extra Arguments
## @item @qcode{"Birnbaum-Saunders"} @tab @qcode{"bbs"} @tab 3
## @item @qcode{"Beta"} @tab @qcode{"beta"} @tab 2
## @item @qcode{"Binomial"} @tab @qcode{"bino"} @tab 2
## @item @qcode{"Burr"} @tab @qcode{"burr"} @tab 3
## @item @qcode{"Cauchy"} @tab @qcode{"cauchy"} @tab 2
## @item @qcode{"Chi-square"} @tab @qcode{"chi2"} @tab 1
## @item @qcode{"Extreme Value"} @tab @qcode{"ev"} @tab 4
## @item @qcode{"Exponential"} @tab @qcode{"exp"} @tab 1
## @item @qcode{"F-Distribution"} @tab @qcode{"f"} @tab 2
## @item @qcode{"Gamma"} @tab @qcode{"gam"} @tab 2
## @item @qcode{"Geometric"} @tab @qcode{"geo"} @tab 1
## @item @qcode{"Generalized Extreme Value"} @tab @qcode{"gev"} @tab 3
## @item @qcode{"Generalized Pareto"} @tab @qcode{"gp"} @tab 3
## @item @qcode{"Hypergeometric"} @tab @qcode{"hyge"} @tab 3
## @item @qcode{"Laplace"} @tab @qcode{"laplace"} @tab 2
## @item @qcode{"Logistic"} @tab @qcode{"logistic"} @tab 2
## @item @qcode{"Lognormal"} @tab @qcode{"logn"} @tab 2
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
## @item @qcode{"Weibull"} @tab @qcode{"wbl"} @tab 2
## @end multitable
##
## @seealso{cdf, pdf, random, bbsinv, betainv, binoinv, burrinv, cauchy_inv,
## chi2inv, evinv, expinv, finv, gaminv, geoinv, gevinv, gpinv, hygeinv,
## laplace_inv, logistic_inv, logninv, nakainv, nbininv, ncfinv, nctinv,
## ncx2inv, norminv, poissinv, raylinv, stdnormal_inv, tinv, triinv, unidinv,
## unifinv, wblinv}
## @end deftypefn

function [retval] = pdf (name, varargin)
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
