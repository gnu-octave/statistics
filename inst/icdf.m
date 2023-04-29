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
## @deftypefn  {statistics} {@var{x} =} icdf (@var{name}, @var{p}, @var{A})
## @deftypefnx {statistics} {@var{x} =} icdf (@var{name}, @var{p}, @var{A}, @var{B})
## @deftypefnx {statistics} {@var{x} =} icdf (@var{name}, @var{p}, @var{A}, @var{B}, @var{C})
##
## Return the inverse CDF of a univariate distribution evaluated at @var{p}.
##
## @code{icdf} is a wrapper for the univariate quantile distribution functions
## (iCDF) available in the statistics package.  See the corresponding functions'
## help to learn the signification of the parameters after @var{p}.
##
## @code{@var{x} = icdf (@var{name}, @var{p}, @var{A})} returns the iCDF for the
## one-parameter distribution family specified by @var{name} and the
## distribution parameter @var{A}, evaluated at the values in @var{p}.
##
## @code{@var{x} = icdf (@var{name}, @var{p}, @var{A}, @var{B})} returns the
## iCDF for the two-parameter distribution family specified by @var{name} and
## the distribution parameters @var{A} and @var{B}, evaluated at the values in
## @var{p}.
##
## @code{@var{x} = icdf (@var{name}, @var{p}, @var{A}, @var{B}, @var{C})}
## returns the iCDF for the three-parameter distribution family specified by
## @var{name} and the distribution parameters @var{A}, @var{B}, and @var{C},
## evaluated at the values in @var{p}.
##
## @var{name} must be a char string of the name or the abbreviation of the
## desired quantile distribution function as listed in the followng table.
## The last column shows the number of required parameters that should be parsed
## after @var{x} to the desired iCDF.
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
## @seealso{cdf, pdf, random, betainv, binoinv, bisainv, burrinv, cauchyinv,
## chi2inv, evinv, expinv, finv, gaminv, geoinv, gevinv, gpinv, hygeinv,
## laplaceinv, logiinv, loglinv, logninv, nakainv, nbininv, ncfinv, nctinv,
## ncx2inv, norminv, poissinv, raylinv, tinv, triinv, unidinv, unifinv, wblinv}
## @end deftypefn

function x = icdf (name, varargin)

  ## implemented functions
  persistent allDF = { ...
    {"beta"     , "Beta"},                      @betainv,      2, ...
    {"bino"     , "Binomial"},                  @binoinv,      2, ...
    {"bisa"     , "Birnbaum-Saunders"},         @bisainv,      2, ...
    {"burr"     , "Burr"},                      @burrinv,      3, ...
    {"cauchy"   , "Cauchy"},                    @cauchyinv,    2, ...
    {"chi2"     , "Chi-squared"},               @chi2inv,      1, ...
    {"ev"       , "Extreme Value"},             @evinv,        2, ...
    {"exp"      , "Exponential"},               @expinv,       1, ...
    {"f"        , "F-Distribution"},            @finv,         2, ...
    {"gam"      , "Gamma"},                     @gaminv,       2, ...
    {"geo"      , "Geometric"},                 @geoinv,       1, ...
    {"gev"      , "Generalized Extreme Value"}, @gevinv,       3, ...
    {"gp"       , "Generalized Pareto"},        @gpinv,        3, ...
    {"hyge"     , "Hypergeometric"},            @hygeinv,      3, ...
    {"laplace"  , "Laplace"},                   @laplaceinv,   2, ...
    {"logi"     , "Logistic"},                  @logiinv,      2, ...
    {"logl"     , "Log-Logistic"},              @logiinv,      2, ...
    {"logn"     , "Lognormal"},                 @logninv,      2, ...
    {"naka"     , "Nakagami"},                  @nakainv,      2, ...
    {"nbin"     , "Negative Binomial"},         @nbininv,      2, ...
    {"ncf"      , "Noncentral F-Distribution"}, @ncfinv,       3, ...
    {"nct"      , "Noncentral Student T"},      @nctinv,       2, ...
    {"ncx2"     , "Noncentral Chi-squared"},    @ncx2inv,      2, ...
    {"norm"     , "Normal"},                    @norminv,      2, ...
    {"poiss"    , "Poisson"},                   @poissinv,     1, ...
    {"rayl"     , "Rayleigh"},                  @raylinv,      1, ...
    {"t"        , "Student T"},                 @tinv,         1, ...
    {"tri"      , "Triangular"},                @triinv,       3, ...
    {"unid"     , "Discrete Uniform"},          @unidinv,      1, ...
    {"unif"     , "Uniform"},                   @unifinv,      2, ...
    {"wbl"      , "Weibull"},                   @wblinv,       2};

  if (numel (varargin) < 1 || ! ischar (name))
    print_usage ();
  endif

  ## Get p-values
  p = varargin{1};
  varargin(1) = [];

  ## Get number of arguments
  nargs = numel (varargin);

  ## Get available functions
  icdfnames = allDF(1:3:end);
  icdfhandl = allDF(2:3:end);
  icdf_args = allDF(3:3:end);

  ## Search for iCDF function
  idx = cellfun (@(x)any(strcmpi (name, x)), icdfnames);

  if (any (idx))

    if (nargs == icdf_args{idx})
      ## Check that all distribution parameters are numeric
      if (! all (cellfun (@(x)isnumeric(x), (varargin))))
        error ("cdf: distribution parameters must be numeric.");
      endif
      ## Call appropriate iCDF
      x = feval (icdfhandl{idx}, p, varargin{:});

    else
      if (icdf_args{idx} == 1)
        error ("icdf: %s requires 1 parameter.", name);
      else
        error ("icdf: %s requires %d parameters.", name, icdf_args{idx});
      endif

    endif

  else
    error ("icdf: %s distribution is not implemented in Statistics.", name);
  endif

endfunction

## Test results
%!test
%! assert(icdf ("norm", 0.05, 0, 1), norminv (0.05, 0, 1))
