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
## @seealso{icdf, pdf, random, betainv, binoinv, bisainv, burrinv, cauchyinv,
## chi2inv, evinv, expinv, finv, gaminv, geoinv, gevinv, gpinv, hygeinv,
## laplaceinv, logiinv, loglinv, logninv, nakainv, nbininv, ncfinv, nctinv,
## ncx2inv, norminv, poissinv, raylinv, tinv, triinv, unidinv, unifinv, wblinv}
## @end deftypefn

function x = icdf (name, p, varargin)

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
    {"gumbel"   , "Gumbel"},                    @gpcdf,        2, ...
    {"hyge"     , "Hypergeometric"},            @hygeinv,      3, ...
    {"laplace"  , "Laplace"},                   @laplaceinv,   2, ...
    {"logi"     , "Logistic"},                  @logiinv,      2, ...
    {"logl"     , "Log-Logistic"},              @loglinv,      2, ...
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

  if (! ischar (name))
    error ("icdf: distribution NAME must a char string.");
  endif

  ## Check P being numeric and real
  if (! isnumeric (p))
    error ("icdf: P must be numeric.");
  elseif (! isreal (p))
    error ("icdf: values in P must be real.");
  endif

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
        error ("icdf: distribution parameters must be numeric.");
      endif
      ## Call appropriate iCDF
      x = feval (icdfhandl{idx}, p, varargin{:});

    else
      if (icdf_args{idx} == 1)
        error ("icdf: %s distribution requires 1 parameter.", name);
      else
        error ("icdf: %s distribution requires %d parameters.", ...
               name, icdf_args{idx});
      endif

    endif

  else
    error ("icdf: %s distribution is not implemented in Statistics.", name);
  endif

endfunction

## Test results
%!shared p
%! p = [0.05:0.05:0.5];
%!assert (icdf ("Beta", p, 5, 2), betainv (p, 5, 2))
%!assert (icdf ("beta", p, 5, 2), betainv (p, 5, 2))
%!assert (icdf ("Binomial", p, 5, 2), binoinv (p, 5, 2))
%!assert (icdf ("bino", p, 5, 2), binoinv (p, 5, 2))
%!assert (icdf ("Birnbaum-Saunders", p, 5, 2), bisainv (p, 5, 2))
%!assert (icdf ("bisa", p, 5, 2), bisainv (p, 5, 2))
%!assert (icdf ("Burr", p, 5, 2, 2), burrinv (p, 5, 2, 2))
%!assert (icdf ("burr", p, 5, 2, 2), burrinv (p, 5, 2, 2))
%!assert (icdf ("Cauchy", p, 5, 2), cauchyinv (p, 5, 2))
%!assert (icdf ("cauchy", p, 5, 2), cauchyinv (p, 5, 2))
%!assert (icdf ("Chi-squared", p, 5), chi2inv (p, 5))
%!assert (icdf ("chi2", p, 5), chi2inv (p, 5))
%!assert (icdf ("Extreme Value", p, 5, 2), evinv (p, 5, 2))
%!assert (icdf ("ev", p, 5, 2), evinv (p, 5, 2))
%!assert (icdf ("Exponential", p, 5), expinv (p, 5))
%!assert (icdf ("exp", p, 5), expinv (p, 5))
%!assert (icdf ("F-Distribution", p, 5, 2), finv (p, 5, 2))
%!assert (icdf ("f", p, 5, 2), finv (p, 5, 2))
%!assert (icdf ("Gamma", p, 5, 2), gaminv (p, 5, 2))
%!assert (icdf ("gam", p, 5, 2), gaminv (p, 5, 2))
%!assert (icdf ("Geometric", p, 5), geoinv (p, 5))
%!assert (icdf ("geo", p, 5), geoinv (p, 5))
%!assert (icdf ("Generalized Extreme Value", p, 5, 2, 2), gevinv (p, 5, 2, 2))
%!assert (icdf ("gev", p, 5, 2, 2), gevinv (p, 5, 2, 2))
%!assert (icdf ("Generalized Pareto", p, 5, 2, 2), gpinv (p, 5, 2, 2))
%!assert (icdf ("gp", p, 5, 2, 2), gpinv (p, 5, 2, 2))
%!assert (icdf ("Gumbel", p, 5, 2), gumbelinv (p, 5, 2))
%!assert (icdf ("gumbel", p, 5, 2), gumbelinv (p, 5, 2))
%!assert (icdf ("Hypergeometric", p, 5, 2, 2), hygeinv (p, 5, 2, 2))
%!assert (icdf ("hyge", p, 5, 2, 2), hygeinv (p, 5, 2, 2))
%!assert (icdf ("Laplace", p, 5, 2), laplaceinv (p, 5, 2))
%!assert (icdf ("laplace", p, 5, 2), laplaceinv (p, 5, 2))
%!assert (icdf ("Logistic", p, 5, 2), logiinv (p, 5, 2))
%!assert (icdf ("logi", p, 5, 2), logiinv (p, 5, 2))
%!assert (icdf ("Log-Logistic", p, 5, 2), loglinv (p, 5, 2))
%!assert (icdf ("logl", p, 5, 2), loglinv (p, 5, 2))
%!assert (icdf ("Lognormal", p, 5, 2), logninv (p, 5, 2))
%!assert (icdf ("logn", p, 5, 2), logninv (p, 5, 2))
%!assert (icdf ("Nakagami", p, 5, 2), nakainv (p, 5, 2))
%!assert (icdf ("naka", p, 5, 2), nakainv (p, 5, 2))
%!assert (icdf ("Negative Binomial", p, 5, 2), nbininv (p, 5, 2))
%!assert (icdf ("nbin", p, 5, 2), nbininv (p, 5, 2))
%!assert (icdf ("Noncentral F-Distribution", p, 5, 2, 2), ncfinv (p, 5, 2, 2))
%!assert (icdf ("ncf", p, 5, 2, 2), ncfinv (p, 5, 2, 2))
%!assert (icdf ("Noncentral Student T", p, 5, 2), nctinv (p, 5, 2))
%!assert (icdf ("nct", p, 5, 2), nctinv (p, 5, 2))
%!assert (icdf ("Noncentral Chi-Squared", p, 5, 2), ncx2inv (p, 5, 2))
%!assert (icdf ("ncx2", p, 5, 2), ncx2inv (p, 5, 2))
%!assert (icdf ("Normal", p, 5, 2), norminv (p, 5, 2))
%!assert (icdf ("norm", p, 5, 2), norminv (p, 5, 2))
%!assert (icdf ("Poisson", p, 5), poissinv (p, 5))
%!assert (icdf ("poiss", p, 5), poissinv (p, 5))
%!assert (icdf ("Rayleigh", p, 5), raylinv (p, 5))
%!assert (icdf ("rayl", p, 5), raylinv (p, 5))
%!assert (icdf ("Student T", p, 5), tinv (p, 5))
%!assert (icdf ("t", p, 5), tinv (p, 5))
%!assert (icdf ("Triangular", p, 5, 2, 2), triinv (p, 5, 2, 2))
%!assert (icdf ("tri", p, 5, 2, 2), triinv (p, 5, 2, 2))
%!assert (icdf ("Discrete Uniform", p, 5), unidinv (p, 5))
%!assert (icdf ("unid", p, 5), unidinv (p, 5))
%!assert (icdf ("Uniform", p, 5, 2), unifinv (p, 5, 2))
%!assert (icdf ("unif", p, 5, 2), unifinv (p, 5, 2))
%!assert (icdf ("Weibull", p, 5, 2), wblinv (p, 5, 2))
%!assert (icdf ("wbl", p, 5, 2), wblinv (p, 5, 2))

## Test input validation
%!error<icdf: distribution NAME must a char string.> icdf (1)
%!error<icdf: distribution NAME must a char string.> icdf ({"beta"})
%!error<icdf: P must be numeric.> icdf ("beta", {[1 2 3 4 5]})
%!error<icdf: P must be numeric.> icdf ("beta", "text")
%!error<icdf: values in P must be real.> icdf ("beta", 1+i)
%!error<icdf: distribution parameters must be numeric.> ...
%! icdf ("Beta", p, "a", 2)
%!error<icdf: distribution parameters must be numeric.> ...
%! icdf ("Beta", p, 5, "")
%!error<icdf: distribution parameters must be numeric.> ...
%! icdf ("Beta", p, 5, {2})
%!error<icdf: chi2 distribution requires 1 parameter.> icdf ("chi2", p)
%!error<icdf: Beta distribution requires 2 parameters.> icdf ("Beta", p, 5)
%!error<icdf: Burr distribution requires 3 parameters.> icdf ("Burr", p, 5)
%!error<icdf: Burr distribution requires 3 parameters.> icdf ("Burr", p, 5, 2)
