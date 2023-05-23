## Copyright (C) 2013 Pantxo Diribarne
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
## @deftypefn  {statistics} {@var{p} =} cdf (@var{name}, @var{x}, @var{A})
## @deftypefnx {statistics} {@var{p} =} cdf (@var{name}, @var{x}, @var{A}, @var{B})
## @deftypefnx {statistics} {@var{p} =} cdf (@var{name}, @var{x}, @var{A}, @var{B}, @var{C})
## @deftypefnx {statistics} {@var{p} =} cdf (@dots{}, @qcode{"upper"})
##
## Return the CDF of a univariate distribution evaluated at @var{x}.
##
## @code{cdf} is a wrapper for the univariate cumulative distribution functions
## available in the statistics package.  See the corresponding functions' help
## to learn the signification of the parameters after @var{x}.
##
## @code{@var{p} = cdf (@var{name}, @var{x}, @var{A})} returns the CDF for the
## one-parameter distribution family specified by @var{name} and the
## distribution parameter @var{A}, evaluated at the values in @var{x}.
##
## @code{@var{p} = cdf (@var{name}, @var{x}, @var{A}, @var{B})} returns the CDF
## for the two-parameter distribution family specified by @var{name} and the
## distribution parameters @var{A} and @var{B}, evaluated at the values in
## @var{x}.
##
## @code{@var{p} = cdf (@var{name}, @var{x}, @var{A}, @var{B}, @var{C})} returns
## the CDF for the three-parameter distribution family specified by @var{name}
## and the distribution parameters @var{A}, @var{B}, and @var{C}, evaluated at
## the values in @var{x}.
##
## @code{@var{p} = cdf (@dots{}, @qcode{"upper"})} returns the complement of the
## CDF using an algorithm that more accurately computes the extreme upper-tail
## probabilities.  @qcode{"upper"} can follow any of the input arguments in the
## previous syntaxes.
##
## @var{name} must be a char string of the name or the abbreviation of the
## desired cumulative distribution function as listed in the followng table.
## The last column shows the number of required parameters that should be parsed
## after @var{x} to the desired CDF.  The optional input argument
## @qcode{"upper"} does not count in the required number of parameters.
##
## @multitable @columnfractions 0.4 0.05 0.2 0.05 0.3
## @headitem Distribution Name @tab @tab Abbreviation @tab @tab Input Parameters
## @item @qcode{"Beta"} @tab @tab @qcode{"beta"} @tab @tab 2
## @item @qcode{"Binomial"} @tab @tab @qcode{"bino"} @tab @tab 2
## @item @qcode{"Birnbaum-Saunders"} @tab @tab @qcode{"bisa"} @tab @tab 2
## @item @qcode{"Burr"} @tab  @tab @qcode{"burr"} @tab  @tab 3
## @item @qcode{"Cauchy"} @tab @tab @qcode{"cauchy"} @tab @tab 2
## @item @qcode{"Chi-squared"} @tab @tab @qcode{"chi2"} @tab @tab 1
## @item @qcode{"Extreme Value"} @tab @tab @qcode{"ev"} @tab @tab 2
## @item @qcode{"Exponential"} @tab @tab @qcode{"exp"} @tab @tab 1
## @item @qcode{"F-Distribution"} @tab @tab @qcode{"f"} @tab @tab 2
## @item @qcode{"Gamma"} @tab @tab @qcode{"gam"} @tab @tab 2
## @item @qcode{"Geometric"} @tab @tab @qcode{"geo"} @tab @tab 1
## @item @qcode{"Generalized Extreme Value"} @tab @tab @qcode{"gev"} @tab @tab 3
## @item @qcode{"Generalized Pareto"} @tab @tab @qcode{"gp"} @tab @tab 3
## @item @qcode{"Gumbel"} @tab @tab @qcode{"gumbel"} @tab @tab 2
## @item @qcode{"Half-normal"} @tab @tab @qcode{"hn"} @tab @tab 2
## @item @qcode{"Hypergeometric"} @tab @tab @qcode{"hyge"} @tab @tab 3
## @item @qcode{"Inverse Gaussian"} @tab @tab @qcode{"invg"} @tab @tab 2
## @item @qcode{"Laplace"} @tab @tab @qcode{"laplace"} @tab @tab 2
## @item @qcode{"Logistic"} @tab @tab @qcode{"logi"} @tab @tab 2
## @item @qcode{"Log-Logistic"} @tab @tab @qcode{"logl"} @tab @tab 2
## @item @qcode{"Lognormal"} @tab @tab @qcode{"logn"} @tab @tab 2
## @item @qcode{"Nakagami"} @tab @tab @qcode{"naka"} @tab @tab 2
## @item @qcode{"Negative Binomial"} @tab @tab @qcode{"nbin"} @tab @tab 2
## @item @qcode{"Noncentral F-Distribution"} @tab @tab @qcode{"ncf"} @tab @tab 3
## @item @qcode{"Noncentral Student T"} @tab @tab @qcode{"nct"} @tab @tab 2
## @item @qcode{"Noncentral Chi-Squared"} @tab @tab @qcode{"ncx2"} @tab @tab 2
## @item @qcode{"Normal"} @tab @tab @qcode{"norm"} @tab @tab 2
## @item @qcode{"Poisson"} @tab @tab @qcode{"poiss"} @tab @tab 1
## @item @qcode{"Rayleigh"} @tab @tab @qcode{"rayl"} @tab @tab 1
## @item @qcode{"Student T"} @tab @tab @qcode{"t"} @tab @tab 1
## @item @qcode{"Triangular"} @tab @tab @qcode{"tri"} @tab @tab 3
## @item @qcode{"Discrete Uniform"} @tab @tab @qcode{"unid"} @tab @tab 1
## @item @qcode{"Uniform"} @tab @tab @qcode{"unif"} @tab @tab 2
## @item @qcode{"Von Mises"} @tab @tab @qcode{"vm"} @tab @tab 2
## @item @qcode{"Weibull"} @tab @tab @qcode{"wbl"} @tab @tab 2
## @end multitable
##
## @seealso{icdf, pdf, cdf, betacdf, binocdf, bisacdf, burrcdf, cauchycdf,
## chi2cdf, evcdf, expcdf, fcdf, gamcdf, geocdf, gevcdf, gpcdf, gumbelcdf,
## hygecdf, laplacecdf, logicdf, loglcdf, logncdf, nakacdf, nbincdf, ncfcdf,
## nctcdf, ncx2cdf, normcdf, poisscdf, raylcdf, tcdf, tricdf, unidcdf, unifcdf,
## wblcdf}
## @end deftypefn

function p = cdf (name, x, varargin)

  ## implemented functions
  persistent allDF = { ...
    {"beta"     , "Beta"},                      @betacdf,      2, ...
    {"bino"     , "Binomial"},                  @binocdf,      2, ...
    {"bisa"     , "Birnbaum-Saunders"},         @bisacdf,      2, ...
    {"burr"     , "Burr"},                      @burrcdf,      3, ...
    {"cauchy"   , "Cauchy"},                    @cauchycdf,    2, ...
    {"chi2"     , "Chi-squared"},               @chi2cdf,      1, ...
    {"ev"       , "Extreme Value"},             @evcdf,        2, ...
    {"exp"      , "Exponential"},               @expcdf,       1, ...
    {"f"        , "F-Distribution"},            @fcdf,         2, ...
    {"gam"      , "Gamma"},                     @gamcdf,       2, ...
    {"geo"      , "Geometric"},                 @geocdf,       1, ...
    {"gev"      , "Generalized Extreme Value"}, @gevcdf,       3, ...
    {"gp"       , "Generalized Pareto"},        @gpcdf,        3, ...
    {"gumbel"   , "Gumbel"},                    @gumbelcdf,    2, ...
    {"hn"       , "Half-normal"},               @hncdf,        2, ...
    {"hyge"     , "Hypergeometric"},            @hygecdf,      3, ...
    {"invg"     , "Inverse Gaussian"},          @invgcdf,      2, ...
    {"laplace"  , "Laplace"},                   @laplacecdf,   2, ...
    {"logi"     , "Logistic"},                  @logicdf,      2, ...
    {"logl"     , "Log-Logistic"},              @loglcdf,      2, ...
    {"logn"     , "Lognormal"},                 @logncdf,      2, ...
    {"naka"     , "Nakagami"},                  @nakacdf,      2, ...
    {"nbin"     , "Negative Binomial"},         @nbincdf,      2, ...
    {"ncf"      , "Noncentral F-Distribution"}, @ncfcdf,       3, ...
    {"nct"      , "Noncentral Student T"},      @nctcdf,       2, ...
    {"ncx2"     , "Noncentral Chi-squared"},    @ncx2cdf,      2, ...
    {"norm"     , "Normal"},                    @normcdf,      2, ...
    {"poiss"    , "Poisson"},                   @poisscdf,     1, ...
    {"rayl"     , "Rayleigh"},                  @raylcdf,      1, ...
    {"t"        , "Student T"},                 @tcdf,         1, ...
    {"tri"      , "Triangular"},                @tricdf,       3, ...
    {"unid"     , "Discrete Uniform"},          @unidcdf,      1, ...
    {"unif"     , "Uniform"},                   @unifcdf,      2, ...
    {"vm"       , "Von Mises"},                 @vmcdf,        2, ...
    {"wbl"      , "Weibull"},                   @wblcdf,       2};

  ## Check NAME being a char string
  if (! ischar (name))
    error ("cdf: distribution NAME must a char string.");
  endif

  ## Check X being numeric and real
  if (! isnumeric (x))
    error ("cdf: X must be numeric.");
  elseif (! isreal (x))
    error ("cdf: values in X must be real.");
  endif

  ## Get number of arguments
  nargs = numel (varargin);

  ## Get available functions
  cdfnames = allDF(1:3:end);
  cdfhandl = allDF(2:3:end);
  cdf_args = allDF(3:3:end);

  ## Search for CDF function
  idx = cellfun (@(x)any(strcmpi (name, x)), cdfnames);

  if (any (idx))

    if (nargs == cdf_args{idx} + 1)
      ## Check for "upper" option
      if (! strcmpi (varargin{nargs}, "upper"))
        error ("cdf: invalid argument for upper tail.");
      else
        ## Check that all remaining distribution parameters are numeric
        if (! all (cellfun (@(x)isnumeric(x), (varargin([1:nargs-1])))))
          error ("cdf: distribution parameters must be numeric.");
        endif
        ## Call appropriate CDF with "upper" flag
        p = feval (cdfhandl{idx}, x, varargin{:});
      endif

    elseif (nargs == cdf_args{idx})
      ## Check that all distribution parameters are numeric
      if (! all (cellfun (@(x)isnumeric(x), (varargin))))
        error ("cdf: distribution parameters must be numeric.");
      endif
      ## Call appropriate CDF without "upper" flag
      p = feval (cdfhandl{idx}, x, varargin{:});

    else
      if (cdf_args{idx} == 1)
        error ("cdf: %s distribution requires 1 parameter.", name);
      else
        error ("cdf: %s distribution requires %d parameters.", ...
               name, cdf_args{idx});
      endif

    endif

  else
    error ("cdf: %s distribution is not implemented in Statistics.", name);
  endif

endfunction

## Test results
%!shared x
%! x = [1:5];
%!assert (cdf ("Beta", x, 5, 2), betacdf (x, 5, 2))
%!assert (cdf ("beta", x, 5, 2, "upper"), betacdf (x, 5, 2, "upper"))
%!assert (cdf ("Binomial", x, 5, 2), binocdf (x, 5, 2))
%!assert (cdf ("bino", x, 5, 2, "upper"), binocdf (x, 5, 2, "upper"))
%!assert (cdf ("Birnbaum-Saunders", x, 5, 2), bisacdf (x, 5, 2))
%!assert (cdf ("bisa", x, 5, 2, "upper"), bisacdf (x, 5, 2, "upper"))
%!assert (cdf ("Burr", x, 5, 2, 2), burrcdf (x, 5, 2, 2))
%!assert (cdf ("burr", x, 5, 2, 2, "upper"), burrcdf (x, 5, 2, 2, "upper"))
%!assert (cdf ("Cauchy", x, 5, 2), cauchycdf (x, 5, 2))
%!assert (cdf ("cauchy", x, 5, 2, "upper"), cauchycdf (x, 5, 2, "upper"))
%!assert (cdf ("Chi-squared", x, 5), chi2cdf (x, 5))
%!assert (cdf ("chi2", x, 5, "upper"), chi2cdf (x, 5, "upper"))
%!assert (cdf ("Extreme Value", x, 5, 2), evcdf (x, 5, 2))
%!assert (cdf ("ev", x, 5, 2, "upper"), evcdf (x, 5, 2, "upper"))
%!assert (cdf ("Exponential", x, 5), expcdf (x, 5))
%!assert (cdf ("exp", x, 5, "upper"), expcdf (x, 5, "upper"))
%!assert (cdf ("F-Distribution", x, 5, 2), fcdf (x, 5, 2))
%!assert (cdf ("f", x, 5, 2, "upper"), fcdf (x, 5, 2, "upper"))
%!assert (cdf ("Gamma", x, 5, 2), gamcdf (x, 5, 2))
%!assert (cdf ("gam", x, 5, 2, "upper"), gamcdf (x, 5, 2, "upper"))
%!assert (cdf ("Geometric", x, 5), geocdf (x, 5))
%!assert (cdf ("geo", x, 5, "upper"), geocdf (x, 5, "upper"))
%!assert (cdf ("Generalized Extreme Value", x, 5, 2, 2), gevcdf (x, 5, 2, 2))
%!assert (cdf ("gev", x, 5, 2, 2, "upper"), gevcdf (x, 5, 2, 2, "upper"))
%!assert (cdf ("Generalized Pareto", x, 5, 2, 2), gpcdf (x, 5, 2, 2))
%!assert (cdf ("gp", x, 5, 2, 2, "upper"), gpcdf (x, 5, 2, 2, "upper"))
%!assert (cdf ("Gumbel", x, 5, 2), gumbelcdf (x, 5, 2))
%!assert (cdf ("gumbel", x, 5, 2, "upper"), gumbelcdf (x, 5, 2, "upper"))
%!assert (cdf ("Half-normal", x, 5, 2), hncdf (x, 5, 2))
%!assert (cdf ("hn", x, 5, 2, "upper"), hncdf (x, 5, 2, "upper"))
%!assert (cdf ("Hypergeometric", x, 5, 2, 2), hygecdf (x, 5, 2, 2))
%!assert (cdf ("hyge", x, 5, 2, 2, "upper"), hygecdf (x, 5, 2, 2, "upper"))
%!assert (cdf ("Inverse Gaussian", x, 5, 2), invgcdf (x, 5, 2))
%!assert (cdf ("invg", x, 5, 2, "upper"), invgcdf (x, 5, 2, "upper"))
%!assert (cdf ("Laplace", x, 5, 2), laplacecdf (x, 5, 2))
%!assert (cdf ("laplace", x, 5, 2, "upper"), laplacecdf (x, 5, 2, "upper"))
%!assert (cdf ("Logistic", x, 5, 2), logicdf (x, 5, 2))
%!assert (cdf ("logi", x, 5, 2, "upper"), logicdf (x, 5, 2, "upper"))
%!assert (cdf ("Log-Logistic", x, 5, 2), loglcdf (x, 5, 2))
%!assert (cdf ("logl", x, 5, 2, "upper"), loglcdf (x, 5, 2, "upper"))
%!assert (cdf ("Lognormal", x, 5, 2), logncdf (x, 5, 2))
%!assert (cdf ("logn", x, 5, 2, "upper"), logncdf (x, 5, 2, "upper"))
%!assert (cdf ("Nakagami", x, 5, 2), nakacdf (x, 5, 2))
%!assert (cdf ("naka", x, 5, 2, "upper"), nakacdf (x, 5, 2, "upper"))
%!assert (cdf ("Negative Binomial", x, 5, 2), nbincdf (x, 5, 2))
%!assert (cdf ("nbin", x, 5, 2, "upper"), nbincdf (x, 5, 2, "upper"))
%!assert (cdf ("Noncentral F-Distribution", x, 5, 2, 2), ncfcdf (x, 5, 2, 2))
%!assert (cdf ("ncf", x, 5, 2, 2, "upper"), ncfcdf (x, 5, 2, 2, "upper"))
%!assert (cdf ("Noncentral Student T", x, 5, 2), nctcdf (x, 5, 2))
%!assert (cdf ("nct", x, 5, 2, "upper"), nctcdf (x, 5, 2, "upper"))
%!assert (cdf ("Noncentral Chi-Squared", x, 5, 2), ncx2cdf (x, 5, 2))
%!assert (cdf ("ncx2", x, 5, 2, "upper"), ncx2cdf (x, 5, 2, "upper"))
%!assert (cdf ("Normal", x, 5, 2), normcdf (x, 5, 2))
%!assert (cdf ("norm", x, 5, 2, "upper"), normcdf (x, 5, 2, "upper"))
%!assert (cdf ("Poisson", x, 5), poisscdf (x, 5))
%!assert (cdf ("poiss", x, 5, "upper"), poisscdf (x, 5, "upper"))
%!assert (cdf ("Rayleigh", x, 5), raylcdf (x, 5))
%!assert (cdf ("rayl", x, 5, "upper"), raylcdf (x, 5, "upper"))
%!assert (cdf ("Student T", x, 5), tcdf (x, 5))
%!assert (cdf ("t", x, 5, "upper"), tcdf (x, 5, "upper"))
%!assert (cdf ("Triangular", x, 5, 2, 2), tricdf (x, 5, 2, 2))
%!assert (cdf ("tri", x, 5, 2, 2, "upper"), tricdf (x, 5, 2, 2, "upper"))
%!assert (cdf ("Discrete Uniform", x, 5), unidcdf (x, 5))
%!assert (cdf ("unid", x, 5, "upper"), unidcdf (x, 5, "upper"))
%!assert (cdf ("Uniform", x, 5, 2), unifcdf (x, 5, 2))
%!assert (cdf ("unif", x, 5, 2, "upper"), unifcdf (x, 5, 2, "upper"))
%!assert (cdf ("Von Mises", x, 5, 2), vmcdf (x, 5, 2))
%!assert (cdf ("vm", x, 5, 2, "upper"), vmcdf (x, 5, 2, "upper"))
%!assert (cdf ("Weibull", x, 5, 2), wblcdf (x, 5, 2))
%!assert (cdf ("wbl", x, 5, 2, "upper"), wblcdf (x, 5, 2, "upper"))

## Test input validation
%!error<cdf: distribution NAME must a char string.> cdf (1)
%!error<cdf: distribution NAME must a char string.> cdf ({"beta"})
%!error<cdf: X must be numeric.> cdf ("beta", {[1 2 3 4 5]})
%!error<cdf: X must be numeric.> cdf ("beta", "text")
%!error<cdf: values in X must be real.> cdf ("beta", 1+i)
%!error<cdf: distribution parameters must be numeric.> ...
%! cdf ("Beta", x, "a", 2)
%!error<cdf: distribution parameters must be numeric.> ...
%! cdf ("Beta", x, 5, "")
%!error<cdf: distribution parameters must be numeric.> ...
%! cdf ("Beta", x, 5, {2})
%!error<cdf: chi2 distribution requires 1 parameter.> cdf ("chi2", x)
%!error<cdf: Beta distribution requires 2 parameters.> cdf ("Beta", x, 5)
%!error<cdf: Burr distribution requires 3 parameters.> cdf ("Burr", x, 5)
%!error<cdf: Burr distribution requires 3 parameters.> cdf ("Burr", x, 5, 2)
