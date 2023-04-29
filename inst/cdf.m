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
## @seealso{icdf, pdf, random, betacdf, binocdf, bisacdf, burrcdf, cauchycdf,
## chi2cdf, evcdf, expcdf, fcdf, gamcdf, geocdf, gevcdf, gpcdf, hygecdf,
## laplacecdf, logicdf, loglcdf, logncdf, nakacdf, nbincdf, ncfcdf, nctcdf,
## ncx2cdf, normcdf, poisscdf, raylcdf, tcdf, tricdf, unidcdf, unifcdf, wblcdf}
## @end deftypefn

function p = cdf (name, varargin)

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
    {"hyge"     , "Hypergeometric"},            @hygecdf,      3, ...
    {"laplace"  , "Laplace"},                   @laplacecdf,   2, ...
    {"logi"     , "Logistic"},                  @logicdf,      2, ...
    {"logl"     , "Log-Logistic"},              @logicdf,      2, ...
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
    {"wbl"      , "Weibull"},                   @wblcdf,       2};

  if (numel (varargin) < 1 || ! ischar (name))
    print_usage ();
  endif

  ## Get data
  x = varargin{1};
  varargin(1) = [];

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
        error ("cdf: %s requires 1 parameter.", name);
      else
        error ("cdf: %s requires %d parameters.", name, cdf_args{idx});
      endif

    endif

  else
    error ("cdf: %s distribution is not implemented in Statistics.", name);
  endif

endfunction

## Test results
%!test
%! assert (cdf ("norm", 1, 0, 1), normcdf (1, 0, 1))
