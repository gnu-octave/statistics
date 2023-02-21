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
## @deftypefn  {statistics} {@var{retval} =} cdf (@var{name}, @var{x}, @dots{})
##
## Return the CDF of @var{name} distribution function for value @var{x}.
##
## This is a wrapper around various @qcode{name}cdf and @qcode{name}_cdf
## functions. See the corresponding functions' help to learn the signification
## of the arguments after @var{x}.
##
## @var{name} must be a char string of the name or the abbreviation of the
## desired cumulative distribution function as listed in the followng table.
## The last column shows the maximum number of extra arguments that can be
## passed to the desired CDF.
##
## @multitable @columnfractions 0.45 0.2 0.35
## @headitem Distribution Name @tab Abbreviation @tab Max. Extra Arguments
## @item @qcode{"Birnbaum-Saunders"} @tab @qcode{"bbs"} @tab 3
## @item @qcode{"Beta"} @tab @qcode{"beta"} @tab 3
## (including @qcode{"upper"})
## @item @qcode{"Binomial"} @tab @qcode{"bino"} @tab 3
## (including @qcode{"upper"})
## @item @qcode{"Burr"} @tab @qcode{"burr"} @tab 3
## @item @qcode{"Bivariate Normal"} @tab @qcode{"bvn"} @tab 2
## @item @qcode{"Cauchy"} @tab @qcode{"cauchy"} @tab 2
## @item @qcode{"Chi-square"} @tab @qcode{"chi2"} @tab 2
## (including @qcode{"upper"})
## @item @qcode{"Copula Family"} @tab @qcode{"copula"} @tab 3
## @item @qcode{"Extreme Value"} @tab @qcode{"ev"} @tab 5
## (including @qcode{"upper"})
## @item @qcode{"Exponential"} @tab @qcode{"exp"} @tab 4
## (including @qcode{"upper"})
## @item @qcode{"F-Distribution"} @tab @qcode{"f"} @tab 3
## (including @qcode{"upper"})
## @item @qcode{"Gamma"} @tab @qcode{"gam"} @tab 5
## (including @qcode{"upper"})
## @item @qcode{"Geometric"} @tab @qcode{"geo"} @tab 2
## (including @qcode{"upper"})
## @item @qcode{"Generalized Extreme Value"} @tab @qcode{"gev"} @tab 4
## (including @qcode{"upper"})
## @item @qcode{"Generalized Pareto"} @tab @qcode{"gp"} @tab 4
## (including @qcode{"upper"})
## @item @qcode{"Hypergeometric"} @tab @qcode{"hyge"} @tab 4
## (including @qcode{"upper"})
## @item @qcode{"Johnson SU"} @tab @qcode{"jsu"} @tab 2
## @item @qcode{"Laplace"} @tab @qcode{"laplace"} @tab 2
## @item @qcode{"Logistic"} @tab @qcode{"logistic"} @tab 2
## @item @qcode{"Lognormal"} @tab @qcode{"logn"} @tab 5
## (including @qcode{"upper"})
## @item @qcode{"Multivariate Normal"} @tab @qcode{"mvn"} @tab 4
## @item @qcode{"Multivariate Student T"} @tab @qcode{"mvt"} @tab 3
## @item @qcode{"Quasi-Monte-Carlo Multivariate Student T"}
## @tab @qcode{"mvtqmc"} @tab 6
## @item @qcode{"Nakagami"} @tab @qcode{"naka"} @tab 2
## @item @qcode{"Negative Binomial"} @tab @qcode{"nbin"} @tab 3
## (including @qcode{"upper"})
## @item @qcode{"Noncentral F-Distribution"} @tab @qcode{"ncf"} @tab 4
## (including @qcode{"upper"})
## @item  @qcode{"Noncentral Student T"} @tab @qcode{"nct"} @tab 3
## (including @qcode{"upper"})
## @item @qcode{"Noncentral Chi-Square"} @tab @qcode{"ncx2"} @tab 3
## (including @qcode{"upper"})
## @item @qcode{"Normal"} @tab @qcode{"norm"} @tab 5
## (including @qcode{"upper"})
## @item @qcode{"Poisson"} @tab @qcode{"poiss"} @tab 2
## (including @qcode{"upper"})
## @item @qcode{"Rayleigh"} @tab @qcode{"rayl"} @tab 2
## (including @qcode{"upper"})
## @item @qcode{"Standard Normal"} @tab @qcode{"stdnormal"} @tab 0
## @item @qcode{"Student T"} @tab @qcode{"t"} @tab 2
## (including @qcode{"upper"})
## @item @qcode{"Triangular"} @tab @qcode{"tri"} @tab 3
## @item @qcode{"Discrete Uniform"} @tab @qcode{"unid"} @tab 2
## (including @qcode{"upper"})
## @item @qcode{"Uniform"} @tab @qcode{"unif"} @tab 3
## (including @qcode{"upper"})
## @item @qcode{"Von Mises"} @tab @qcode{"vm"} @tab 2
## @item @qcode{"Weibull"} @tab @qcode{"wbl"} @tab 5
## (including @qcode{"upper"})
## @end multitable
##
## @seealso{icdf, pdf, random, bbscdf, betacdf, binocdf, burrcdf, bvncdf,
## cauchy_cdf, chi2cdf, copulacdf, evcdf, expcdf, fcdf, gamcdf, geocdf, gevcdf,
## gpcdf, hygecdf, jsucdf, laplace_cdf, logistic_cdf, logncdf, mvncdf, mvtcdf,
## mvtcdfqmc, nakacdf, nbincdf, ncfcdf, nctcdf, ncx2cdf, normcdf, poisscdf,
## raylcdf, stdnormal_cdf, tcdf, tricdf, unidcdf, unifcdf, vmcdf, wblcdf}
## @end deftypefn

function [retval] = cdf (name, varargin)
  ## implemented functions
  persistent allcdf = { ...
    {"bbs"      , "Birnbaum-Saunders"},         @bbscdf,       3, ...
    {"beta"     , "Beta"},                      @betacdf,      3, ... # "upper"
    {"bino"     , "Binomial"},                  @binocdf,      3, ... # "upper"
    {"burr"     , "Burr"},                      @burrcdf,      3, ...
    {"bvn"      , "Bivariate Normal"},          @bvncdf,       2, ...
    {"cauchy"   , "Cauchy"},                    @cauchy_cdf,   2, ...
    {"chi2"     , "Chi-square"},                @chi2cdf,      2, ... # "upper"
    {"copula"   , "Copula Family"},             @copulacdf,    3, ...
    {"ev"       , "Extreme Value"},             @evcdf,        5, ... # "upper"
    {"exp"      , "Exponential"},               @expcdf,       4, ... # "upper"
    {"f"        , "F-Distribution"},            @fcdf,         3, ... # "upper"
    {"gam"      , "Gamma"},                     @gamcdf,       5, ... # "upper"
    {"geo"      , "Geometric"},                 @geocdf,       2, ... # "upper"
    {"gev"      , "Generalized Extreme Value"}, @gevcdf,       4, ... # "upper"
    {"gp"       , "Generalized Pareto"},        @gpcdf,        4, ... # "upper"
    {"hyge"     , "Hypergeometric"},            @hygecdf,      4, ... # "upper"
    {"jsu"      , "Johnson SU"},                @jsucdf,       2, ...
    {"laplace"  , "Laplace"},                   @laplace_cdf,  2, ...
    {"logistic" , "Logistic"},                  @logistic_cdf, 2, ...
    {"logn"     , "Lognormal"},                 @logncdf,      5, ... # "upper"
    {"mvn"      , "Multivariate Normal"},       @mvncdf,       4, ...
    {"mvt"      , "Multivariate Student T"},    @mvtcdf,       3, ...
    {"mvtqmc"   , "Quasi-Monte-Carlo Multivariate Student T"}, @mvtcdfqmc,6, ...
    {"naka"     , "Nakagami"},                  @nakacdf,      2, ...
    {"nbin"     , "Negative Binomial"},         @nbincdf,      3, ... # "upper"
    {"ncf"      , "Noncentral F-Distribution"}, @ncfcdf,       4, ... # "upper"
    {"nct"      , "Noncentral Student T"},      @nctcdf,       3, ... # "upper"
    {"ncx2"     , "Noncentral Chi-Square"},     @ncx2cdf,      3, ... # "upper"
    {"norm"     , "Normal"},                    @normcdf,      5, ... # "upper"
    {"poiss"    , "Poisson"},                   @poisscdf,     2, ... # "upper"
    {"rayl"     , "Rayleigh"},                  @raylcdf,      2, ... # "upper"
    {"stdnormal", "Standard Normal"},           @stdnormal_cdf,0, ...
    {"t"        , "Student T"},                 @tcdf,         2, ... # "upper"
    {"tri"      , "Triangular"},                @tricdf,       3, ...
    {"unid"     , "Discrete Uniform"},          @unidcdf,      2, ... # "upper"
    {"unif"     , "Uniform"},                   @unifcdf,      3, ... # "upper"
    {"vm"       , "Von Mises"},                 @vmcdf,        2, ...
    {"wbl"      , "Weibull"},                   @wblcdf,       5};    # "upper"

  if (numel (varargin) < 1 || ! ischar (name))
    print_usage ();
  endif

  x = varargin{1};
  varargin(1) = [];
  nargs = numel (varargin);

  cdfnames = allcdf(1:3:end);
  cdfhandl = allcdf(2:3:end);
  cdf_args = allcdf(3:3:end);

  idx = cellfun (@(x)any(strcmpi (name, x)), cdfnames);
  ## Add special list
  special = {"copula", "Copula family", "mvt", "multivariate Student"};

  if (any (idx))

    if (nargs > cdf_args{idx})
      if (pdf_args{idx} == 1)
        error ("cdf: %s takes only 1 extra argument.", name);
      else
        error ("cdf: %s takes up to %d extra arguments.", name, cdf_args{idx});
      endif
    endif

    if (! any (strcmpi (name, special)))
      retval = feval (cdfhandl{idx}, x, varargin{:});

    elseif (any (strcmpi (name, special(1:2))))
      retval = feval (cdfhandl{idx}, varargin{1}, x, varargin{2:end});

    elseif (any (strcmpi (name, special(3:4))))

      if (nargs == 2)
        retval = feval (cdfhandl{idx}, x, varargin{:});
      else
        retval = feval (cdfhandl{idx}, varargin{1}, x, varargin{2:end});
      endif

    endif

  else
    error ("cdf: %s distribution is not implemented in Statistics.", name);
  endif

endfunction

%!test
%! assert (cdf ("norm", 1, 0, 1), normcdf (1, 0, 1))
%!test
%! x = [0.2:0.2:0.6; 0.2:0.2:0.6];
%! theta = [1; 2];
%! assert (cdf ("copula", x, "Clayton", theta), copulacdf ("Clayton", x, theta))
%!test
%! x = [-1, 0, 1, 2, Inf];
%! assert (cdf ("bbs", x, ones (1,5), ones (1,5), zeros (1,5)), ...
%!         bbscdf (x, ones (1,5), ones (1,5), zeros (1,5)))
%!test
%! x = [1 2];
%! mu = [0.5 1.5];
%! sigma = [1.0 0.5; 0.5 1.0];
%! assert (cdf ("multivariate normal", x, mu, sigma), ...
%!         mvncdf (x, mu, sigma), 0.01)
%!test
%! x = 1.2;
%! mu = [1.0 0.5; 0.5 1.0];
%! sigma = [3];
%! assert (cdf ("norm", x, mu, sigma), normcdf (x, mu, sigma), 0.01)

%!error pdf ("Birnbaum-Saunders", 1, 1, 2, 3, 4)
%!error pdf ("some", 1, 2)
