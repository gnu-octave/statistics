## Copyright (C) 2013 Pantxo Diribarne
## Copyright (C) 2022 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## desired cumulative distribution functionas listed in the followng table.
## The last column shows the maximum number of extra arguments that can be
## passed to the desired CDF.
##
## @multitable @columnfractions 0.02 0.43 0.2 0.35
## @headitem @tab Distribution Name @tab Abbreviation @tab Max. Extra Arguments
## @item @tab @qcode{"Birnbaum-Saunders"} @tab @qcode{"bbs"} @tab 3
## @item @tab @qcode{"Beta"} @tab @qcode{"beta"} @tab 3
## (including @qcode{"upper"})
## @item @tab @qcode{"Binomial"} @tab @qcode{"bino"} @tab 3
## (including @qcode{"upper"})
## @item @tab @qcode{"Burr"} @tab @qcode{"burr"} @tab 2
## @item @tab @qcode{"Bivariate"} @tab @qcode{"bvn"} @tab 2
## @item @tab @qcode{"Cauchy"} @tab @qcode{"cauchy"} @tab 2
## @item @tab @qcode{"Chi-square"} @tab @qcode{"chi2"} @tab 2
## (including @qcode{"upper"})
## @item @tab @qcode{"Copula family"} @tab @qcode{"copula"} @tab 3
## @item @tab @qcode{"Extreme value"} @tab @qcode{"ev"} @tab 5
## (including @qcode{"upper"})
## @item @tab @qcode{"Exponential"} @tab @qcode{"exp"} @tab 4
## (including @qcode{"upper"})
## @item @tab @qcode{"F-Distribution"} @tab @qcode{"f"} @tab 3
## (including @qcode{"upper"})
## @item @tab @qcode{"Gamma"} @tab @qcode{"gam"} @tab 5
## (including @qcode{"upper"})
## @item @tab @qcode{"Geometric"} @tab @qcode{"geo"} @tab 2
## (including @qcode{"upper"})
## @item @tab @qcode{"Generalized extreme value"} @tab @qcode{"gev"} @tab 4
## (including @qcode{"upper"})
## @item @tab @qcode{"Generalized Pareto"} @tab @qcode{"gp"} @tab 4
## (including @qcode{"upper"})
## @item @tab @qcode{"Hypergeometric"} @tab @qcode{"hyge"} @tab 4
## (including @qcode{"upper"})
## @item @tab @qcode{"Johnson SU"} @tab @qcode{"jsu"} @tab 2
## @item @tab @qcode{"Laplace"} @tab @qcode{"laplace"} @tab 2
## @item @tab @qcode{"Logistic"} @tab @qcode{"logistic"} @tab 2
## @item @tab @qcode{"Lognormal"} @tab @qcode{"logn"} @tab 5
## (including @qcode{"upper"})
## @item @tab @qcode{"Multivariate normal"} @tab @qcode{"mvn"} @tab 4
## @item @tab @qcode{"Multivariate Student T"} @tab @qcode{"mvt"} @tab 3
## @item @tab @qcode{"Quasi-Monte-Carlo Multivariate Student T"}
## @tab @qcode{"mvtqmc"} @tab 6
## @item @tab @qcode{"Nakagami"} @qcode{"naka"} @tab @tab 2
## @item @tab @qcode{"Negative binomial"} @tab @qcode{"nbin"} @tab 3
## (including @qcode{"upper"})
## @item @tab @qcode{"Noncentral F-Distribution"} @tab @qcode{"ncf"} @tab 4
## (including @qcode{"upper"})
## @item @tab @qcode{"Noncentral Student T"} @tab @qcode{"nct"} @tab 3
## (including @qcode{"upper"})
## @item @tab @qcode{"Noncentral Chi-Square"} @tab @qcode{"ncx2"} @tab 3
## (including @qcode{"upper"})
## @item @tab @qcode{"Normal"} @qcode{"norm"} @tab @tab 5
## (including @qcode{"upper"})
## @item @tab @qcode{"Poisson"} @qcode{"poiss"} @tab @tab 2
## (including @qcode{"upper"})
## @item @tab @qcode{"Rayleigh"} @qcode{"rayl"} @tab @tab 2
## (including @qcode{"upper"})
## @item @tab @qcode{"Standard normal"} @tab @qcode{"stdnormal"} @tab 0
## @item @tab @qcode{"Student T"} @tab @qcode{"t"} @tab 2
## (including @qcode{"upper"})
## @item @tab @qcode{"Triangular"} @tab @qcode{"tri"} @tab 3
## @item @tab @qcode{"Discrete uniform"} @tab @qcode{"unid"} @tab 2
## (including @qcode{"upper"})
## @item @tab @qcode{"Uniform"} @tab @qcode{"unif"} @tab 3
## (including @qcode{"upper"})
## @item @tab @qcode{"Von Mises"} @tab @qcode{"vm"} @tab 2
## @item @tab @qcode{"Weibull"} @tab @qcode{"wbl"} @tab 5
## (including @qcode{"upper"})
## @end multitable
##
## @seealso{pdf, rnd}
## @end deftypefn

function [retval] = cdf (name, varargin)
  ## implemented functions
  persistent allcdf = { ...
    {"bbs"      , "Birnbaum-Saunders"},         @bbscdf,       3, ...
    {"beta"     , "Beta"},                      @betacdf,      3, ... # "upper"
    {"bino"     , "Binomial"},                  @binocdf,      3, ... # "upper"
    {"burr"     , "Burr"},                      @burrcdf,      2, ...
    {"bvn"      , "Bivariate"},                 @bvncdf,       2, ...
    {"cauchy"   , "Cauchy"},                    @cauchy_cdf,   2, ...
    {"chi2"     , "chi-square"},                @chi2cdf,      2, ... # "upper"
    {"copula"   , "Copula family"},             @copulacdf,    3, ...
    {"ev"       , "Extreme value"},             @evcdf,        5, ... # "upper"
    {"exp"      , "Exponential"},               @expcdf,       4, ... # "upper"
    {"f"        , "F-Distribution"},            @fcdf,         3, ... # "upper"
    {"gam"      , "Gamma"},                     @gamcdf,       5, ... # "upper"
    {"geo"      , "Geometric"},                 @geocdf,       2, ... # "upper"
    {"gev"      , "Generalized extreme value"}, @gevcdf,       4, ... # "upper"
    {"gp"       , "Generalized Pareto"},        @gpcdf,        4, ... # "upper"
    {"hyge"     , "Hypergeometric"},            @hygecdf,      4, ... # "upper"
    {"jsu"      , "Johnson SU"},                @jsucdf,       2, ...
    {"laplace"  , "Laplace"},                   @laplace_cdf,  2, ...
    {"logistic" , "Logistic"},                  @logistic_cdf, 2, ...
    {"logn"     , "Lognormal"},                 @logncdf,      5, ... # "upper"
    {"mvn"      , "Multivariate normal"},       @mvncdf,       4, ...
    {"mvt"      , "Multivariate Student T"},    @mvtcdf,       3, ...
    {"mvtqmc"   , "Quasi-Monte-Carlo Multivariate Student T"}, @mvtcdfqmc,6, ...
    {"naka"     , "Nakagami"},                  @nakacdf,      2, ...
    {"nbin"     , "Negative binomial"},         @nbincdf,      3, ... # "upper"
    {"ncf"      , "Noncentral F-Distribution"}, @ncfcdf,       4, ... # "upper"
    {"nct"      , "Noncentral Student T"},      @nctcdf,       3, ... # "upper"
    {"ncx2"     , "Noncentral Chi-Square"},     @ncx2cdf,      3, ... # "upper"
    {"norm"     , "Normal"},                    @normcdf,      5, ... # "upper"
    {"poiss"    , "Poisson"},                   @poisscdf,     2, ... # "upper"
    {"rayl"     , "Rayleigh"},                  @raylcdf,      2, ... # "upper"
    {"stdnormal", "Standard normal"},           @stdnormal_cdf,0, ...
    {"t"        , "Student T"},                 @tcdf,         2, ... # "upper"
    {"tri"      , "Triangular"},                @tricdf,       3, ...
    {"unid"     , "Discrete uniform"},          @unidcdf,      2, ... # "upper"
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
  cdfargs  = allcdf(3:3:end);

  idx = cellfun (@(x)any(strcmpi (name, x)), cdfnames);
  ## Add special list
  special = {"copula", "Copula family", "mvt", "multivariate Student"};

  if (any (idx))

    if (nargs > cdfargs{idx})
      error ("cdf: %s takes up to %d extra arguments.", name, cdfargs{idx})
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
