## Copyright (C) 2007 Soren Hauberg <soren@hauberg.org>
## Copyright (C) 2023 Andreas Bertsatos <abertsatos@biol.uoa.gr>
##
## This file is part of the statistics package for GNU Octave.
##
## This program is free software; you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation; either version 3 of the License, or (at your option) any later
## version.
##
## This program is distributed in the hope that it will be useful, but WITHOUT
## ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
## FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
## details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn  {statistics} {@var{r} =} random(@var{name}, @var{arg1})
## @deftypefnx {statistics} {@var{r} =} random(@var{name}, @var{arg1}, @var{arg2})
## @deftypefnx {statistics} {@var{r} =} random(@var{name}, @var{arg1}, @var{arg2}, @var{arg3})
## @deftypefnx {statistics} {@var{r} =} random(@var{name}, @dots{}, @var{rows}, @var{cols})
## @deftypefnx {statistics} {@var{r} =} random(@var{name}, @dots{}, @var{rows}, @var{cols}, @dots{})
## @deftypefnx {statistics} {@var{r} =} random(@var{name}, @dots{}, [@var{sz}])
##
## Random arrays from from a given one-, two-, or three-parameter distribution.
##
## The variable @var{name} must be a string that names the distribution from
## which to sample.  If this distribution is a one-parameter distribution
## @var{arg1} should be supplied, if it is a two-paramter distribution
## @var{arg2} must also be supplied, and if it is a three-parameter distribution
## @var{arg3} must also be present.  Any arguments following the distribution
## parameters will determine the size of the result.
##
## When called with a single size argument, return a square matrix with
## the dimension specified.  When called with more than one scalar argument the
## first two arguments are taken as the number of rows and columns and any
## further arguments specify additional matrix dimensions.  The size may also
## be specified with a vector of dimensions @var{sz}.
##
## As an example, the following code generates a 10 by 20 matrix containing
## random numbers from a normal distribution with mean 5 and standard deviation
## 2.
## @example
## R = random("normal", 5, 2, [10, 20]);
## @end example
##
## @var{name} must be a char string of the name or the abbreviation of the
## desired probability distribution function as listed in the followng table.
## The last column shows the required number of parameters that must be passed
## passed to the desired @qcode{*rnd} distribution function.
##
## @multitable @columnfractions 0.4 0.05 0.2 0.05 0.3
## @headitem Distribution Name @tab @tab Abbreviation @tab @tab
## Required Parameters
## @item @qcode{"Birnbaum-Saunders"} @tab @tab @qcode{"bbs"} @tab @tab 3
## @item @qcode{"Beta"} @tab @tab @qcode{"beta"} @tab @tab 2
## @item @qcode{"Binomial"} @tab @tab @qcode{"bino"} @tab @tab 2
## @item @qcode{"Burr"} @tab @tab @qcode{"burr"} @tab @tab 3
## @item @qcode{"Cauchy"} @tab @tab @qcode{"cauchy"} @tab @tab 2
## @item @qcode{"Chi-square"} @tab @tab @qcode{"chi2"} @tab @tab 1
## @item @qcode{"Copula Family"} @tab @tab @qcode{"copula"} @tab @tab
## 2-4 no size args
## @item @qcode{"Extreme Value"} @tab @tab @qcode{"ev"} @tab @tab 2
## @item @qcode{"Exponential"} @tab @tab @qcode{"exp"} @tab @tab 1
## @item @qcode{"F-Distribution"} @tab @tab @qcode{"f"} @tab @tab 2
## @item @qcode{"Gamma"} @tab @tab @qcode{"gam"} @tab @tab 2
## @item @qcode{"Geometric"} @tab @tab @qcode{"geo"} @tab @tab 1
## @item @qcode{"Generalized Extreme Value"} @tab @tab @qcode{"gev"} @tab @tab 3
## @item @qcode{"Generalized Pareto"} @tab @tab @qcode{"gp"} @tab @tab 3
## @item @qcode{"Hypergeometric"} @tab @tab @qcode{"hyge"} @tab @tab 3
## @item @qcode{"Inverse Wishart"} @tab @tab @qcode{"iwish"} @tab @tab
## 3-4 no size args
## @item @qcode{"Laplace"} @tab @tab @qcode{"laplace"} @tab @tab 2
## @item @qcode{"Logistic"} @tab @tab @qcode{"logistic"} @tab @tab 2
## @item @qcode{"Lognormal"} @tab @tab @qcode{"logn"} @tab @tab 2
## @item @qcode{"Multinomial"} @tab @tab @qcode{"mn"} @tab @tab 2
## @item @qcode{"Multivariate Normal"} @tab @tab @qcode{"mvn"} @tab @tab
## 2-4 no size args
## @item @qcode{"Multivariate Student T"} @tab @tab @qcode{"mvt"} @tab @tab 2
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
## @item @qcode{"Von Mises"} @tab @tab @qcode{"vm"} @tab @tab 2
## @item @qcode{"Weibull"} @tab @tab @qcode{"wbl"} @tab @tab 2
## @item @qcode{"Wiener Process"} @tab @tab @qcode{"wien"} @tab @tab
## 1-3 no size args
## @item @qcode{"Wishart"} @tab @tab @qcode{"wish"} @tab @tab
## 3-4 no size args
## @end multitable
##
## @seealso{cdf, icdf, pdf, bbsrnd, betarnd, binornd, burrrnd, cauchy_rnd,
## chi2rnd, copularnd, evrnd, exprnd, frnd, gamrnd, geornd, gevrnd, gprnd,
## hygernd, iwishrnd, laplace_rnd, logistic_rnd, lognrnd, mnrnd, mvnrnd,
## mvtrnd, nakarnd, nbinrnd, ncfrnd, nctrnd, ncx2rnd, normrnd, poissrnd,
## raylrnd, stdnormal_rnd, trnd, trirnd, unidrnd, unifrnd, vmrnd, wblrnd,
## wienrnd, wishrnd}
## @end deftypefn

function retval = random (name, varargin)

  ## implemented functions
  persistent allpdf = { ...
    {"bbs"      , "Birnbaum-Saunders"},         @bbsrnd,             3, ...
    {"beta"     , "Beta"},                      @betarnd,            2, ...
    {"bino"     , "Binomial"},                  @binornd,            2, ...
    {"burr"     , "Burr"},                      @burrrnd,            3, ...
    {"cauchy"   , "Cauchy"},                    @cauchy_rnd,         2, ...
    {"chi2"     , "Chi-square"},                @chi2rnd,            1, ...
    {"copula"   , "Copula Family"},             @copularnd,    [2 3 4], ...
    {"ev"       , "Extreme Value"},             @evrnd,              2, ...
    {"exp"      , "Exponential"},               @exprnd,             1, ...
    {"f"        , "F-Distribution"},            @frnd,               2, ...
    {"gam"      , "Gamma"},                     @gamrnd,             2, ...
    {"geo"      , "Geometric"},                 @geornd,             1, ...
    {"gev"      , "Generalized Extreme Value"}, @gevrnd,             3, ...
    {"gp"       , "Generalized Pareto"},        @gprnd,              3, ...
    {"hyge"     , "Hypergeometric"},            @hygernd,            4, ...
    {"iwish"    , "Inverse Wishart"},           @iwishrnd,       [3 4], ...
    {"laplace"  , "Laplace"},                   @laplace_rnd,        2, ...
    {"logistic" , "Logistic"},                  @logistic_rnd,       2, ...
    {"logn"     , "Lognormal"},                 @lognrnd,            2, ...
    {"mn"       , "Multinomial"},               @mnrnd,              2, ...
    {"mvn"      , "Multivariate Normal"},       @mvnrnd,       [2 3 4], ...
    {"mvt"      , "Multivariate Student T"},    @mvtrnd,             2, ...
    {"naka"     , "Nakagami"},                  @nakarnd,            2, ...
    {"nbin"     , "Negative Binomial"},         @nbinrnd,            2, ...
    {"ncf"      , "Noncentral F-Distribution"}, @ncfrnd,             3, ...
    {"nct"      , "Noncentral Student T"},      @nctrnd,             2, ...
    {"ncx2"     , "Noncentral Chi-Square"},     @ncx2rnd,            2, ...
    {"norm"     , "Normal"},                    @normrnd,            2, ...
    {"poiss"    , "Poisson"},                   @poissrnd,           1, ...
    {"rayl"     , "Rayleigh"},                  @raylrnd,            1, ...
    {"stdnormal", "Standard Normal"},           @stdnormal_rnd,      0, ...
    {"t"        , "Student T"},                 @trnd,               1, ...
    {"tri"      , "Triangular"},                @trirnd,             3, ...
    {"unid"     , "Discrete Uniform"},          @unidrnd,            1, ...
    {"unif"     , "Uniform"},                   @unifrnd,            2, ...
    {"vm"       , "Von Mises"},                 @vmrnd,              2, ...
    {"wbl"      , "Weibull"},                   @wblrnd,             2, ...
    {"wien"     , "Wiener Process"},            @wienrnd,      [1 2 3], ...
    {"wish"     , "Wishart"},                   @wishrnd,        [3 4]};

  if (numel (varargin) < 1 || ! ischar (name))
    print_usage ();
  endif

  nargs = numel (varargin);

  rndnames = allpdf(1:3:end);
  rndhandl = allpdf(2:3:end);
  rnd_args = allpdf(3:3:end);

  idx = cellfun (@(x)any(strcmpi (name, x)), rndnames);

  if (any (idx))

    ## Check rnd functions with no size arguments
    if (numel (rnd_args{idx}) > 1)

      if (any (nargs == rnd_args{idx}))
        retval = feval (rndhandl{idx}, varargin{:});
      else
        if (numel (rndargs{idx}) == 2)
          error ("random: %s requires %d or %d parameters.", ...
                 name, rnd_args{idx});
        else
          error ("random: %s requires %d, %d or %d parameters.", ...
                 name, rnd_args{idx});
        endif
      endif

    else
      retval = feval (rndhandl{idx}, varargin{:});
    endif

  else
    error ("pdf: %s distribution is not implemented in Statistics.", name);
  endif

endfunction

%!assert (size (random ("Birnbaum-Saunders", 5, 2, 2, 10)), size (bbsrnd (5, 2, 2, 10)))
%!assert (size (random ("normal", 5, 2, [10, 20])), size (normrnd (5, 2, 10, 20)))
%!assert (size (random ("beta", 5, 2, [10, 20])), size (betarnd (5, 2, 10, 20)))

%!error random ("copula", 1)
%!error random ("copula", 1, 2, 3, 4, 5)
%!error random ("iwish", 1)
%!error random ("iwish", 1, 2)
%!error random ("iwish", 1, 2, 3, 4, 5)
%!error random ("mvn", 1)
%!error random ("mvn", 1, 2, 3, 4, 5)
%!error random ("wien", 1, 2, 3, 4)
%!error random ("wish", 1)
%!error random ("wish", 1, 2)
%!error random ("wish", 1, 2, 3, 4, 5)
%!error random ("some", 1, 2, 3, 4, 5)
