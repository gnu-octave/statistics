## Copyright (C) 2007 Soren Hauberg <soren@hauberg.org>
## Copyright (C) 2023-2024 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {statistics} {@var{r} =} random (@var{name}, @var{A})
## @deftypefnx {statistics} {@var{r} =} random (@var{name}, @var{A}, @var{B})
## @deftypefnx {statistics} {@var{r} =} random (@var{name}, @var{A}, @var{B}, @var{C})
## @deftypefnx {statistics} {@var{r} =} random (@var{name}, @dots{}, @var{rows}, @var{cols})
## @deftypefnx {statistics} {@var{r} =} random (@var{name}, @dots{}, @var{rows}, @var{cols}, @dots{})
## @deftypefnx {statistics} {@var{r} =} random (@var{name}, @dots{}, [@var{sz}])
##
## Random arrays from a given one-, two-, or three-parameter distribution.
##
## The variable @var{name} must be a string with the name of the distribution to
## sample from.  If this distribution is a one-parameter distribution, @var{A}
## must be supplied, if it is a two-parameter distribution, @var{B} must also be
## supplied, and if it is a three-parameter distribution, @var{C} must also be
## supplied.  Any arguments following the distribution parameters will determine
## the size of the result.
##
## When called with a single size argument, return a square matrix with the
## dimension specified.  When called with more than one scalar argument the
## first two arguments are taken as the number of rows and columns and any
## further arguments specify additional matrix dimensions.  The size may also
## be specified with a vector of dimensions @var{sz}.
##
## @var{name} must be a char string of the name or the abbreviation of the
## desired probability distribution function as listed in the followng table.
## The last column shows the required number of parameters that must be passed
## passed to the desired @qcode{*rnd} distribution function.
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
## @item @qcode{"Rician"} @tab @tab @qcode{"rice"} @tab @tab 2
## @item @qcode{"Student T"} @tab @tab @qcode{"t"} @tab @tab 1
## @item @qcode{"Triangular"} @tab @tab @qcode{"tri"} @tab @tab 3
## @item @qcode{"Discrete Uniform"} @tab @tab @qcode{"unid"} @tab @tab 1
## @item @qcode{"Uniform"} @tab @tab @qcode{"unif"} @tab @tab 2
## @item @qcode{"Von Mises"} @tab @tab @qcode{"vm"} @tab @tab 2
## @item @qcode{"Weibull"} @tab @tab @qcode{"wbl"} @tab @tab 2
## @end multitable
##
## @seealso{cdf, icdf, pdf, betarnd, binornd, bisarnd, burrrnd, cauchyrnd,
## chi2rnd, evrnd, exprnd, frnd, gamrnd, geornd, gevrnd, gprnd, gumbelrnd,
## hnrnd, hygernd, invgrnd, laplacernd, logirnd, loglrnd, lognrnd, nakarnd,
## nbinrnd, ncfrnd, nctrnd, ncx2rnd, normrnd, poissrnd, raylrnd, ricernd, trnd,
## trirnd, unidrnd, unifrnd, vmrnd, wblrnd}
## @end deftypefn

function r = random (name, varargin)

  ## implemented functions
  persistent allDF = { ...
    {"beta"     , "Beta"},                      @betarnd,      2, ...
    {"bino"     , "Binomial"},                  @binornd,      2, ...
    {"bisa"     , "Birnbaum-Saunders"},         @bisarnd,      2, ...
    {"burr"     , "Burr"},                      @burrrnd,      3, ...
    {"cauchy"   , "Cauchy"},                    @cauchyrnd,    2, ...
    {"chi2"     , "Chi-squared"},               @chi2rnd,      1, ...
    {"ev"       , "Extreme Value"},             @evrnd,        2, ...
    {"exp"      , "Exponential"},               @exprnd,       1, ...
    {"f"        , "F-Distribution"},            @frnd,         2, ...
    {"gam"      , "Gamma"},                     @gamrnd,       2, ...
    {"geo"      , "Geometric"},                 @geornd,       1, ...
    {"gev"      , "Generalized Extreme Value"}, @gevrnd,       3, ...
    {"gp"       , "Generalized Pareto"},        @gprnd,        3, ...
    {"gumbel"   , "Gumbel"},                    @gumbelrnd,    2, ...
    {"hn"       , "Half-normal"},               @hnrnd,        2, ...
    {"hyge"     , "Hypergeometric"},            @hygernd,      3, ...
    {"invg"     , "Inverse Gaussian"},          @invgrnd,      2, ...
    {"laplace"  , "Laplace"},                   @laplacernd,   2, ...
    {"logi"     , "Logistic"},                  @logirnd,      2, ...
    {"logl"     , "Log-Logistic"},              @loglrnd,      2, ...
    {"logn"     , "Lognormal"},                 @lognrnd,      2, ...
    {"naka"     , "Nakagami"},                  @nakarnd,      2, ...
    {"nbin"     , "Negative Binomial"},         @nbinrnd,      2, ...
    {"ncf"      , "Noncentral F-Distribution"}, @ncfrnd,       3, ...
    {"nct"      , "Noncentral Student T"},      @nctrnd,       2, ...
    {"ncx2"     , "Noncentral Chi-squared"},    @ncx2rnd,      2, ...
    {"norm"     , "Normal"},                    @normrnd,      2, ...
    {"poiss"    , "Poisson"},                   @poissrnd,     1, ...
    {"rayl"     , "Rayleigh"},                  @raylrnd,      1, ...
    {"rice"     , "Rician"},                    @ricernd,      2, ...
    {"t"        , "Student T"},                 @trnd,         1, ...
    {"tri"      , "Triangular"},                @trirnd,       3, ...
    {"unid"     , "Discrete Uniform"},          @unidrnd,      1, ...
    {"unif"     , "Uniform"},                   @unifrnd,      2, ...
    {"vm"       , "Von Mises"},                 @vmrnd,        2, ...
    {"wbl"      , "Weibull"},                   @wblrnd,       2};

  if (! ischar (name))
    error ("random: distribution NAME must a char string.");
  endif

  ## Get number of arguments
  nargs = numel (varargin);

  ## Get available functions
  rndnames = allDF(1:3:end);
  rndhandl = allDF(2:3:end);
  rnd_args = allDF(3:3:end);

  ## Search for RND function
  idx = cellfun (@(x)any(strcmpi (name, x)), rndnames);

  if (any (idx))

    if (nargs == rnd_args{idx})
      ## Check that all distribution parameters are numeric
      if (! all (cellfun (@(x)isnumeric(x), (varargin))))
        error ("random: distribution parameters must be numeric.");
      endif
      ## Call appropriate RND
      r = feval (rndhandl{idx}, varargin{:});

    elseif (nargs > rnd_args{idx})
      ## Check that all distribution parameters are numeric
      if (! all (cellfun (@(x)isnumeric(x), (varargin(1:rnd_args{idx})))))
        error ("random: distribution parameters must be numeric.");
      endif
      ## Call appropriate RND. SIZE arguments are checked by the RND function.
      r = feval (rndhandl{idx}, varargin{:});

    else
      if (rnd_args{idx} == 1)
        error ("random: %s distribution requires 1 parameter.", name);
      else
        error ("random: %s distribution requires %d parameters.", ...
               name, rnd_args{idx});
      endif

    endif

  else
    error ("random: %s distribution is not implemented in Statistics.", name);
  endif

endfunction

## Test results
%!assert (size (random ("Beta", 5, 2, 2, 10)), size (betarnd (5, 2, 2, 10)))
%!assert (size (random ("beta", 5, 2, 2, 10)), size (betarnd (5, 2, 2, 10)))
%!assert (size (random ("Binomial", 5, 2, [10, 20])), size (binornd (5, 2, 10, 20)))
%!assert (size (random ("bino", 5, 2, [10, 20])), size (binornd (5, 2, 10, 20)))
%!assert (size (random ("Birnbaum-Saunders", 5, 2, [10, 20])), size (bisarnd (5, 2, 10, 20)))
%!assert (size (random ("bisa", 5, 2, [10, 20])), size (bisarnd (5, 2, 10, 20)))
%!assert (size (random ("Burr", 5, 2, 2, [10, 20])), size (burrrnd (5, 2, 2, 10, 20)))
%!assert (size (random ("burr", 5, 2, 2, [10, 20])), size (burrrnd (5, 2, 2, 10, 20)))
%!assert (size (random ("Cauchy", 5, 2, [10, 20])), size (cauchyrnd (5, 2, 10, 20)))
%!assert (size (random ("cauchy", 5, 2, [10, 20])), size (cauchyrnd (5, 2, 10, 20)))
%!assert (size (random ("Chi-squared", 5, [10, 20])), size (chi2rnd (5, 10, 20)))
%!assert (size (random ("chi2", 5, [10, 20])), size (chi2rnd (5, 10, 20)))
%!assert (size (random ("Extreme Value", 5, 2, [10, 20])), size (evrnd (5, 2, 10, 20)))
%!assert (size (random ("ev", 5, 2, [10, 20])), size (evrnd (5, 2, 10, 20)))
%!assert (size (random ("Exponential", 5, [10, 20])), size (exprnd (5, 10, 20)))
%!assert (size (random ("exp", 5, [10, 20])), size (exprnd (5, 10, 20)))
%!assert (size (random ("F-Distribution", 5, 2, [10, 20])), size (frnd (5, 2, 10, 20)))
%!assert (size (random ("f", 5, 2, [10, 20])), size (frnd (5, 2, 10, 20)))
%!assert (size (random ("Gamma", 5, 2, [10, 20])), size (gamrnd (5, 2, 10, 20)))
%!assert (size (random ("gam", 5, 2, [10, 20])), size (gamrnd (5, 2, 10, 20)))
%!assert (size (random ("Geometric", 5, [10, 20])), size (geornd (5, 10, 20)))
%!assert (size (random ("geo", 5, [10, 20])), size (geornd (5, 10, 20)))
%!assert (size (random ("Generalized Extreme Value", 5, 2, 2, [10, 20])), size (gevrnd (5, 2, 2, 10, 20)))
%!assert (size (random ("gev", 5, 2, 2, [10, 20])), size (gevrnd (5, 2, 2, 10, 20)))
%!assert (size (random ("Generalized Pareto", 5, 2, 2, [10, 20])), size (gprnd (5, 2, 2, 10, 20)))
%!assert (size (random ("gp", 5, 2, 2, [10, 20])), size (gprnd (5, 2, 2, 10, 20)))
%!assert (size (random ("Gumbel", 5, 2, [10, 20])), size (gumbelrnd (5, 2, 10, 20)))
%!assert (size (random ("gumbel", 5, 2, [10, 20])), size (gumbelrnd (5, 2, 10, 20)))
%!assert (size (random ("Half-normal", 5, 2, [10, 20])), size (hnrnd (5, 2, 10, 20)))
%!assert (size (random ("hn", 5, 2, [10, 20])), size (hnrnd (5, 2, 10, 20)))
%!assert (size (random ("Hypergeometric", 5, 2, 2, [10, 20])), size (hygernd (5, 2, 2, 10, 20)))
%!assert (size (random ("hyge", 5, 2, 2, [10, 20])), size (hygernd (5, 2, 2, 10, 20)))
%!assert (size (random ("Inverse Gaussian", 5, 2, [10, 20])), size (invgrnd (5, 2, 10, 20)))
%!assert (size (random ("invg", 5, 2, [10, 20])), size (invgrnd (5, 2, 10, 20)))
%!assert (size (random ("Laplace", 5, 2, [10, 20])), size (laplacernd (5, 2, 10, 20)))
%!assert (size (random ("laplace", 5, 2, [10, 20])), size (laplacernd (5, 2, 10, 20)))
%!assert (size (random ("Logistic", 5, 2, [10, 20])), size (logirnd (5, 2, 10, 20)))
%!assert (size (random ("logi", 5, 2, [10, 20])), size (logirnd (5, 2, 10, 20)))
%!assert (size (random ("Log-Logistic", 5, 2, [10, 20])), size (loglrnd (5, 2, 10, 20)))
%!assert (size (random ("logl", 5, 2, [10, 20])), size (loglrnd (5, 2, 10, 20)))
%!assert (size (random ("Lognormal", 5, 2, [10, 20])), size (lognrnd (5, 2, 10, 20)))
%!assert (size (random ("logn", 5, 2, [10, 20])), size (lognrnd (5, 2, 10, 20)))
%!assert (size (random ("Nakagami", 5, 2, [10, 20])), size (nakarnd (5, 2, 10, 20)))
%!assert (size (random ("naka", 5, 2, [10, 20])), size (nakarnd (5, 2, 10, 20)))
%!assert (size (random ("Negative Binomial", 5, 2, [10, 20])), size (nbinrnd (5, 2, 10, 20)))
%!assert (size (random ("nbin", 5, 2, [10, 20])), size (nbinrnd (5, 2, 10, 20)))
%!assert (size (random ("Noncentral F-Distribution", 5, 2, 2, [10, 20])), size (ncfrnd (5, 2, 2, 10, 20)))
%!assert (size (random ("ncf", 5, 2, 2, [10, 20])), size (ncfrnd (5, 2, 2, 10, 20)))
%!assert (size (random ("Noncentral Student T", 5, 2, [10, 20])), size (nctrnd (5, 2, 10, 20)))
%!assert (size (random ("nct", 5, 2, [10, 20])), size (nctrnd (5, 2, 10, 20)))
%!assert (size (random ("Noncentral Chi-Squared", 5, 2, [10, 20])), size (ncx2rnd (5, 2, 10, 20)))
%!assert (size (random ("ncx2", 5, 2, [10, 20])), size (ncx2rnd (5, 2, 10, 20)))
%!assert (size (random ("Normal", 5, 2, [10, 20])), size (normrnd (5, 2, 10, 20)))
%!assert (size (random ("norm", 5, 2, [10, 20])), size (normrnd (5, 2, 10, 20)))
%!assert (size (random ("Poisson", 5, [10, 20])), size (poissrnd (5, 10, 20)))
%!assert (size (random ("poiss", 5, [10, 20])), size (poissrnd (5, 10, 20)))
%!assert (size (random ("Rayleigh", 5, [10, 20])), size (raylrnd (5, 10, 20)))
%!assert (size (random ("rayl", 5, [10, 20])), size (raylrnd (5, 10, 20)))
%!assert (size (random ("Rician", 5, 1, [10, 20])), size (ricernd (5, 1, 10, 20)))
%!assert (size (random ("rice", 5, 1, [10, 20])), size (ricernd (5, 1, 10, 20)))
%!assert (size (random ("Student T", 5, [10, 20])), size (trnd (5, 10, 20)))
%!assert (size (random ("t", 5, [10, 20])), size (trnd (5, 10, 20)))
%!assert (size (random ("Triangular", 5, 2, 2, [10, 20])), size (trirnd (5, 2, 2, 10, 20)))
%!assert (size (random ("tri", 5, 2, 2, [10, 20])), size (trirnd (5, 2, 2, 10, 20)))
%!assert (size (random ("Discrete Uniform", 5, [10, 20])), size (unidrnd (5, 10, 20)))
%!assert (size (random ("unid", 5, [10, 20])), size (unidrnd (5, 10, 20)))
%!assert (size (random ("Uniform", 5, 2, [10, 20])), size (unifrnd (5, 2, 10, 20)))
%!assert (size (random ("unif", 5, 2, [10, 20])), size (unifrnd (5, 2, 10, 20)))
%!assert (size (random ("Von Mises", 5, 2, [10, 20])), size (vmrnd (5, 2, 10, 20)))
%!assert (size (random ("vm", 5, 2, [10, 20])), size (vmrnd (5, 2, 10, 20)))
%!assert (size (random ("Weibull", 5, 2, [10, 20])), size (wblrnd (5, 2, 10, 20)))
%!assert (size (random ("wbl", 5, 2, [10, 20])), size (wblrnd (5, 2, 10, 20)))

## Test input validation
%!error<random: distribution NAME must a char string.> random (1)
%!error<random: distribution NAME must a char string.> random ({"beta"})
%!error<random: distribution parameters must be numeric.> ...
%! random ("Beta", "a", 2)
%!error<random: distribution parameters must be numeric.> ...
%! random ("Beta", 5, "")
%!error<random: distribution parameters must be numeric.> ...
%! random ("Beta", 5, {2})
%!error<random: distribution parameters must be numeric.> ...
%! random ("Beta", "a", 2, 2, 10)
%!error<random: distribution parameters must be numeric.> ...
%! random ("Beta", 5, "", 2, 10)
%!error<random: distribution parameters must be numeric.> ...
%! random ("Beta", 5, {2}, 2, 10)
%!error<random: distribution parameters must be numeric.> ...
%! random ("Beta", 5, "", 2, 10)
%!error<random: chi2 distribution requires 1 parameter.> random ("chi2")
%!error<random: Beta distribution requires 2 parameters.> random ("Beta", 5)
%!error<random: Burr distribution requires 3 parameters.> random ("Burr", 5)
%!error<random: Burr distribution requires 3 parameters.> random ("Burr", 5, 2)
