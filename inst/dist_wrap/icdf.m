# Copyright (C) 2023 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## desired quantile distribution function as listed in the following table.
## The last column shows the number of required parameters that should be parsed
## after @var{x} to the desired iCDF.
##
## @multitable @columnfractions 0.4 0.2 0.3
## @headitem Distribution Name @tab Abbreviation @tab Input Parameters
## @item @qcode{'Beta'} @tab @qcode{'beta'} @tab 2
## @item @qcode{'Binomial'} @tab @qcode{'bino'} @tab 2
## @item @qcode{'Birnbaum-Saunders'} @tab @qcode{'bisa'} @tab 2
## @item @qcode{'Burr'} @tab @qcode{'burr'} @tab 3
## @item @qcode{'Cauchy'} @tab @qcode{'cauchy'} @tab 2
## @item @qcode{'Chi-squared'} @tab @qcode{'chi2'} @tab 1
## @item @qcode{'Extreme Value'} @tab @qcode{'ev'} @tab 2
## @item @qcode{'Exponential'} @tab @qcode{'exp'} @tab 1
## @item @qcode{'F-Distribution'} @tab @qcode{'f'} @tab 2
## @item @qcode{'Gamma'} @tab @qcode{'gam'} @tab 2
## @item @qcode{'Geometric'} @tab @qcode{'geo'} @tab 1
## @item @qcode{'Generalized Extreme Value'} @tab @qcode{'gev'} @tab 3
## @item @qcode{'Generalized Pareto'} @tab @qcode{'gp'} @tab 3
## @item @qcode{'Gumbel'} @tab @qcode{'gumbel'} @tab 2
## @item @qcode{'Half-normal'} @tab @qcode{'hn'} @tab 2
## @item @qcode{'Hypergeometric'} @tab @qcode{'hyge'} @tab 3
## @item @qcode{'Inverse Gaussian'} @tab @qcode{'invg'} @tab 2
## @item @qcode{'Laplace'} @tab @qcode{'laplace'} @tab 2
## @item @qcode{'Logistic'} @tab @qcode{'logi'} @tab 2
## @item @qcode{'Log-Logistic'} @tab @qcode{'logl'} @tab 2
## @item @qcode{'Lognormal'} @tab @qcode{'logn'} @tab 2
## @item @qcode{'Nakagami'} @tab @qcode{'naka'} @tab 2
## @item @qcode{'Negative Binomial'} @tab @qcode{'nbin'} @tab 2
## @item @qcode{'Noncentral F-Distribution'} @tab @qcode{'ncf'} @tab 3
## @item @qcode{'Noncentral Student T'} @tab @qcode{'nct'} @tab 2
## @item @qcode{'Noncentral Chi-Squared'} @tab @qcode{'ncx2'} @tab 2
## @item @qcode{'Normal'} @tab @qcode{'norm'} @tab 2
## @item @qcode{'Poisson'} @tab @qcode{'poiss'} @tab 1
## @item @qcode{'Rayleigh'} @tab @qcode{'rayl'} @tab 1
## @item @qcode{'Rician'} @tab @qcode{'rice'} @tab 2
## @item @qcode{'Student T'} @tab @qcode{'t'} @tab 1
## @item @qcode{'location-scale T'} @tab @qcode{'tls'} @tab 3
## @item @qcode{'Triangular'} @tab @qcode{'tri'} @tab 3
## @item @qcode{'Discrete Uniform'} @tab @qcode{'unid'} @tab 1
## @item @qcode{'Uniform'} @tab @qcode{'unif'} @tab 2
## @item @qcode{'Von Mises'} @tab @qcode{'vm'} @tab 2
## @item @qcode{'Weibull'} @tab @qcode{'wbl'} @tab 2
## @end multitable
##
## @seealso{icdf, pdf, random, betainv, binoinv, bisainv, burrinv, cauchyinv,
## chi2inv, evinv, expinv, finv, gaminv, geoinv, gevinv, gpinv, gumbelinv,
## hninv, hygeinv, invginv, laplaceinv, logiinv, loglinv, logninv, nakainv,
## nbininv, ncfinv, nctinv, ncx2inv, norminv, poissinv, raylinv, riceinv, tinv,
## tlsinv, triinv, unidinv, unifinv, vminv, wblinv}
## @end deftypefn

function x = icdf (name, p, varargin)

  ## implemented functions
  persistent allDF = { ...
    {'beta'     , 'Beta'},                      @betainv,      2, ...
    {'bino'     , 'Binomial'},                  @binoinv,      2, ...
    {'bisa'     , 'Birnbaum-Saunders'},         @bisainv,      2, ...
    {'burr'     , 'Burr'},                      @burrinv,      3, ...
    {'cauchy'   , 'Cauchy'},                    @cauchyinv,    2, ...
    {'chi2'     , 'Chi-squared'},               @chi2inv,      1, ...
    {'ev'       , 'Extreme Value'},             @evinv,        2, ...
    {'exp'      , 'Exponential'},               @expinv,       1, ...
    {'f'        , 'F-Distribution'},            @finv,         2, ...
    {'gam'      , 'Gamma'},                     @gaminv,       2, ...
    {'geo'      , 'Geometric'},                 @geoinv,       1, ...
    {'gev'      , 'Generalized Extreme Value'}, @gevinv,       3, ...
    {'gp'       , 'Generalized Pareto'},        @gpinv,        3, ...
    {'gumbel'   , 'Gumbel'},                    @gumbelinv,    2, ...
    {'hn'       , 'Half-normal'},               @hninv,        2, ...
    {'hyge'     , 'Hypergeometric'},            @hygeinv,      3, ...
    {'invg'     , 'Inverse Gaussian'},          @invginv,      2, ...
    {'laplace'  , 'Laplace'},                   @laplaceinv,   2, ...
    {'logi'     , 'Logistic'},                  @logiinv,      2, ...
    {'logl'     , 'Log-Logistic'},              @loglinv,      2, ...
    {'logn'     , 'Lognormal'},                 @logninv,      2, ...
    {'naka'     , 'Nakagami'},                  @nakainv,      2, ...
    {'nbin'     , 'Negative Binomial'},         @nbininv,      2, ...
    {'ncf'      , 'Noncentral F-Distribution'}, @ncfinv,       3, ...
    {'nct'      , 'Noncentral Student T'},      @nctinv,       2, ...
    {'ncx2'     , 'Noncentral Chi-squared'},    @ncx2inv,      2, ...
    {'norm'     , 'Normal'},                    @norminv,      2, ...
    {'poiss'    , 'Poisson'},                   @poissinv,     1, ...
    {'rayl'     , 'Rayleigh'},                  @raylinv,      1, ...
    {'rice'     , 'Rician'},                    @riceinv,      2, ...
    {'t'        , 'Student T'},                 @tinv,         1, ...
    {'tls'      , 'location-scale T'},          @tlsinv,       3, ...
    {'tri'      , 'Triangular'},                @triinv,       3, ...
    {'unid'     , 'Discrete Uniform'},          @unidinv,      1, ...
    {'unif'     , 'Uniform'},                   @unifinv,      2, ...
    {'vm'       , 'Von Mises'},                 @vminv,        2, ...
    {'wbl'      , 'Weibull'},                   @wblinv,       2};

  if (! ischar (name))
    error ("icdf: distribution NAME must be a char string.");
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
  idx = cellfun (@(x)any (strcmpi (name, x)), icdfnames);

  if (any (idx))

    if (nargs == icdf_args{idx})
      ## Check that all distribution parameters are numeric
      if (! all (cellfun (@(x)isnumeric (x), (varargin))))
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
%!assert_equal (icdf ('Beta', p, 5, 2), betainv (p, 5, 2))
%!assert_equal (icdf ('beta', p, 5, 2), betainv (p, 5, 2))
%!assert_equal (icdf ('Binomial', p, 5, 2), binoinv (p, 5, 2))
%!assert_equal (icdf ('bino', p, 5, 2), binoinv (p, 5, 2))
%!assert_equal (icdf ('Birnbaum-Saunders', p, 5, 2), bisainv (p, 5, 2))
%!assert_equal (icdf ('bisa', p, 5, 2), bisainv (p, 5, 2))
%!assert_equal (icdf ('Burr', p, 5, 2, 2), burrinv (p, 5, 2, 2))
%!assert_equal (icdf ('burr', p, 5, 2, 2), burrinv (p, 5, 2, 2))
%!assert_equal (icdf ('Cauchy', p, 5, 2), cauchyinv (p, 5, 2))
%!assert_equal (icdf ('cauchy', p, 5, 2), cauchyinv (p, 5, 2))
%!assert_equal (icdf ('Chi-squared', p, 5), chi2inv (p, 5))
%!assert_equal (icdf ('chi2', p, 5), chi2inv (p, 5))
%!assert_equal (icdf ('Extreme Value', p, 5, 2), evinv (p, 5, 2))
%!assert_equal (icdf ('ev', p, 5, 2), evinv (p, 5, 2))
%!assert_equal (icdf ('Exponential', p, 5), expinv (p, 5))
%!assert_equal (icdf ('exp', p, 5), expinv (p, 5))
%!assert_equal (icdf ('F-Distribution', p, 5, 2), finv (p, 5, 2))
%!assert_equal (icdf ('f', p, 5, 2), finv (p, 5, 2))
%!assert_equal (icdf ('Gamma', p, 5, 2), gaminv (p, 5, 2))
%!assert_equal (icdf ('gam', p, 5, 2), gaminv (p, 5, 2))
%!assert_equal (icdf ('Geometric', p, 5), geoinv (p, 5))
%!assert_equal (icdf ('geo', p, 5), geoinv (p, 5))
%!assert_equal (icdf ('Generalized Extreme Value', p, 5, 2, 2), gevinv (p, 5, 2, 2))
%!assert_equal (icdf ('gev', p, 5, 2, 2), gevinv (p, 5, 2, 2))
%!assert_equal (icdf ('Generalized Pareto', p, 5, 2, 2), gpinv (p, 5, 2, 2))
%!assert_equal (icdf ('gp', p, 5, 2, 2), gpinv (p, 5, 2, 2))
%!assert_equal (icdf ('Gumbel', p, 5, 2), gumbelinv (p, 5, 2))
%!assert_equal (icdf ('gumbel', p, 5, 2), gumbelinv (p, 5, 2))
%!assert_equal (icdf ('Half-normal', p, 5, 2), hninv (p, 5, 2))
%!assert_equal (icdf ('hn', p, 5, 2), hninv (p, 5, 2))
%!assert_equal (icdf ('Hypergeometric', p, 5, 2, 2), hygeinv (p, 5, 2, 2))
%!assert_equal (icdf ('hyge', p, 5, 2, 2), hygeinv (p, 5, 2, 2))
%!assert_equal (icdf ('Inverse Gaussian', p, 5, 2), invginv (p, 5, 2))
%!assert_equal (icdf ('invg', p, 5, 2), invginv (p, 5, 2))
%!assert_equal (icdf ('Laplace', p, 5, 2), laplaceinv (p, 5, 2))
%!assert_equal (icdf ('laplace', p, 5, 2), laplaceinv (p, 5, 2))
%!assert_equal (icdf ('Logistic', p, 5, 2), logiinv (p, 5, 2))
%!assert_equal (icdf ('logi', p, 5, 2), logiinv (p, 5, 2))
%!assert_equal (icdf ('Log-Logistic', p, 5, 2), loglinv (p, 5, 2))
%!assert_equal (icdf ('logl', p, 5, 2), loglinv (p, 5, 2))
%!assert_equal (icdf ('Lognormal', p, 5, 2), logninv (p, 5, 2))
%!assert_equal (icdf ('logn', p, 5, 2), logninv (p, 5, 2))
%!assert_equal (icdf ('Nakagami', p, 5, 2), nakainv (p, 5, 2))
%!assert_equal (icdf ('naka', p, 5, 2), nakainv (p, 5, 2))
%!assert_equal (icdf ('Negative Binomial', p, 5, 2), nbininv (p, 5, 2))
%!assert_equal (icdf ('nbin', p, 5, 2), nbininv (p, 5, 2))
%!assert_equal (icdf ('Noncentral F-Distribution', p, 5, 2, 2), ncfinv (p, 5, 2, 2))
%!assert_equal (icdf ('ncf', p, 5, 2, 2), ncfinv (p, 5, 2, 2))
%!assert_equal (icdf ('Noncentral Student T', p, 5, 2), nctinv (p, 5, 2))
%!assert_equal (icdf ('nct', p, 5, 2), nctinv (p, 5, 2))
%!assert_equal (icdf ('Noncentral Chi-Squared', p, 5, 2), ncx2inv (p, 5, 2))
%!assert_equal (icdf ('ncx2', p, 5, 2), ncx2inv (p, 5, 2))
%!assert_equal (icdf ('Normal', p, 5, 2), norminv (p, 5, 2))
%!assert_equal (icdf ('norm', p, 5, 2), norminv (p, 5, 2))
%!assert_equal (icdf ('Poisson', p, 5), poissinv (p, 5))
%!assert_equal (icdf ('poiss', p, 5), poissinv (p, 5))
%!assert_equal (icdf ('Rayleigh', p, 5), raylinv (p, 5))
%!assert_equal (icdf ('rayl', p, 5), raylinv (p, 5))
%!assert_equal (icdf ('Rician', p, 5, 1), riceinv (p, 5, 1))
%!assert_equal (icdf ('rice', p, 5, 1), riceinv (p, 5, 1))
%!assert_equal (icdf ('Student T', p, 5), tinv (p, 5))
%!assert_equal (icdf ('t', p, 5), tinv (p, 5))
%!assert_equal (icdf ('location-scale T', p, 5, 1, 2), tlsinv (p, 5, 1, 2))
%!assert_equal (icdf ('tls', p, 5, 1, 2), tlsinv (p, 5, 1, 2))
%!assert_equal (icdf ('Triangular', p, 5, 2, 2), triinv (p, 5, 2, 2))
%!assert_equal (icdf ('tri', p, 5, 2, 2), triinv (p, 5, 2, 2))
%!assert_equal (icdf ('Discrete Uniform', p, 5), unidinv (p, 5))
%!assert_equal (icdf ('unid', p, 5), unidinv (p, 5))
%!assert_equal (icdf ('Uniform', p, 5, 2), unifinv (p, 5, 2))
%!assert_equal (icdf ('unif', p, 5, 2), unifinv (p, 5, 2))
%!assert_equal (icdf ('Von Mises', p, 5, 2), vminv (p, 5, 2))
%!assert_equal (icdf ('vm', p, 5, 2), vminv (p, 5, 2))
%!assert_equal (icdf ('Weibull', p, 5, 2), wblinv (p, 5, 2))
%!assert_equal (icdf ('wbl', p, 5, 2), wblinv (p, 5, 2))

## Test input validation
%!error<icdf: distribution NAME must be a char string.> icdf (1)
%!error<icdf: distribution NAME must be a char string.> icdf ({'beta'})
%!error<icdf: P must be numeric.> icdf ('beta', {[1 2 3 4 5]})
%!error<icdf: P must be numeric.> icdf ('beta', 'text')
%!error<icdf: values in P must be real.> icdf ('beta', 1+i)
%!error<icdf: distribution parameters must be numeric.> ...
%! icdf ('Beta', p, 'a', 2)
%!error<icdf: distribution parameters must be numeric.> ...
%! icdf ('Beta', p, 5, '')
%!error<icdf: distribution parameters must be numeric.> ...
%! icdf ('Beta', p, 5, {2})
%!error<icdf: chi2 distribution requires 1 parameter.> icdf ('chi2', p)
%!error<icdf: Beta distribution requires 2 parameters.> icdf ('Beta', p, 5)
%!error<icdf: Burr distribution requires 3 parameters.> icdf ('Burr', p, 5)
%!error<icdf: Burr distribution requires 3 parameters.> icdf ('Burr', p, 5, 2)
