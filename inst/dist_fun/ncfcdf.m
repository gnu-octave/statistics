## Copyright (C) 2022-2023 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {statistics} {@var{p} =} ncfcdf (@var{x}, @var{df1}, @var{df2}, @var{lambda})
## @deftypefnx {statistics} {@var{p} =} ncfcdf (@var{x}, @var{df1}, @var{df2}, @var{lambda}, @qcode{"upper"})
##
## Noncentral F cumulative distribution function (CDF).
##
## For each element of @var{x}, compute the cumulative distribution function
## (CDF) of the noncentral F distribution with @var{df1} and @var{df2} degrees
## of freedom and noncentrality parameter @var{lambda}.  The size of @var{p} is
## the common size of @var{x}, @var{df1}, @var{df2}, and @var{lambda}.  A scalar
## input functions as a constant matrix of the same size as the other inputs.
##
## @code{@var{p} = ncfcdf (@var{x}, @var{df1}, @var{df2}, @var{lambda}, "upper")}
## computes the upper tail probability of the noncentral F distribution with
## parameters @var{df1}, @var{df2}, and @var{lambda}, at the values in @var{x}.
##
## Further information about the noncentral F distribution can be found at
## @url{https://en.wikipedia.org/wiki/Noncentral_F-distribution}
##
## @seealso{ncfinv, ncfpdf, ncfrnd, ncfstat, fcdf}
## @end deftypefn

function p = ncfcdf (x, df1, df2, lambda, uflag)

  ## Check for valid number of input arguments
  if (nargin <  4)
    error ("ncfcdf: function called with too few input arguments.");
  endif

  ## Check for valid "upper" flag
  if (nargin > 4)
    if (! strcmpi (uflag, "upper"))
      error ("ncfcdf: invalid argument for upper tail.");
    else
      uflag = true;
    endif
  else
    uflag = false;
  endif

  ## Check for common size of X, DF1, DF2, and LAMBDA
  [err, x, df1, df2, lambda] = common_size (x, df1, df2, lambda);
  if (err > 0)
    error ("ncfcdf: X, DF1, DF2, and LAMBDA must be of common size or scalars.");
  endif

  ## Check for X, DF1, DF2, and LAMBDA being reals
  if (iscomplex (x) || iscomplex (df1) || iscomplex (df2) || iscomplex (lambda))
    error ("ncfcdf: X, DF1, DF2, and LAMBDA must not be complex.");
  endif

  ## Check for class type
  if (isa (x, "single") || isa (df1, "single") || ...
      isa (df2, "single") || isa (lambda, "single"))
    p = zeros (size (x), "single");
    c_eps = eps ("single") .^ (3/4);
  else
    p = zeros (size (x));
    c_eps = eps .^ (3/4);
  endif

  ## Find NaNs in input arguments (if any) and propagate them to p
  is_nan = isnan (x) | isnan (df1) | isnan (df1) | isnan (lambda);
  p(is_nan) = NaN;

  ## For "upper" option, force p = 1 for x <= 0, and p = 0 for x == Inf,
  ## otherwise, force p = 1 for x == Inf.
  if (uflag)
    p(x == Inf & ! is_nan) = 0;
    p(x <= 0 & ! is_nan) = 1;
  else
    p(x == Inf & ! is_nan) = 1;
  endif

  ## Find invalid values of parameters and propagate them to p as NaN
  k = (df1 <= 0 | df2 <= 0 | lambda < 0);
  p(k) = NaN;

  ## Compute central distribution (lambda == 0)
  k0 = (lambda==0);
  if (any (k0(:)))
    if (uflag)
      p(k0) = fcdf (x(k0), df1(k0), df2(k0), "upper");
    else
      p(k0) = fcdf (x(k0), df1(k0), df2(k0));
    endif
  endif

  ## Check if there are remaining elements and reset variables
  k1 = ! (k0 | k | x == Inf | x <= 0 | is_nan);
  if (! any (k1(:)))
    return;
  else
    x = x(k1);
    df1 = df1(k1);
    df2 = df2(k1);
    lambda = lambda(k1);
  endif

  ## Prepare variables
  x = x(:);
  df1 = df1(:) / 2;
  df2 = df2(:) / 2;
  lambda = lambda(:) / 2;

  ## Value passed to Beta distribution function.
  tmp = df1 .* x ./ (df2 + df1 .* x);
  logtmp = log (tmp);
  nu2const = df2 .* log (1 - tmp) - localgammaln (df2);

  ## Sum the series.  The general idea is that we are going to sum terms
  ## of the form 'poisspdf(j,lambda) .* betacdf(tmp,j+df1,df2)'
  j0 = floor (lambda(:));

  ## Compute Poisson pdf and beta cdf at the starting point
  if (uflag)
    bcdf0 = betainc (tmp, j0 + df1, df2, "upper");
  else
    bcdf0 = betacdf (tmp, j0 + df1, df2);
  endif
  ppdf0 = exp (-lambda + j0 .* log (lambda) - localgammaln (j0 + 1));

  ## Set up for loop over values less than j0
  y = ppdf0 .* bcdf0;
  ppdf = ppdf0;
  bcdf = bcdf0;
  olddy = zeros (size (lambda));
  delty = zeros (size (lambda));
  j = j0 - 1;
  ok = j >= 0;
  while (any (ok))
    ## Use recurrence relation to compute new pdf and cdf
    ppdf(ok) = ppdf(ok) .* (j(ok) + 1) ./ lambda(ok);
    if (uflag)
      bcdf(ok) = betainc (tmp(ok), j(ok) + df1(ok), df2(ok), "upper");
    else
      db = exp ((j + df1) .* logtmp + nu2const + ...
           localgammaln (j + df1 + df2) - localgammaln (j + df1 + 1));
      bcdf(ok) = bcdf(ok) + db(ok);
    endif
    delty(ok) = ppdf(ok) .* bcdf(ok);
    y(ok) = y(ok) + delty(ok);
    ## Convergence test:  change must be small and not increasing
    ok = ok & (delty > y*c_eps | abs (delty) > olddy);
    j = j - 1;
    ok = ok & j >= 0;
    olddy(ok) = abs (delty(ok));
  endwhile

  ## Set up again for loop upward from j0
  ppdf = ppdf0;
  bcdf = bcdf0;
  olddy = zeros (size (lambda));
  j = j0 + 1;
  ok = true(size(j));
  ## Set up for loop to avoid endless loop
  for jj = 1:5000
    ppdf = ppdf .* lambda ./ j;
    if (uflag)
      bcdf = betainc (tmp, j + df1, df2, "upper");
    else
      bcdf = bcdf - exp ((j + df1 - 1) .* logtmp + nu2const + ...
             localgammaln (j + df1 + df2 - 1) - localgammaln (j + df1));
    endif
    delty = ppdf.*bcdf;
    ## ok = indices not converged
    y(ok) = y(ok) + delty(ok);
    ## Convergence test:  change must be small and not increasing
    ok = ok & (delty>y*c_eps | abs(delty)>olddy);
    ## Break if all indices converged
    if (! any (ok))
      break;
    endif
    olddy(ok) = abs (delty(ok));
    j = j + 1;
  endfor
  if (jj == 5000)
    warning ("ncfcdf: no convergence.");
  endif

  ## Save returning p-value
  p(k1) = y;

endfunction

function x = localgammaln (y)
  x = Inf (size (y), class (y));
  x(! (y < 0)) = gammaln (y(! (y < 0)));
endfunction

%!demo
%! ## Plot various CDFs from the noncentral F distribution
%! x = 0:0.01:5;
%! p1 = ncfcdf (x, 2, 5, 1);
%! p2 = ncfcdf (x, 2, 5, 2);
%! p3 = ncfcdf (x, 5, 10, 1);
%! p4 = ncfcdf (x, 10, 20, 10);
%! plot (x, p1, "-r", x, p2, "-g", x, p3, "-k", x, p4, "-m")
%! grid on
%! xlim ([0, 5])
%! legend ({"df1 = 2, df2 = 5, λ = 1", "df1 = 2, df2 = 5, λ = 2", ...
%!          "df1 = 5, df2 = 10, λ = 1", "df1 = 10, df2 = 20, λ = 10"}, ...
%!         "location", "southeast")
%! title ("Noncentral F CDF")
%! xlabel ("values in x")
%! ylabel ("probability")

%!demo
%! ## Compare the noncentral F CDF with LAMBDA = 10 to the F CDF with the
%! ## same number of numerator and denominator degrees of freedom (5, 20)
%!
%! x = 0.01:0.1:10.01;
%! p1 = ncfcdf (x, 5, 20, 10);
%! p2 = fcdf (x, 5, 20);
%! plot (x, p1, "-", x, p2, "-");
%! grid on
%! xlim ([0, 10])
%! legend ({"Noncentral F(5,20,10)", "F(5,20)"}, "location", "southeast")
%! title ("Noncentral F vs F CDFs")
%! xlabel ("values in x")
%! ylabel ("probability")

## Test output
%!test
%! x = -2:0.1:2;
%! p = ncfcdf (x, 10, 1, 3);
%! assert (p([1:21]), zeros (21, 1), 1e-76);
%! assert (p(22), 0.004530737275319753, 1e-14);
%! assert (p(30), 0.255842099135669, 1e-14);
%! assert (p(41), 0.4379890998457305, 1e-14);
%!test
%! p = ncfcdf (12, 10, 3, 2);
%! assert (p, 0.9582287900447416, 1e-14);
%!test
%! p = ncfcdf (2, 3, 2, 1);
%! assert (p, 0.5731985522994989, 1e-14);
%!test
%! p = ncfcdf (2, 3, 2, 1, "upper");
%! assert (p, 0.4268014477004823, 1e-14);
%!test
%! p = ncfcdf ([3, 6], 3, 2, 5, "upper");
%! assert (p, [0.530248523596927, 0.3350482341323044], 1e-14);

## Test input validation
%!error<ncfcdf: function called with too few input arguments.> ncfcdf ()
%!error<ncfcdf: function called with too few input arguments.> ncfcdf (1)
%!error<ncfcdf: function called with too few input arguments.> ncfcdf (1, 2)
%!error<ncfcdf: function called with too few input arguments.> ncfcdf (1, 2, 3)
%!error<ncfcdf: invalid argument for upper tail.> ncfcdf (1, 2, 3, 4, "tail")
%!error<ncfcdf: invalid argument for upper tail.> ncfcdf (1, 2, 3, 4, 5)
%!error<ncfcdf: X, DF1, DF2, and LAMBDA must be of common size or scalars.> ...
%! ncfcdf (ones (3), ones (2), ones (2), ones (2))
%!error<ncfcdf: X, DF1, DF2, and LAMBDA must be of common size or scalars.> ...
%! ncfcdf (ones (2), ones (3), ones (2), ones (2))
%!error<ncfcdf: X, DF1, DF2, and LAMBDA must be of common size or scalars.> ...
%! ncfcdf (ones (2), ones (2), ones (3), ones (2))
%!error<ncfcdf: X, DF1, DF2, and LAMBDA must be of common size or scalars.> ...
%! ncfcdf (ones (2), ones (2), ones (2), ones (3))
%!error<ncfcdf: X, DF1, DF2, and LAMBDA must not be complex.> ncfcdf (i, 2, 2, 2)
%!error<ncfcdf: X, DF1, DF2, and LAMBDA must not be complex.> ncfcdf (2, i, 2, 2)
%!error<ncfcdf: X, DF1, DF2, and LAMBDA must not be complex.> ncfcdf (2, 2, i, 2)
%!error<ncfcdf: X, DF1, DF2, and LAMBDA must not be complex.> ncfcdf (2, 2, 2, i)
