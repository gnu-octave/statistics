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
## @deftypefn  {statistics} {@var{x} =} ncfinv (@var{p}, @var{df1}, @var{df2}, @var{lambda})
##
## Inverse of the noncentral F cumulative distribution function (iCDF).
##
## For each element of @var{p}, compute the quantile (the inverse of the CDF) of
## the noncentral F distribution with @var{df1} and @var{df2} degrees of freedom
## and noncentrality parameter @var{lambda}.  The size of @var{x} is the common
## size of @var{p}, @var{df1}, @var{df2}, and @var{lambda}.  A scalar input
## functions as a constant matrix of the same size as the other inputs.
##
## @code{ncfinv} uses Newton's method to converge to the solution.
##
## Further information about the noncentral F distribution can be found at
## @url{https://en.wikipedia.org/wiki/Noncentral_F-distribution}
##
## @seealso{ncfcdf, ncfpdf, ncfrnd, ncfstat, finv}
## @end deftypefn

function x = ncfinv (p, df1, df2, lambda)

  ## Check for valid number of input arguments
  if (nargin <  4)
    error ("ncfinv: function called with too few input arguments.");
  endif

  ## Check for common size of P, DF1, DF2, and LAMBDA
  [err, p, df1, df2, lambda] = common_size (p, df1, df2, lambda);
  if (err > 0)
    error ("ncfinv: P, DF1, DF2, and LAMBDA must be of common size or scalars.");
  endif

  ## Check for P, DF1, DF2, and LAMBDA being reals
  if (iscomplex (p) || iscomplex (df1) || iscomplex (df2) || iscomplex (lambda))
    error ("ncfinv: P, DF1, DF2, and LAMBDA must not be complex.");
  endif

  ## Check for class type
  if (isa (p, "single") || isa (df1, "single") || ...
      isa (df2, "single") || isa (lambda, "single"))
    x = NaN (size (p), "single");
    crit = sqrt (eps ("single"));
  else
    x = NaN (size (p), "double");
    crit = sqrt (eps ("double"));
  endif

  ## For lambda == 0, call finv
  d0 = lambda == 0;
  if (any (d0(:)))
    x(d0) = finv (p(d0), df1(d0), df2(d0));
  endif

  ## For lambda > 0 and valid dfs
  valid = df1 > 0 & df2 > 0 & lambda > 0;
  ## Force x = 0 for p == 0 ax = Inf for p ==1
  x(p == 0 & valid) = 0;
  x(p == 1 & valid) = Inf;
  ## Find remaining valid cases within the range of 0 < p < 1
  k = find (p > 0 & p < 1 & valid);
  ## Return if nothing left
  if isempty(k)
    return;
  endif

  ## Reset input variables to remaining cases
  p = p(k);
  df1 = df1(k);
  df2 = df2(k);
  lambda = lambda(k);

  ## Initialize counter
  count_limit = 100;
  count = 0;

  ## Start at the mean (if it exists)
  mu0 = df2.*(df1+lambda) ./ (df1.*max(1,df2-2));
  next = mu0;
  prev = 0;
  F = ncfcdf (mu0, df1, df2, lambda);
  while(count < count_limit)
    count += 1;
    next = (F - p) ./ ncfpdf (mu0, df1, df2, lambda);

    ## Prevent oscillations
    if (length (next) == length (prev))
      t = sign (next) == -sign (prev);
      next(t) = sign (next(t)) .* min (abs (next(t)), abs (prev(t))) / 2;
    endif

    ## Prepare for next step
    mu1 = max (mu0 / 5, min (5 * mu0, mu0 - next));

    ## Check that next step improves, otherwise abort
    F1 = ncfcdf (mu1, df1, df2, lambda);
    while (true)
      worse = (abs (F1-p) > abs (F - p) * (1 + crit)) & ...
              (abs (mu0 - mu1) > crit * mu0);
      if (! any (worse))
        break;
      endif
      mu1(worse) = 0.5 * (mu1(worse) + mu0(worse));
      F1(worse) = ncfcdf (mu1(worse), df1(worse), df2(worse), lambda(worse));
    endwhile
    x(k) = mu1;

    ## Find elements that are not converged yet
    next = mu0 - mu1;
    mask = (abs (next) > crit * abs (mu0));
    if (! any (mask))
      break;
    endif

    ## Save parameters for these elements only
    F = F1(mask);
    mu0 = mu1(mask);
    prev = next(mask);
    if (! all(mask))
      df1 = df1(mask);
      df2 = df2(mask);
      lambda = lambda(mask);
      p = p(mask);
      k = k(mask);
    endif
  endwhile

  if (count == count_limit)
    warning ("ncfinv: did not converge.");
  endif

endfunction

%!demo
%! ## Plot various iCDFs from the noncentral F distribution
%! p = 0.001:0.001:0.999;
%! x1 = ncfinv (p, 2, 5, 1);
%! x2 = ncfinv (p, 2, 5, 2);
%! x3 = ncfinv (p, 5, 10, 1);
%! x4 = ncfinv (p, 10, 20, 10);
%! plot (p, x1, "-r", p, x2, "-g", p, x3, "-k", p, x4, "-m")
%! grid on
%! ylim ([0, 5])
%! legend ({"df1 = 2, df2 = 5, λ = 1", "df1 = 2, df2 = 5, λ = 2", ...
%!          "df1 = 5, df2 = 10, λ = 1", "df1 = 10, df2 = 20, λ = 10"}, ...
%!         "location", "northwest")
%! title ("Noncentral F iCDF")
%! xlabel ("probability")
%! ylabel ("values in x")

%!demo
%! ## Compare the noncentral F iCDF with LAMBDA = 10 to the F iCDF with the
%! ## same number of numerator and denominator degrees of freedom (5, 20)
%!
%! p = 0.001:0.001:0.999;
%! x1 = ncfinv (p, 5, 20, 10);
%! x2 = finv (p, 5, 20);
%! plot (p, x1, "-", p, x2, "-");
%! grid on
%! ylim ([0, 10])
%! legend ({"Noncentral F(5,20,10)", "F(5,20)"}, "location", "northwest")
%! title ("Noncentral F vs F quantile functions")
%! xlabel ("probability")
%! ylabel ("values in x")

## Test output
%!test
%! x = [0,0.1775,0.3864,0.6395,0.9564,1.3712,1.9471,2.8215,4.3679,8.1865,Inf];
%! assert (ncfinv ([0:0.1:1], 2, 3, 1), x, 1e-4);
%!test
%! x = [0,0.7492,1.3539,2.0025,2.7658,3.7278,5.0324,6.9826,10.3955,18.7665,Inf];
%! assert (ncfinv ([0:0.1:1], 2, 3, 5), x, 1e-4);
%!test
%! x = [0,0.2890,0.8632,1.5653,2.4088,3.4594,4.8442,6.8286,10.0983,17.3736,Inf];
%! assert (ncfinv ([0:0.1:1], 1, 4, 3), x, 1e-4);
%!test
%! x = [0.078410, 0.212716, 0.288618, 0.335752, 0.367963, 0.391460];
%! assert (ncfinv (0.05, [1, 2, 3, 4, 5, 6], 10, 3), x, 1e-6);
%!test
%! x = [0.2574, 0.2966, 0.3188, 0.3331, 0.3432, 0.3507];
%! assert (ncfinv (0.05, 5, [1, 2, 3, 4, 5, 6], 3), x, 1e-4);
%!test
%! x = [1.6090, 1.8113, 1.9215, 1.9911, NaN, 2.0742];
%! assert (ncfinv (0.05, 1, [1, 2, 3, 4, -1, 6], 10), x, 1e-4);
%!test
%! assert (ncfinv (0.996, 3, 5, 8), 58.0912074080671, 2e-13);

## Test input validation
%!error<ncfinv: function called with too few input arguments.> ncfinv ()
%!error<ncfinv: function called with too few input arguments.> ncfinv (1)
%!error<ncfinv: function called with too few input arguments.> ncfinv (1, 2)
%!error<ncfinv: function called with too few input arguments.> ncfinv (1, 2, 3)
%!error<ncfinv: P, DF1, DF2, and LAMBDA must be of common size or scalars.> ...
%! ncfinv (ones (3), ones (2), ones (2), ones (2))
%!error<ncfinv: P, DF1, DF2, and LAMBDA must be of common size or scalars.> ...
%! ncfinv (ones (2), ones (3), ones (2), ones (2))
%!error<ncfinv: P, DF1, DF2, and LAMBDA must be of common size or scalars.> ...
%! ncfinv (ones (2), ones (2), ones (3), ones (2))
%!error<ncfinv: P, DF1, DF2, and LAMBDA must be of common size or scalars.> ...
%! ncfinv (ones (2), ones (2), ones (2), ones (3))
%!error<ncfinv: P, DF1, DF2, and LAMBDA must not be complex.> ncfinv (i, 2, 2, 2)
%!error<ncfinv: P, DF1, DF2, and LAMBDA must not be complex.> ncfinv (2, i, 2, 2)
%!error<ncfinv: P, DF1, DF2, and LAMBDA must not be complex.> ncfinv (2, 2, i, 2)
%!error<ncfinv: P, DF1, DF2, and LAMBDA must not be complex.> ncfinv (2, 2, 2, i)
