## Copyright (C) 2022 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn {Function File} @var{x} = ncfinv (@var{p}, @var{df1}, @var{df2}, @var{delta})
##
## Inverse of the non-central F cumulative distribution function (cdf).
##
## @code{@var{x} = ncfinv (@var{p}, @var{df1}, @var{df2}, @var{delta})}
## the inverse of the noncentral F distribution with @var{df1} numerator degrees
## of freedom, @var{df2} denumerator degrees of freedom, and noncentrality
## parameter @var{delta}, at the probabilities of @var{p}.
##
## The size of @var{x} is the common size of @var{df} and @var{delta}.  A scalar
## input functions as a constant matrix of the same size as the other inputs.
##
## @code{ncfinv} uses Newton's method to converge to the solution.
##
## @seealso{ncfcdf, ncfpdf, ncfrnd, ncfstat}
## @end deftypefn

function x = ncfinv (p, df1, df2, delta)

  ## Check for valid input arguments
  if (nargin <  4)
    error ("ncfinv: too few input arguments.");
  endif

  ## Check and fix size of input arguments
  [err, p, df1, df2, delta] = common_size (p, df1, df2, delta);
  if (err > 0)
    error ("ncfinv: input size mismatch.");
  endif

  ## Initialize x
  if (isa (p, "single") || isa (df1, "single") || ...
      isa (df2, "single") || isa (delta, "single"))
    x = NaN (size (p), "single");
    crit = sqrt (eps ("single"));
  else
    x = NaN (size (p), "double");
    crit = sqrt (eps ("double"));
  endif

  ## For delta == 0, call finv
  d0 = delta == 0;
  if (any (d0(:)))
    x(d0) = finv (p(d0), df1(d0), df2(d0));
  endif

  ## For delta > 0 and valid dfs
  valid = df1 > 0 & df2 > 0 & delta > 0;
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
  delta = delta(k);

  ## Initialize counter
  count_limit = 100;
  count = 0;

  ## Start at the mean (if it exists)
  mu0 = df2.*(df1+delta) ./ (df1.*max(1,df2-2));
  next = mu0;
  prev = 0;
  F = ncfcdf (mu0, df1, df2, delta);
  while(count < count_limit)
    count += 1;
    next = (F - p) ./ ncfpdf (mu0, df1, df2, delta);

    ## Prevent oscillations
    if (length (next) == length (prev))
      t = sign (next) == -sign (prev);
      next(t) = sign (next(t)) .* min (abs (next(t)), abs (prev(t))) / 2;
    endif

    ## Prepare for next step
    mu1 = max (mu0 / 5, min (5 * mu0, mu0 - next));

    ## Check that next step improves, otherwise abort
    F1 = ncfcdf (mu1, df1, df2, delta);
    while (true)
      worse = (abs (F1-p) > abs (F - p) * (1 + crit)) & ...
              (abs (mu0 - mu1) > crit * mu0);
      if (! any (worse))
        break;
      endif
      mu1(worse) = 0.5 * (mu1(worse) + mu0(worse));
      F1(worse) = ncfcdf (mu1(worse), df1(worse), df2(worse), delta(worse));
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
      delta = delta(mask);
      p = p(mask);
      k = k(mask);
    endif
  endwhile

  if (count == count_limit)
    warning ("ncfinv: did not converge.");
    fprintf ("ncfinv: Last Step: %13.8f\n", next);
  endif

endfunction

## Input validation tests
%!error<ncfinv: too few input arguments.> p = ncfinv ();
%!error<ncfinv: too few input arguments.> p = ncfinv (1);
%!error<ncfinv: too few input arguments.> p = ncfinv (1, 2);
%!error<ncfinv: too few input arguments.> p = ncfinv (1, 2, 3);
%!error<ncfinv: input size mismatch.> p = ncfinv (1, [4, 3], [3, 4, 5], 3);

## Output validation tests
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
