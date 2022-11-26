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
## @deftypefn {Function File} @var{x} = ncx2inv (@var{p}, @var{df}, @var{delta})
##
## Inverse of the non-central T cumulative distribution function (cdf).
##
## @code{@var{x} = nctinv (@var{p}, @var{df}, @var{delta})} returns the inverse
## of the noncentral T cdf with @var{df} degrees of freedom and noncentrality
## parameter @var{delta}, at the probabilities of @var{p}.
##
## The size of @var{x} is the common size of @var{df} and @var{delta}.  A scalar
## input functions as a constant matrix of the same size as the other inputs.
##
## @code{nctinv} uses Newton's method to converge to the solution.
##
## @seealso{nctcdf, nctpdf, nctrnd, nctstat}
## @end deftypefn

function x = nctinv (p, df, delta)

  ## Check for valid input arguments
  if (nargin <  3)
    error ("nctinv: too few input arguments.");
  endif

  ## Check and fix size of input arguments
  [err, p, df, delta] = common_size (p, df, delta);
  if (err > 0)
    error ("nctinv: input size mismatch.");
  endif

  ## Initialize x
  if (isa (p, "single") || isa (df, "single") || isa (delta, "single"))
    x = NaN (size (p), "single");
    crit = sqrt (eps ("single"));
  else
    x = NaN (size (p), "double");
    crit = sqrt (eps ("double"));
  endif

  ## For delta == 0, call chi2inv
  d0 = delta == 0;
  if (any (d0(:)))
    x(d0) = tinv (p(d0), df(d0));
    ## If delta == 0 for all entries, then return
    if (all (d0(:)))
      return;
    endif
  endif

  ## For all valid entries
  valid = df > 0 & ! isnan (delta) & ! isinf (delta);
  ## Force x = -Inf for p == 0 and x = Inf for p == 1
  x(p == 0 & valid) = -Inf;
  x(p == 1 & valid) = Inf;
  ## Find valid samples within the range of 0 < p < 1
  k = find (p > 0 & p < 1 & valid);
  p_k = p(k);
  df_k = df(k);
  delta_k = delta(k);

  ## Initialize counter
  count_limit = 100;
  count = 0;

  ## Supply a starting guess for the iteration with norminv
  x_k = norminv (p_k, delta_k, 1);
  h_k = ones (size (x_k), class (x_k));

  ## Start iteration with a break out loop
  F =  nctcdf (x_k, df_k, delta_k);
  while (any (abs (h_k) > crit * abs (x_k)) && ...
         max (abs (h_k)) > crit && count < count_limit)
    count = count + 1;
    h_k = (F - p_k) ./ nctpdf (x_k, df_k, delta_k);
    ## Prevent Infs - NaNs
    infnan = isinf(h_k) | isnan(h_k);
    if (any (infnan(:)))
      h_k(infnan) = x_k(infnan) / 10;
    endif
    ## Prepare for next step
    xnew = max (-5 * abs (x_k), min (5 * abs (x_k), x_k - h_k));
    ## Check that next step improves, otherwise abort
    Fnew = nctcdf (xnew, df_k, delta_k);
    while (true)
       worse = (abs (Fnew - p_k) > abs (F - p_k) * (1 + crit)) & ...
               (abs (x_k - xnew) > crit * abs (x_k));
       if (! any (worse))
         break;
       endif
       xnew(worse) = 0.5 * (xnew(worse) + x_k(worse));
       Fnew(worse) = nctcdf (xnew(worse), df_k(worse), delta_k(worse));
    endwhile
    x_k = xnew;
    F = Fnew;
  endwhile

  ## Return the converged value(s).
  x(k) = x_k;

  if (count == count_limit)
    warning ("nctinv: did not converge.");
    fprintf ("nctinv: Last Step: %13.8f\n", max (abs (h_k)));
  endif

endfunction

## Input validation tests
%!error<nctinv: too few input arguments.> p = nctinv ();
%!error<nctinv: too few input arguments.> p = nctinv (1);
%!error<nctinv: input size mismatch.> p = nctinv (1, [4, 3], [3, 4, 5]);

## Output validation tests
%!test
%! x = [-Inf,-0.3347,0.1756,0.5209,0.8279,1.1424,1.5021,1.9633,2.6571,4.0845,Inf];
%! assert (nctinv ([0:0.1:1], 2, 1), x, 1e-4);
%!test
%! x = [-Inf,1.5756,2.0827,2.5343,3.0043,3.5406,4.2050,5.1128,6.5510,9.6442,Inf];
%! assert (nctinv ([0:0.1:1], 2, 3), x, 1e-4);
%!test
%! x = [-Inf,2.2167,2.9567,3.7276,4.6464,5.8455,7.5619,10.3327,15.7569,31.8159,Inf];
%! assert (nctinv ([0:0.1:1], 1, 4), x, 1e-4);
%!test
%! x = [1.7791   1.9368   2.0239   2.0801   2.1195   2.1489];
%! assert (nctinv (0.05, [1, 2, 3, 4, 5, 6], 4), x, 1e-4);
%!test
%! x = [-0.7755, 0.3670, 1.2554, 2.0239, 2.7348, 3.4154];
%! assert (nctinv (0.05, 3, [1, 2, 3, 4, 5, 6]), x, 1e-4);
%!test
%! x = [-0.7183, 0.3624, 1.2878, 2.1195, -3.5413, 3.6430];
%! assert (nctinv (0.05, 5, [1, 2, 3, 4, -1, 6]), x, 1e-4);
%!test
%! assert (nctinv (0.996, 5, 8), 30.02610554063658, 1e-14);
