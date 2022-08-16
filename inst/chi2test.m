## Copyright (C) 2022 Andreas Bertsatos <abertsatos@biol.uoa.gr>
## Copyright (C) 1995-2017 Kurt Hornik
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
## @deftypefn {Function File} @var{pval} = chi2test (@var{x})
## @deftypefnx {Function File} [@var{pval}, @var{chisq}] = chi2test (@var{x})
## @deftypefnx {Function File} [@var{pval}, @var{chisq}, @var{dF}] = chi2test (@var{x})
##
## Perform a chi-squared test (for independence or homogeneity).
##
##
##
## @end deftypefn

function [pval, chisq, df] = chi2test (x)
  ## Check input argument
  if (nargin != 1)
    print_usage ();
  endif
  if (isvector (x))
    error ("chi2test: X must be a matrix.");
  endif
  if (! isreal (x))
    error ("chi2test: values in X must be real numbers.");
  endif
  if (any (isnan (x(:))))
    error ("chi2test: X must not have missing values (NaN).");
  endif
  ## Get size of contigency table
  [row, col] = size (x);
  ## Calaculate degrees of freedom
  df = (row - 1) * (col -1);
  ## Calculate total sample size
  n = sum (sum (x));
  ## Calaculate expected values
  E = sum (x')' * sum (x) / n;
  ## Check that expected values are >= 5
  if (any (any (E < 5)))
    warning ("chi2test: Expected values less than 5.");
  endif
  ## Calculate chi-squared and p-value
  cells = ((x - E) .^2) ./ E;
  chisq = sum (sum (cells));
  pval = 1 - chi2cdf (chisq, df);
  ## Print results if no output requested
  if (nargout == 0)
    printf ("p-val = %f with chi^2 statistic = %f and d.f. = %d.\n", ...
            pval, chisq, df);
  endif
endfunction

## Test input
%!error chi2test ();
%!error chi2test ([1,2,3,4,5]);
%!error chi2test ([1,2;2;1+3i]);
%!error chi2test ([NaN,6;34,12]);
%!warning p = chi2test (ones (3,3));
%!test
%! x = [11, 3, 8; 2, 9, 14; 12, 13, 28];
%! p = chi2test (x);
%! assert (p, 0.017787, 1e-6);
%!test
%! x = [11, 3, 8; 2, 9, 14; 12, 13, 28];
%! [p, chisq] = chi2test (x);
%! assert (chisq, 11.9421, 1e-4);
%!test
%! x = [11, 3, 8; 2, 9, 14; 12, 13, 28];
%! [p, chisq, df] = chi2test (x);
%! assert (df, 4);

