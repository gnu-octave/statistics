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
## @deftypefn {statistics} [@var{yCDF}, @var{xCDF}, @var{n}, @var{emsg}, @var{eid}] = cdfcalc (@var{x})
##
## Calculate an empirical cumulative distribution function.
##
## @code{[@var{yCDF}, @var{xCDF}] = cdfcalc (@var{x})} calculates an empirical
## cumulative distribution function (CDF) of the observations in the data sample
## vector @var{x}.  @var{x} may be a row or column vector, and represents a
## random sample of observations from some underlying distribution.  On return
## @var{xCDF} is the set of @var{x} values at which the CDF increases.
## At XCDF(i), the function increases from YCDF(i) to YCDF(i+1).
##
## @code{[@var{yCDF}, @var{xCDF}, @var{n}] = cdfcalc (@var{x})} also returns
## @var{n}, the sample size.
##
## @code{[@var{yCDF}, @var{xCDF}, @var{n}, @var{emsg}, @var{eid}] = cdfcalc
## (@var{x})} also returns an error message and error id if @var{x} is not a
## vector or if it contains no values other than NaN.
##
## @seealso{cdfplot}
## @end deftypefn

function [yCDF, xCDF, n, emsg, eid] = cdfcalc (x)
  ## Check number of input and output argument
  narginchk (1,1);
  nargoutchk (2,5);
  ## Add defaults
  yCDF = [];
  xCDF = [];
  n = 0;
  ## Check that x is a vector
  if (! isvector (x))
    warning ("cdfcalc: vector required as input.");
    emsg = "VectorRequired";
    eid = "VectorRequired";
    return
  endif
  ## Remove NaNs and check if there are remaining data to calculate ecdf
  x = x(! isnan (x));
  n = length (x);
  if (n == 0)
    warning ("cdfcalc: not enough data.");
    emsg = "NotEnoughData";
    eid = "NotEnoughData";
    return
  endif
  ## Sort data in ascending order
  x = sort (x(:));
  ## Get cumulative sums
  yCDF = (1:n)' / n;
  ## Remove duplicates, keep the last one
  keep_idx = ([diff(x(:)); 1] > 0);
  xCDF = x(keep_idx);
  yCDF = [0; yCDF(keep_idx)];
  emsg = '';
  eid = '';
endfunction

%!test
%! x = [2, 4, 3, 2, 4, 3, 2, 5, 6, 4];
%! [yCDF, xCDF, n, emsg, eid] = cdfcalc (x);
%! assert (yCDF, [0, 0.3, 0.5, 0.8, 0.9, 1]');
%! assert (xCDF, [2, 3, 4, 5, 6]');
%! assert (n, 10);
%!shared x
%! x = [2, 4, 3, 2, 4, 3, 2, 5, 6, 4];
%!error yCDF = cdfcalc (x);
%!error [yCDF, xCDF] = cdfcalc ();
%!error [yCDF, xCDF] = cdfcalc (x, x);
%!warning [yCDF, xCDF] = cdfcalc (ones(10,2));
