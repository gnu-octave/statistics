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
## @deftypefn {Function File} @var{ndim} = barttest (@var{x})
## @deftypefnx {Function File} @var{ndim} = barttest (@var{x}, @var{alpha})
## @deftypefnx {Function File} [@var{ndim}, @var{pval}] = barttest (@var{x}, @var{alpha})
## @deftypefnx {Function File} [@var{ndim}, @var{pval}, @var{chisq}] = barttest (@var{x}, @var{alpha})
##
## Bartlett's test of sphericity for correlation.
##
## It compares an observed correlation matrix to the identity matrix in order to
## check if there is a certain redundancy between the variables that we can
## summarize with a few number of factors.  A statistically significant test
## shows that the variables (columns) in @var{x} are correlated, thus it makes
## sense to perform some dimensionality reduction of the data in @var{x}.
##
## @code{@var{ndim} = barttest (@var{x}, @var{alpha})} returns the number of
## dimensions necessary to explain the nonrandom variation in the data matrix
## @var{x} at the @var{alpha} significance level. @var{alpha} is an optional
## input argument and, when not provided, it is 0.05 by default.
##
## @code{[@var{ndim}, @var{pval}, @var{chisq}] = barttest (@dots{})} also
## returns the significance values @var{pval} for the hypothesis test for each
## dimension as well as the associated chi^2 values in @var{chisq}
##
## @end deftypefn

function [ndim, pval, chisq] = barttest (x, alpha);

  ## Check for valid number of input arguments
  if (nargin < 1 || nargin >2)
    error ("barttest: invalid number of input arguments.");
  endif
  ## Check for NaN values in X
  if (any (isnan( x(:))))
    error ("barttest: NaN values in input are not allowed.");
  endif
  ## Add default value for alpha if not supplied
  if (nargin == 1)
    alpha = 0.05;
  endif
  ## Check for valid value of alpha
  if (! isscalar (alpha) || ! isnumeric (alpha) || alpha <= 0 || alpha >= 1)
     error ("barttest: wrong value for alpha.");
  endif
  ## Check size of data
  [row, col] = size (x);
  if (col <= 1 || row <= 1)
    error ("barttest: not enough data in X.");
  endif
  ## Compute the eigenvalues of X in a more efficient way
  latent = sort ((svd (x - repmat (mean (x, 1), row, 1)) .^ 2) / (row - 1));
  ## The degrees of ffredom should be N-1, where N is the sample size
  row -= 1;
  k = (0:col - 2)';
  pk = col - k;
  loglatent = flipud (cumsum (log (latent)));
  ## Compute the chi-square statistic
  logsum = log (flipud ((latent(1) + cumsum (latent(2:col))) ./ flipud(pk)));
  chisq = (pk .* logsum - loglatent(1:col - 1)) * row;
  ## Calculate the degrees of freedom
  df = (pk - 1) .* (pk + 2) / 2;
  ## Find the corresponding p-values
  pval = 1 - chi2cdf (chisq, df);
  ## Get ndim
  dim  = min (find (pval > alpha));
  if (isempty (dim))
   ndim = col;
   return;
  endif
  if (dim == 1)
   ndim = NaN;
   warning ("barttest: heuristics are violated.");
  else
   ndim = dim - 1;
  endif
endfunction

## Test input validation
%!error<barttest: invalid number of input arguments.> barttest ()
%!error<barttest: NaN values in input are not allowed.> barttest ([2,NaN;3,4])
%!error<barttest: wrong value for alpha.> barttest (ones (30, 4), "alpha")
%!error<barttest: wrong value for alpha.> barttest (ones (30, 4), 0)
%!error<barttest: wrong value for alpha.> barttest (ones (30, 4), 1.2)
%!error<barttest: wrong value for alpha.> barttest (ones (30, 4), [0.2, 0.05])
%!error<barttest: not enough data in X.> barttest (ones (30, 1))
%!error<barttest: not enough data in X.> barttest (ones (30, 1), 0.05)
## Test results
%!test
%! x = [2, 3, 4, 5, 6, 7, 8, 9; 1, 2, 3, 4, 5, 6, 7, 8]';
%! [ndim, pval, chisq] = barttest (x);
%! assert (ndim, 2);
%! assert (pval, 0);
%! ## assert (chisq, 512.0558, 1e-4); Result differs between octave 6 and 7 ?
%!test
%! x = [0.53767,  0.62702,   -0.10224,   -0.25485,   1.4193,   1.5237  ; ...
%!      1.8339,   1.6452,    -0.24145,   -0.23444,   0.29158,  0.1634  ; ...
%!     -2.2588,  -2.1351,     0.31286,    0.39396,   0.19781,  0.20995 ; ...
%!      0.86217,  1.0835,     0.31286,    0.46499,   1.5877,   1.495   ; ...
%!      0.31877,  0.38454,   -0.86488,   -0.63839,  -0.80447, -0.7536  ; ...
%!     -1.3077,  -1.1487,    -0.030051,  -0.017629,  0.69662,  0.60497 ; ...
%!     -0.43359, -0.32672,   -0.16488,   -0.37364,   0.83509,  0.89586 ; ...
%!      0.34262,  0.29639,    0.62771,    0.51672,  -0.24372, -0.13698 ; ...
%!      3.5784,   3.5841,     1.0933,     0.93258,   0.21567,  0.455   ; ...
%!      2.7694,   2.6307,     1.1093,     1.4298,   -1.1658,  -1.1816  ; ...
%!     -1.3499,  -1.2111,    -0.86365,   -0.94186,  -1.148,   -1.4381  ; ...
%!      3.0349,   2.8428,     0.077359,   0.18211,   0.10487, -0.014613; ...
%!      0.7254,   0.56737,   -1.2141,    -1.2291,    0.72225,  0.90612 ; ...
%!     -0.063055,-0.17662,   -1.1135,    -0.97701,   2.5855,   2.4084  ; ...
%!      0.71474,  0.29225,   -0.0068493, -0.11468,  -0.66689, -0.52466 ; ...
%!     -0.20497, -7.8874e-06, 1.5326,     1.3195,    0.18733,  0.20296 ; ...
%!     -0.12414, -0.077029,  -0.76967,   -0.96262,  -0.082494, 0.121   ; ...
%!      1.4897,   1.3683,     0.37138,    0.43653,  -1.933,   -2.1903  ; ...
%!      1.409,    1.5882,    -0.22558,   -0.24835,  -0.43897, -0.46247 ; ...
%!      1.4172,   1.1616,     1.1174,     1.0785,   -1.7947,  -1.9471  ];
%! [ndim, pval, chisq] = barttest (x);
%! assert (ndim, 3);
%! assert (pval, [0; 0; 0; 0.52063; 0.34314], 1e-5);
%! chisq_out = [251.6802; 210.2670; 153.1773; 4.2026; 2.1392];
%! assert (chisq, chisq_out, 1e-4);
