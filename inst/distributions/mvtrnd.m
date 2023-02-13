## Copyright (C) 2012  Arno Onken <asnelt@asnelt.org>, Iñigo Urteaga
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
## @deftypefn  {statistics} {@var{r} =} mvtrnd (@var{sigma}, @var{df})
## @deftypefnx {statistics} {@var{r} =} mvtrnd (@var{sigma}, @var{df}, @var{n})
##
## Random vectors from the multivariate Student's t distribution.
##
## @subheading Arguments
##
## @itemize @bullet
## @item
## @var{sigma} is the matrix of correlation coefficients. If there are any
## non-unit diagonal elements then @var{sigma} will be normalized, so that the
## resulting covariance of the obtained samples @var{r} follows:
## @code{cov (r) = df/(df-2) * sigma ./ (sqrt (diag (sigma) * diag (sigma)))}.
## In order to obtain samples distributed according to a standard multivariate
## t-distribution, @var{sigma} must be equal to the identity matrix. To generate
## multivariate t-distribution samples @var{r} with arbitrary covariance matrix
## @var{sigma}, the following scaling might be used:
## @code{r = mvtrnd (sigma, df, n) * diag (sqrt (diag (sigma)))}.
##
## @item
## @var{df} is the degrees of freedom for the multivariate t-distribution.
## @var{df} must be a vector with the same number of elements as samples to be
## generated or be scalar.
##
## @item
## @var{n} is the number of rows of the matrix to be generated. @var{n} must be
## a non-negative integer and corresponds to the number of samples to be
## generated.
## @end itemize
##
## @subheading Return values
##
## @itemize @bullet
## @item
## @var{r} is a matrix of random samples from the multivariate t-distribution
## with @var{n} row samples.
## @end itemize
##
## @subheading Examples
##
## @example
## @group
## sigma = [1, 0.5; 0.5, 1];
## df = 3;
## n = 10;
## r = mvtrnd (sigma, df, n);
## @end group
##
## @group
## sigma = [1, 0.5; 0.5, 1];
## df = [2; 3];
## n = 2;
## r = mvtrnd (sigma, df, 2);
## @end group
## @end example
##
## @subheading References
##
## @enumerate
## @item
## Wendy L. Martinez and Angel R. Martinez. @cite{Computational Statistics
## Handbook with MATLAB}. Appendix E, pages 547-557, Chapman & Hall/CRC, 2001.
##
## @item
## Samuel Kotz and Saralees Nadarajah. @cite{Multivariate t Distributions and
## Their Applications}. Cambridge University Press, Cambridge, 2004.
## @end enumerate
##
## @seealso{mvtcdf, mvtcdfqmc, mvtpdf}
## @end deftypefn

function r = mvtrnd (sigma, df, n)

  # Check arguments
  if (nargin < 2)
    print_usage ();
  endif
  [jnk, p] = cholcov (sigma); # This is a more robust check for positive definite
  if (! ismatrix (sigma) || any (any (sigma != sigma')) || (p != 0))
    error ("mvtrnd: SIGMA must be a positive definite matrix.");
  endif

  if (!isvector (df) || any (df <= 0))
    error ("mvtrnd: DF must be a positive scalar or vector.");
  endif
  df = df(:);

  if (nargin > 2)
    if (! isscalar (n) || n < 0 | round (n) != n)
      error ("mvtrnd: N must be a non-negative integer.")
    endif
    if (isscalar (df))
      df = df * ones (n, 1);
    else
      if (length (df) != n)
        error ("mvtrnd: N must match the length of DF.")
      endif
    endif
  else
    n = length (df);
  endif

  # Normalize sigma
  if (any (diag (sigma) != 1))
    sigma = sigma ./ sqrt (diag (sigma) * diag (sigma)');
  endif

  # Dimension
  d = size (sigma, 1);
  # Draw samples
  y = mvnrnd (zeros (1, d), sigma, n);
  u = repmat (chi2rnd (df), 1, d);
  r = y .* sqrt (repmat (df, 1, d) ./ u);
endfunction

%!test
%! sigma = [1, 0.5; 0.5, 1];
%! df = 3;
%! n = 10;
%! r = mvtrnd (sigma, df, n);
%! assert (size (r), [10, 2]);

%!test
%! sigma = [1, 0.5; 0.5, 1];
%! df = [2; 3];
%! n = 2;
%! r = mvtrnd (sigma, df, 2);
%! assert (size (r), [2, 2]);
