## Copyright (C) 2015 Nir Krakauer <nkrakauer@ccny.cuny.edu>
## Copyright (C) 2023 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {statistics} {@var{y} =} mvtpdf (@var{x}, @var{rho}, @var{df})
##
## Multivariate Student's t probability density function (PDF).
##
## @subheading Arguments
##
## @itemize @bullet
## @item
## @var{x} are the points at which to find the probability, where each row
## corresponds to an observation. (@var{n} by @var{d} matrix)
##
## @item
## @var{rho} is the correlation matrix. (@var{d} by @var{d} symmetric positive
## definite matrix)
##
## @item
## @var{df} is the degrees of freedom. (scalar or @var{n} vector)
##
## @end itemize
##
## The distribution is assumed to be centered (zero mean).
##
## @subheading Return values
##
## @itemize @bullet
## @item
## @var{y} is the probability density for each row of @var{x}.
## (@var{n} by 1 vector)
## @end itemize
##
## @subheading Examples
##
## @example
## @group
## x = [1 2];
## rho = [1.0 0.5; 0.5 1.0];
## df = 4;
## y = mvtpdf (x, rho, df)
## @end group
## @end example
##
## @subheading References
##
## @enumerate
## @item
## Michael Roth, On the Multivariate t Distribution, Technical report from
## Automatic Control at Linkoepings universitet,
## @url{http://users.isy.liu.se/en/rt/roth/student.pdf}
## @end enumerate
##
## @seealso{mvtcdf, mvtcdfqmc, mvtrnd}
## @end deftypefn

function y = mvtpdf (x, rho, df)

  if (nargin != 3)
    print_usage ();
  endif

  # Dimensions
  d = size (rho, 1);
  n = size (x, 1);

  # Check parameters
  if (size (x, 2) != d)
    error ("mvtpdf: x must have the same number of columns as rho.");
  endif
  if (! isscalar (df) && (! isvector (df) || numel (df) != n))
    error (strcat (["mvtpdf: DF must be a scalar or a vector with the"], ...
                   [" same number of rows as X."]));
  endif
  if (d < 1 || size (rho, 2) != d || ! issymmetric (rho))
    error ("mvtpdf: SIGMA must be nonempty and symmetric.");
  endif

  try
    U = chol (rho);
  catch
    error ("mvtpdf: rho must be positive definite");
  end_try_catch

  df = df(:);
  sqrt_det_sigma = prod(diag(U)); #square root of determinant of rho

  ## Scale factor for PDF
  c = (gamma((df+d)/2) ./ gamma(df/2)) ./ (sqrt_det_sigma * (df*pi).^(d/2));
  #note: sumsq(U' \ x') is equivalent to the quadratic form x*inv(rho)*x'
  y = c ./ ((1 + sumsq(U' \ x') ./ df') .^ ((df' + d)/2))';


endfunction

%!demo
%! ## Compute the pdf of a multivariate t distribution with correlation
%! ## parameters rho = [1 .4; .4 1] and 2 degrees of freedom.
%!
%! rho = [1, 0.4; 0.4, 1];
%! df = 2;
%! [X1, X2] = meshgrid (linspace (-2, 2, 25)', linspace (-2, 2, 25)');
%! X = [X1(:), X2(:)];
%! y = mvtpdf (X, rho, df);
%! surf (X1, X2, reshape (y, 25, 25));
%! title ("Bivariate Student's t probability density function");

## Test results verified with R mvtnorm package dmvt function
## dmvt(x = c(0,0), rho = diag(2), log = FALSE)
%!assert (mvtpdf ([0 0], eye(2), 1), 0.1591549, 1E-7)
## dmvt(x = c(1,0), rho = matrix(c(1, 0.5, 0.5, 1), nrow=2, ncol=2), df = 2, log = FALSE)
%!assert (mvtpdf ([1 0], [1 0.5; 0.5 1], 2), 0.06615947, 1E-7)
## dmvt(x = c(1,0.4,0), rho = matrix(c(1, 0.5, 0.3, 0.5, 1, 0.6, 0.3, 0.6, ...
## 1), nrow=3, ncol=3), df = 5, log = FALSE); dmvt(x = c(1.2,0.5,0.5), ...
## rho = matrix(c(1, 0.5, 0.3, 0.5, 1, 0.6, 0.3, 0.6, 1), nrow=3, ncol=3), ...
## df = 6, log = FALSE); dmvt(x = c(1.4,0.6,1), rho = matrix(c(1, 0.5, 0.3,...
## 0.5, 1, 0.6, 0.3, 0.6, 1), nrow=3, ncol=3), df = 7, log = FALSE)
%!assert (mvtpdf ([1 0.4 0; 1.2 0.5 0.5; 1.4 0.6 1], ...
%! [1 0.5 0.3; 0.5 1 0.6; 0.3 0.6 1], [5 6 7]), ...
%! [0.04713313 0.03722421 0.02069011]', 1E-7)
