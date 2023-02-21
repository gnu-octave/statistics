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
## @deftypefn  {statistics} {@var{y} =} mvnpdf (@var{x}, @var{mu}, @var{sigma})
##
## Multivariate normal probability density function (PDF).
##
## @code{@var{y} = mvnpdf (@var{x})} returns the probability density of the
## multivariate normal distribution with zero mean and identity covariance
## matrix, evaluated at each row of @var{x}.  Rows of the N-by-D matrix @var{x}
## correspond to observations orpoints, and columns correspond to variables or
## coordinates.  @var{y} is an N-by-1 vector.
##
## @code{@var{y} = mvnpdf (@var{x}, @var{mu})} returns the density of the
## multivariate normal distribution with mean MU and identity covariance matrix,
## evaluated at each row of @var{x}.  @var{mu} is a 1-by-D vector, or an N-by-D
## matrix, in which case the density is evaluated for each row of @var{x} with
## the corresponding row of @var{mu}.  @var{mu} can also be a scalar value,
## which MVNPDF replicates to match the size of @var{x}.
##
## @code{@var{y} = mvnpdf (@var{x}, @var{mu}, @var{sigma})} returns the density
## of the multivariate normal distribution with mean @var{mu} and covariance
## @var{sigma}, evaluated at each row of @var{x}.  @var{sigma} is a D-by-D
## matrix, or an D-by-D-by-N array, in which case the density is evaluated for
## each row of @var{x} with the corresponding page of @var{sigma}, i.e.,
## @code{mvnpdf} computes @var{y(i)} using @var{x(i,:)} and @var{sigma(:,:,i)}.
## If the covariance matrix is diagonal, containing variances along the diagonal
## and zero covariances off the diagonal, @var{sigma} may also be specified as a
## 1-by-D matrix or a 1-by-D-by-N array, containing just the diagonal. Pass in
## the empty matrix for @var{mu} to use its default value when you want to only
## specify @var{sigma}.
##
## If @var{x} is a 1-by-D vector, @code{mvnpdf} replicates it to match the
## leading dimension of @var{mu} or the trailing dimension of @var{sigma}.
##
## @seealso{mvncdf, mvnrnd}
## @end deftypefn

function y = mvnpdf (x, mu, sigma)

  ## Check for valid number of input arguments
  if (nargin < 1)
    error ("mvnpdf: too few input arguments.");
  endif
  if (nargin < 3)
    sigma = [];
  endif

  ## Check for valid size of data
  [row, col] = size (x);
  if (col < 1)
    error ("mvnpdf: too few dimensions in X.");
  elseif (ndims (x) != 2)
    error ("mvnpdf: wrong dimensions in X.");
  endif

  ## Check for second input argument or assume zero mean
  if (nargin < 2 || isempty (mu))
    xc = x;                 # already centered
  elseif (numel (mu) == 1)
    xc = x - mu;            # mu is a scalar
  elseif (ndims (mu) == 2)
    [rm, cm] = size (mu);   # mu is a vector
    if (cm != col)
      error ("mvnpdf: columns in X and MU mismatch.");
    elseif (rm == row)
      xc = x - mu;
    elseif (rm == 1 || row == 1)
      xc = bsxfun (@minus, x, mu);
    else
      error ("mvnpdf: rows in X and MU mismatch.");
    endif
  else
    error ("mvnpdf: wrong size of MU.");
  endif
  [row, col] = size (xc);

  ## Check for third input argument or assume identity covariance
  if (nargin < 2 || isempty (sigma))
    ## already standardized
    if (col == 1 && row > 1)
      xRinv = xc';    # make row vector
      col == row;
    else
      xRinv = xc;
      col == row;
    endif
    lnSDS = 0;
  elseif (ndims (sigma) == 2)
    ## Single covariance matrix
    [rs, cs] = size (sigma);
    if (rs == 1 && cs > 1)
      rs = cs;        # sigma passed as a diagonal
      is_diag = true;
    else
      is_diag = false;
    endif
    if (col == 1 && row > 1 && rs == row)
      xc = xc';       # make row vector
      col = row;
    endif
    ## Check sigma for correct size
    if (rs != cs)
      error ("mvnpdf: bad covariance matrix.");
    elseif (rs != col)
      error ("mvnpdf: covariance matrix mismatch.");
    else
      if (is_diag)
        ## Check sigma for invalid values
        if (any (sigma <= 0))
          error ("mvnpdf: sigma diagonal contains negative or zero values.");
        endif
        R = sqrt (sigma);
        xRinv = bsxfun (@rdivide, xc, R);
        lnSDS = sum (log (R));
      else
        ## Check for valid covariance matrix
        [R, err] = cholcov (sigma, 0);
        if (err != 0)
          error ("mvnpdf: invalid covariance matrix.");
        endif
        xRinv = xc / R;
        lnSDS = sum (log (diag (R)));
      endif
    endif
  elseif (ndims (sigma) == 3)
    ## Multiple covariance matrices
    sd = size (sigma);
    if (sd(1) == 1 && sd(2) > 1)
      sd(1) = sd(2);  # sigma passed as a diagonal
      sigma = reshape (sigma, sd(2), sd(3))';
      is_diag = true;
    else
      is_diag = false;
    endif
    if (col == 1 && row > 1 && sd(1) == row)
      xc = xc';       # make row vector
      [row, col] = size (xc);
    endif
    ## If X and MU are row vectors, match them with covariance
    if (row == 1)
      row = sd(3);
      xc = repmat (xc, row, 1);
    endif
    ## Check sigma for correct size
    if (sd(1) != sd(2))
      error ("mvnpdf: bad multiple covariance matrix.");
    elseif (sd(1) != col || sd(2) != col)
      error ("mvnpdf: multiple covariance matrix mismatch.");
    elseif (sd(3) != row)
      error ("mvnpdf: multiple covariance pages mismatch.");
    else
      if (is_diag)
        ## Check sigma for invalid values
        if (any (any (sigma <= 0)))
          error ("mvnpdf: sigma diagonals contain negative or zero values.");
        endif
        R = sqrt (sigma);
        xRinv = xc ./ R;
        lnSDS = sum (log (R), 2);
      else
        ## Create arrays according to class type
        if (isa (x, "single") || isa (mu, "single") || isa (sigma, "single"))
          xRinv = zeros (row, col," single");
          lnSDS = zeros (row, 1, "single");
        else
          xRinv = zeros (row, col);
          lnSDS = zeros (row, 1);
        endif
        for i = 1:row
          ## Check for valid covariance matrices
          [R, err] = cholcov (sigma (:,:,i), 0);
          if (err != 0)
            error ("mvnpdf:invalid multiple covariance matrix.");
          endif
          xRinv(i,:) = xc(i,:) / R;
          lnSDS(i) = sum(log(diag(R)));
        endfor
      endif
    endif
  else
    error ("mvnpdf: wrong dimensions in covariance matrix.");
  endif

  ## Compute the PDF
  y = exp (-0.5 * sum (xRinv .^ 2, 2) - lnSDS - col * log (2 * pi) / 2);

endfunction

%!demo
%! mu = [1, -1];
%! sigma = [0.9, 0.4; 0.4, 0.3];
%! [X1, X2] = meshgrid (linspace (-1, 3, 25)', linspace (-3, 1, 25)');
%! x = [X1(:), X2(:)];
%! p = mvnpdf (x, mu, sigma);
%! surf (X1, X2, reshape (p, 25, 25));

## Input validation tests
%!error<mvnpdf: too few input arguments.> y = mvnpdf ();
%!error<mvnpdf: too few dimensions in X.> y = mvnpdf ([]);
%!error<mvnpdf: wrong dimensions in X.> y = mvnpdf (ones (3,3,3));
%!error<mvnpdf: columns in X and MU mismatch.> ...
%! y = mvnpdf (ones (10, 2), [4, 2, 3]);
%!error<mvnpdf: rows in X and MU mismatch.> ...
%! y = mvnpdf (ones (10, 2), [4, 2; 3, 2]);
%!error<mvnpdf: wrong size of MU.> ...
%! y = mvnpdf (ones (10, 2), ones (3, 3, 3));

## Output validation tests
%!shared x, mu, sigma
%! x = [1, 2, 5, 4, 6];
%! mu = [2, 0, -1, 1, 4];
%! sigma = [2, 2, 2, 2, 2];
%!assert (mvnpdf (x), 1.579343404440977e-20, 1e-30);
%!assert (mvnpdf (x, mu), 1.899325144348102e-14, 1e-25);
%!assert (mvnpdf (x, mu, sigma), 2.449062307156273e-09, 1e-20);
