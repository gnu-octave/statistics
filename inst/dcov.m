## Copyright (C) 2014 - Maria L. Rizzo and Gabor J. Szekely
## Copyright (C) 2014 Juan Pablo Carbajal
## Copyright (C) 2023 Andreas Bertsatos <abertsatos@biol.uoa.gr>
##
## This file is part of the statistics package for GNU Octave.
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program. If not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn  {statistics} {[@var{dCor}, @var{dCov}, @var{dVarX}, @var{dVarY}] =} dcov (@var{x}, @var{y})
##
## Distance correlation, covariance and correlation statistics.
##
## It returns the distance correlation (@var{dCor}) and the distance covariance
## (@var{dCov}) between @var{x} and @var{y}, the distance variance of @var{x}
## in (@var{dVarX}) and the distance variance of @var{y} in (@var{dVarY}).
##
## @var{x} and @var{y} must have the same number of observations (rows) but they
## can have different number of dimensions (columns).  Rows with missing values
## (@qcode{NaN}) in either @var{x} or @var{y} are omitted.
##
## The Brownian covariance is the same as the distance covariance:
##
## @tex
## $$ cov_W (X, Y) = dCov(X, Y) $$
##
## @end tex
## @ifnottex
## @math{cov_W (@var{x}, @var{y}) = dCov (@var{x}, @var{y})}
## @end ifnottex
##
## and thus Brownian correlation is the same as distance correlation.
##
## @seealso{corr, cov}
## @end deftypefn

function [dCor, dCov, dVarX, dVarY] = dcov (x, y)

  ## Validate input size
  if (size (x, 1) != size (y, 1))
    error ("dcov: Sample sizes (rows) in X and Y must agree.");
  endif

  ## Exclude missing values
  is_nan = any ([isnan(x) isnan(y)], 2);
  x(is_nan,:) = [];
  y(is_nan,:) = [];

  ## Calculate double centered distance
  A = pdist2 (x, x);
  A_col = mean (A, 1);
  A_row = mean (A, 2);
  Acbar = ones (size (A_row)) * A_col;
  Arbar = A_row * ones (size (A_col));
  A_bar = mean (A(:)) * ones (size (A));
  A = A - Acbar - Arbar + A_bar;

  B = pdist2 (y, y);
  B_col = mean (B, 1);
  B_row = mean (B, 2);
  Bcbar = ones (size (B_row)) * B_col;
  Brbar = B_row * ones (size (B_col));
  B_bar = mean (B(:)) * ones (size (B));
  B = B - Bcbar - Brbar + B_bar;

  ## Calculate distance covariance and variances
  dCov  = sqrt (mean (A(:) .* B(:)));
  dVarX = sqrt (mean (A(:) .^ 2));
  dVarY = sqrt (mean (B(:) .^ 2));

  ## Calculate distance correlation
  V = sqrt (dVarX .* dVarY);

  if V > 0
    dCor = dCov / V;
  else
    dCor = 0;
  end

endfunction

%!demo
%! base=@(x) (x- min(x))./(max(x)-min(x));
%! N = 5e2;
%! x = randn (N,1); x = base (x);
%! z = randn (N,1); z = base (z);
%! # Linear relations
%! cy = [1 0.55 0.3 0 -0.3 -0.55 -1];
%! ly = x .* cy;
%! ly(:,[1:3 5:end]) = base (ly(:,[1:3 5:end]));
%! # Correlated Gaussian
%! cz = 1 - abs (cy);
%! gy = base ( ly + cz.*z);
%! # Shapes
%! sx      = repmat (x,1,7);
%! sy      = zeros (size (ly));
%! v       = 2 * rand (size(x,1),2) - 1;
%! sx(:,1) = v(:,1); sy(:,1) = cos(2*pi*sx(:,1)) + 0.5*v(:,2).*exp(-sx(:,1).^2/0.5);
%! R       =@(d) [cosd(d) sind(d); -sind(d) cosd(d)];
%! tmp     = R(35) * v.';
%! sx(:,2) = tmp(1,:); sy(:,2) = tmp(2,:);
%! tmp     = R(45) * v.';
%! sx(:,3) = tmp(1,:); sy(:,3) = tmp(2,:);
%! sx(:,4) = v(:,1); sy(:,4) = sx(:,4).^2 + 0.5*v(:,2);
%! sx(:,5) = v(:,1); sy(:,5) = 3*sign(v(:,2)).*(sx(:,5)).^2  + v(:,2);
%! sx(:,6) = cos (2*pi*v(:,1)) + 0.5*(x-0.5);
%! sy(:,6) = sin (2*pi*v(:,1)) + 0.5*(z-0.5);
%! sx(:,7) = x + sign(v(:,1)); sy(:,7) = z + sign(v(:,2));
%! sy      = base (sy);
%! sx      = base (sx);
%! # scaled shape
%! sc  = 1/3;
%! ssy = (sy-0.5) * sc + 0.5;
%! n = size (ly,2);
%! ym = 1.2;
%! xm = 0.5;
%! fmt={'horizontalalignment','center'};
%! ff = "% .2f";
%! figure (1)
%! for i=1:n
%!   subplot(4,n,i);
%!   plot (x, gy(:,i), '.b');
%!   axis tight
%!   axis off
%!   text (xm,ym,sprintf (ff, dcov (x,gy(:,i))),fmt{:})
%!
%!   subplot(4,n,i+n);
%!   plot (x, ly(:,i), '.b');
%!   axis tight
%!   axis off
%!   text (xm,ym,sprintf (ff, dcov (x,ly(:,i))),fmt{:})
%!
%!   subplot(4,n,i+2*n);
%!   plot (sx(:,i), sy(:,i), '.b');
%!   axis tight
%!   axis off
%!   text (xm,ym,sprintf (ff, dcov (sx(:,i),sy(:,i))),fmt{:})
%!   v = axis ();
%!
%!   subplot(4,n,i+3*n);
%!   plot (sx(:,i), ssy(:,i), '.b');
%!   axis (v)
%!   axis off
%!   text (xm,ym,sprintf (ff, dcov (sx(:,i),ssy(:,i))),fmt{:})
%! endfor

%!error dcov (randn (30, 5), randn (25,5))
