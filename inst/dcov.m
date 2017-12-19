## Copyright (C) 2014 - Maria L. Rizzo and Gabor J. Szekely
## Copyright (C) 2014 Juan Pablo Carbajal
## This work is derived from the R energy package. It was adapted 
## for Octave by Juan Pablo Carbajal.
## 
## This progrm is free software; you can redistribute it and/or modify
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

## Author: Juan Pablo Carbajal <ajuanpi+dev@gmail.com>

## -*- texinfo -*-
## @deftypefn {Function File} {[@var{dCor}, @var{dCov}, @var{dVarX}, @var{dVarY}] =} dcov (@var{x}, @var{y}, @var{index}=1)
## Distance correlation, covariance and correlation statistics.
##
## It returns distace correlation (@var{dCor}), 
## distance covariance (@var{dCov}), diatance variace on x (@var{dVarX}) and
## distance variance on y (@var{dVarY}).
##
## Reference: https://en.wikipedia.org/wiki/Distance_correlation
##
## @seealso{corr, cov}
## @end deftypefn

function [dCor, dCov, dVarX, dVarY] = dcov (x,y,index=1.0)
  %x = abs(x - x.');
  %y = abs(y - y.');
  x = abs (bsxfun (@minus, x, x.'));
  y = abs (bsxfun (@minus, y, y.'));
  
  [n nc] = size (x);
  [m mc] = size (y);
  if (n != m) 
    error ("Octave:invalid-input-arg", "Sample sizes must agree.");
  endif
  
  if  any (isnan (x) | isnan (y))
      error ("Octave:invalid-input-arg","Data contains missing or infinite values.");
  endif
  
  if index < 0 || index > 2
    warning ("Octave:invalid-input-arg","index must be in [0,2), using default index=1");
    index = 1.0;
  endif

  A = Akl (x, index);
  B = Akl (y, index);

  dCov  = sqrt (mean (A(:) .* B(:)));
  dVarX = sqrt (mean (A(:).^2) );
  dVarY = sqrt (mean (B(:).^2) );
  V     = sqrt (dVarX .* dVarY);

  if V > 0
    dCor = dCov / V;
  else 
    dCor = 0;
  end

endfunction

function c = Akl (x, index)
# Double centered distance
        d = x .^ index;
        rm = mean (d, 2); # row mean
        gm = mean (d(:)); # grand mean
        c  = d - bsxfun (@plus, rm, rm.') + gm;
endfunction

%!demo
%! x = randn (1e3,1); x = (x- mean(x))./(max(x)-mean(x));
%! z = randn (1e3,1);
%! # Linear relations
%! ly = x.*[1 0.55 0.3 0 -0.3 -0.55 -1];
%! # Correlated Gaussian
%! gy =  ly + [0 0.45 0.7 1 0.7 0.45 0].*z;
%! # Shapes
%! sx = repmat (x,7);
%! sy = zeros (size (ly));
%! sy(:,1) = cos(pi*x) + z;
%! R =@(d) [cosd(d) sind(d); -sind(d) cosd(d)];
%! tmp     = R(35) * [x.';z.'];
%! sx(:,2) = tmp(1,:) - mean (tmp(1,:)); sy(:,2) = tmp(2,:);
%! tmp     = R(45) * [x.';z.'];
%! sx(:,3) = tmp(1,:) - mean (tmp(1,:)); sy(:,3) = tmp(2,:);
%! sy(:,4) = x.^2 + 2*abs(z);
%! sy(:,5) = x.^2 .* z;
%! sx(:,6) = cos (pi*x);
%! sy(:,6) = sin (pi*x);
%! tf      = [(x>=0  & z>=0) (x<0 & z>=0) (x>=0 & z<0) (x<0 & z<0)];
%! sy(:,7) = gy(:,4) .* sum((2*tf-1),2);
%!
%! n = size (ly,2);
%! ym = max ([gy(:) ly(:) sy(:)]);
%! xm = 0;
%! fmt={'horizontalalignment','center'};
%! figure (1)
%! for i=1:n
%!  subplot(3,n,i);
%!  plot (x, gy(:,i), '.b');
%!  axis equal
%!  axis off
%!  text (xm,ym(1)*1.5,sprintf ("%.1f", dcov (x,gy(:,i))),fmt{:})
%!
%!  subplot(3,n,i+n);
%!  plot (x, ly(:,i), '.b');
%!  axis equal
%!  axis off
%!  text (xm,ym(2)*1.5,sprintf ("%.1f", dcov (x,ly(:,i))),fmt{:})
%!
%!  subplot(3,n,i+2*n);
%!  plot (x, sy(:,i), '.b');
%!  axis equal
%!  axis off
%!  text (xm,ym(3)*1.5,sprintf ("%.1f", dcov (x,sy(:,i))),fmt{:})
%! endfor


