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
