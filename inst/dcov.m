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
## Distance covariance and correlation statistics.
##
## It returns distace correlation (@var{dCor}), 
## distance covariance (@var{dCov}), diatance variace on x (@var{dVarX}) and
## distance variance on y (@var{dVarY}).
##
## Reference: https://en.wikipedia.org/wiki/Distance_correlation
##
## @seealso{cov}
## @end deftypefn

function [dCov, dCor, dVarX, dVarY] = dcov (x,y,index=1.0)
  %x = abs(x - x.');
  %y = abs(y - y.');
  x = bsxfun (@minus, x, x.');
  y = bsxfun (@minus, y, y.');
  
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
        d = x .^ index;
        m = mean (d, 2);
        M = mean (d(:));
        %c = d - m - m.' + M;
        c = d - bsxfun (@plus, m, m.') + M;
endfunction
