## Copyright (C) 2014 - Juan Pablo Carbajal
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
## @deftypefn {Function File} {@var{z} =} qrandn (@var{q}, @var{r},@var{c})
## @deftypefnx {Function File} {@var{z} =} qrandn (@var{q}, [@var{r},@var{c}])
## Returns random deviates drawn from a q-Gaussian distribution.
##
## Parameter @var{q} charcterizes the q-Gaussian distribution.
## The result has the size indicated by @var{s}.
##
## Reference:
## W. Thistleton, J. A. Marsh, K. Nelson, C. Tsallis (2006)
## "Generalized Box-Muller method for generating q-Gaussian random deviates"
## arXiv:cond-mat/0605570 http://arxiv.org/abs/cond-mat/0605570
##
## @seealso{rand, randn}
## @end deftypefn

function z = qrandn(q,R,C=[]) 
  if !isscalar (q)
    error ('Octave:invalid-input-arg', 'The parameter q must be a scalar.')
  endif
  
  # Check that q < 3 
  if q > 3
    error ('Octave:invalid-input-arg', 'The parameter q must be lower than 3.');
  endif

  if numel (R) > 1
  S = R;
  elseif numel (R) ==1 && isempty (C)
  S = [R,1];
  elseif numel (R) ==1 && !isempty (C)
  S = [R,C];
  endif

  # Calaulate the q to be used on the q-log 
  qGen = (1 + q) / (3 - q); 

  # Initialize the output vector 
  z = sqrt (-2 * log_q (rand (S),qGen)) .* sin (2*pi*rand (S)); 

endfunction 
 
function a = log_q (x,q) 
  # 
  # Returns the q-log of x, using q 
  # 
  dq = 1 - q;
  # Check to see if q = 1 (to double precision) 
  if abs (dq) < 10*eps 
    # If q is 1, use the usual natural logarithm 
    a = log (x); 
  else 
    # If q differs from 1, use the definition of the q-log 
    a = ( x .^ dq - 1 ) ./ dq; 
  endif
   
endfunction

%!demo
%! z = qrandn (-5, 5e6);
%! [c x] = hist (z,linspace(-1.5,1.5,200),1);
%! figure(1)
%! plot(x,c,"r."); axis tight; axis([-1.5,1.5]);
%!
%! z = qrandn (-0.14286, 5e6);
%! [c x] = hist (z,linspace(-2,2,200),1);
%! figure(2)
%! plot(x,c,"r."); axis tight; axis([-2,2]);
%!
%! z = qrandn (2.75, 5e6);
%! [c x] = hist (z,linspace(-1e3,1e3,1e3),1);
%! figure(3)
%! semilogy(x,c,"r."); axis tight; axis([-100,100]);
%!
%! # ---------
%! # Figures from the reference paper.
