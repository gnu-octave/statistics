## Copyright (C) 2014 - Nir Krakauer
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

## Author: Nir Krakauer <nkrakauer@ccny.cuny.com>

## -*- texinfo -*-
## @deftypefn {Function File} {@var{y} =} randsample (@var{v}, @var{k}, @var{replacement}=false [, @var{w}])
## Elements sampled from a vector.
##
## Returns @var{k} random elements from a vector @var{v} with @var{n} elements, sampled without or with @var{replacement}.
##
## If @var{v} is a scalar, samples from 1:@var{v}.
##
## If sampling with replacement, can specify a weight vector @var{w} of the same size as @var{v} such that the probablility of each element being sampled is proportional to @var{w}.
##
## Randomization is performed using rand().
##
## @seealso{randperm}
## @end deftypefn

function y = randsample(v,k,replacement=false,w=[])

  if isscalar (v)
    n = v;
    vector_v = false;
  elseif isvector (v)
    n = numel (v);
    vector_v = true;
  else
    error ('Octave:invalid-input-arg', 'The input v must be a vector or positive integer.');
  endif
   
  if k < 0 || k > n
    error ('Octave:invalid-input-arg', 'The input k must be an integer between 0 and n.');
  endif


  if replacement #sample with replacement
    if isempty (w) #all elements are equally likely to be sampled
      y = round (n * rand(1, k) + 0.5);
    else
      w = w / sum(w);
      w = [0 cumsum(w(:))'];
      y = arrayfun(@(x) find(w <= x, 1, "last"), rand (1, k)); #distribute k uniform random deviates based on the given weighting
    endif
  else #sample without replacement
    y = randperm (n, k);
  endif
  
  if vector_v
    y = v(y);
  endif
  
endfunction 
 

%!test
%! n = 20;
%! k = 5;
%! x = randsample(n, k);
%! assert (size(x), [1 k]);
%! x = randsample(n, k, true);
%! assert (size(x), [1 k]);
%! x = randsample(n, k, false);
%! assert (size(x), [1 k]);
%! x = randsample(n, k, true, ones(n, 1));
%! assert (size(x), [1 k]);
%! x = randsample(1:n, k);
%! assert (size(x), [1 k]);
%! x = randsample(1:n, k, true);
%! assert (size(x), [1 k]);
%! x = randsample(1:n, k, false);
%! assert (size(x), [1 k]);
%! x = randsample(1:n, k, true, ones(n, 1));
%! assert (size(x), [1 k]);
%! x = randsample((1:n)', k);
%! assert (size(x), [k 1]);
%! x = randsample((1:n)', k, true);
%! assert (size(x), [k 1]);
%! x = randsample((1:n)', k, false);
%! assert (size(x), [k 1]);
%! x = randsample((1:n)', k, true, ones(n, 1));
%! assert (size(x), [k 1]);
