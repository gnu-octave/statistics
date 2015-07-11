## Copyright (C) 2014 - Nir Krakauer
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

## Author: Nir Krakauer <nkrakauer@ccny.cuny.com>

## -*- texinfo -*-
## @deftypefn {Function File} {@var{y} =} randsample (@var{v}, @var{k}, @var{replacement}=false [, @var{w}])
## Elements sampled from a vector.
##
## Returns @var{k} random elements from a vector @var{v} with @var{n} elements, sampled without or with @var{replacement}.
##
## If @var{v} is a scalar, samples from 1:@var{v}.
##
## If a weight vector @var{w} of the same size as @var{v} is specified, the probablility of each element being sampled is proportional to @var{w}.  Unlike Matlab's function of the same name, this can be done for sampling with or without replacement.
##
## Randomization is performed using rand().
##
## @seealso{randperm}
## @end deftypefn

function y = randsample(v,k,replacement=false,w=[])

  if (isscalar (v) && isreal (v))
    n = v;
    vector_v = false;
  elseif (isvector (v))
    n = numel (v);
    vector_v = true;
  else
    error ('Octave:invalid-input-arg', 'randsample: The input v must be a vector or positive integer.');
  endif
   
  if k < 0 || ( k > n && !replacement )
    error ('Octave:invalid-input-arg', 'randsample: The input k must be a non-negative integer. Sampling without replacement needs k <= n.');
  endif

  if (all (length (w) != [0, n]))
    error ('Octave:invalid-input-arg', 'randsample: the size w (%d) must match the first argument (%d)', length(w), n);
  endif


  if (replacement)               # sample with replacement
    if (isempty (w))             # all elements are equally likely to be sampled
      y = round (n * rand(1, k) + 0.5);
    else
      y = weighted_replacement (k, w);
    endif
   else                           # sample without replacement
     if (isempty (w))             # all elements are equally likely to be sampled
       y = randperm (n, k);
     else                         # use "accept-reject"-like sampling
       y = weighted_replacement (k, w);
       while (1)
         [yy, idx] = sort (y);    # Note: sort keeps order of equal elements.
         Idup = [false, (diff (yy)==0)];
         if !any (Idup)
           break
         else
           Idup(idx) = Idup;      # find duplicates in original vector
           w(y) = 0;              # don't permit resampling
                      # remove duplicates, then sample again
           y = [y(~Idup), (weighted_replacement (sum (Idup), w))];
         endif
       endwhile
     endif
  endif
  
  if vector_v
    y = v(y);
  endif
  
endfunction 
 
function y = weighted_replacement (k, w)
  w = w / sum(w);
  w = [0 cumsum(w(:))'];
        # distribute k uniform random deviates based on the given weighting
  y = arrayfun (@(x) find (w <= x, 1, "last"), rand (1, k));
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
%! n = 10;
%! k = 100;
%! x = randsample(n, k, true, 1:n);
%! assert (size(x), [1 k]);
%! x = randsample((1:n)', k, true);
%! assert (size(x), [k 1]);
%! x = randsample(k, k, false, 1:k);
%! assert (size(x), [1 k]);
