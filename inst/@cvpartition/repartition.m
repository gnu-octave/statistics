## Copyright (C) 2014 Nir Krakauer
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
## along with this program; If not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn{Function File}{@var{Cnew} =} repartition (@var{C})
## Return a new cvpartition object.
##
## @var{C} should be a cvpartition object. @var{Cnew} will use the same partition_type as @var{C} but redo any randomization performed (currently, only the HoldOut type uses randomization).
##
## @seealso{cvpartition}
## @end deftypefn

## Author: Nir Krakauer

function Cnew = repartition (C)

  if (nargin < 1 || nargin > 2)
    print_usage ();
  endif

  Cnew = C;

  switch C.Type
    case 'kfold'
    case 'given'
    case 'holdout' #currently, only the HoldOut method uses randomization
      n = C.NumObservations;
      k = C.TestSize;
      n_classes = C.n_classes;
      if k < 1
        f = k; #target fraction to sample
        k = round (k * n); #number of samples
      else
        f = k / n;
      endif
      inds = zeros (n, 1, "logical");
      if n_classes == 1
        inds(randsample(n, k)) = true; #indices for test set
      else #sample from each class
        j = C.classes; #integer class labels
        n_per_class = accumarray (j, 1);
        n_classes = numel (n_per_class);
        k_check = 0;
        for i = 1:n_classes
          ki = round(f*n_per_class(i));
          inds(find(j == i)(randsample(n_per_class(i), ki))) = true;
          k_check += ki;
        endfor
        if k_check < k #add random elements to test set to make it k
          inds(find(!inds)(randsample(n - k_check, k - k_check))) = true;
        elseif k_check > k #remove random elements from test set
          inds(find(inds)(randsample(k_check, k_check - k))) = false;
        endif
      endif
      Cnew.inds = inds;
    case 'leaveout'
    case 'resubstitution'
  endswitch
