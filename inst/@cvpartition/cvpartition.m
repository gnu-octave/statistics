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
## @deftypefn{Function File}{@var{C} =} cvpartition (@var{X}, [@var{partition_type}, [@var{k}]])
## Create a partition object for cross validation.
##
## @var{X} may be a positive integer, interpreted as the number of values @var{n} to partition, or a vector of length @var{n} containing class designations for the elements, in which case the partitioning types @var{KFold} and @var{HoldOut} attempt to ensure each partition represents the classes proportionately.
##
## @var{partition_type} must be one of the following:
##
## @table @asis
## @item @samp{KFold}
## Divide set into @var{k} equal-size subsets (this is the default, with @var{k}=10).
## @item @samp{HoldOut}
## Divide set into two subsets, "training" and "validation". If @var{k} is a fraction, that is the fraction of values put in the validation subset; if it is a positive integer, that is the number of values in the validation subset (by default @var{k}=0.1).
## @item @samp{LeaveOut}
## Leave-one-out partition (each element is placed in its own subset).
## @item @samp{resubstitution}
## Training and validation subsets that both contain all the original elements.
## @item @samp{Given}
## Subset indices are as given in @var{X}.
## @end table
##
## The following fields are defined for the @samp{cvpartition} class:
## 
## @table @asis
## @item @samp{classes}
## Class designations for the elements.
## @item @samp{inds}
## Subset indices for the elements.
## @item @samp{n_classes}
## Number of different classes.
## @item @samp{NumObservations}
## @var{n}, number of elements in data set.
## @item @samp{NumTestSets}
## Number of testing subsets.
## @item @samp{TestSize}
## Number of elements in (each) testing subset.
## @item @samp{TrainSize}
## Number of elements in (each) training subset.
## @item @samp{Type}
## Partition type.
## @end table
##
## @seealso{crossval}
## @end deftypefn

## Author: Nir Krakauer

function C = cvpartition (X, partition_type = 'KFold', k = [])

  if (nargin < 1 || nargin > 3 || !isvector(X))
    print_usage ();
  endif
  
  if isscalar (X)
    n = X;
    n_classes = 1;
  else
    n = numel (X);
  endif

  switch tolower(partition_type)
    case {'kfold' 'holdout' 'leaveout' 'resubstitution' 'given'}
    otherwise
      warning ('unrecognized type, using KFold')
      partition_type = 'KFold';
  endswitch

  switch tolower(partition_type)
    case {'kfold' 'holdout' 'given'}  
      if !isscalar (X)  
        [y, ~, j] = unique (X(:));
        n_per_class = accumarray (j, 1);
        n_classes = numel (n_per_class);
      endif
  endswitch

  C = struct ("classes", [], "inds", [], "n_classes", [], "NumObservations", [], "NumTestSets", [], "TestSize", [], "TrainSize", [], "Type", []);
  #The non-Matlab fields classes, inds, n_classes are only useful for some methods

  switch tolower(partition_type)
    case 'kfold'
      if isempty (k)
        k = 10;
      endif
      if n_classes == 1
        inds = floor((0:(n-1))' * (k / n)) + 1;
      else
        inds = nan(n, 1);
        for i = 1:n_classes
          if mod (i, 2) #alternate ordering over classes so that the subsets are more nearly the same size
            inds(j == i) = floor((0:(n_per_class(i)-1))' * (k / n_per_class(i))) + 1;
          else
            inds(j == i) = floor(((n_per_class(i)-1):-1:0)' * (k / n_per_class(i))) + 1;
          endif
        endfor
      endif
      C.inds = inds;
      C.NumTestSets = k;
      [~, ~, jj] = unique (inds);
      n_per_subset = accumarray (jj, 1);     
      C.TrainSize = n - n_per_subset;
      C.TestSize = n_per_subset;
    case 'given'
      C.inds = j;
      C.NumTestSets = n_classes;
      C.TrainSize = n - n_per_class;
      C.TestSize = n_per_class;          
    case 'holdout'
      if isempty (k)
        k = 0.1;
      endif
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
        C.classes = j;
      endif        
      C.n_classes = n_classes;
      C.TrainSize = n - k;
      C.TestSize = k;
      C.NumTestSets = 1;
      C.inds = inds;    
    case 'leaveout'
      C.TrainSize = ones (n, 1);
      C.TestSize = (n-1)  * ones (n, 1);
      C.NumTestSets = n;
    case 'resubstitution'
      C.TrainSize = C.TestSize = n;
      C.NumTestSets = 1;
  endswitch
  
  C.NumObservations = n;
  C.Type = tolower (partition_type);

  C = class (C, "cvpartition");

endfunction


%!demo
%! # Partition with Fisher iris dataset (n = 150)
%! # Stratified by species
%! load fisheriris.txt
%! y = fisheriris(:, 1);
%! # 10-fold cross-validation partition
%! c = cvpartition (y, 'KFold', 10)
%! # leave-10-out partition
%! c1 = cvpartition (y, 'HoldOut', 10)
%! idx1 = test (c, 2);
%! idx2 = training (c, 2);
%! # another leave-10-out partition
%! c2 = repartition (c1)
#plot(struct(c).inds, '*')

