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
## @deftypefn {Function File} {@var{results} =} crossval (@var{f}, @var{X}, @var{y}[, @var{params}])
## Perform cross validation on given data.
##
## @var{f} should be a function that takes 4 inputs @var{xtrain}, @var{ytrain},
## @var{xtest}, @var{ytest}, fits a model based on @var{xtrain}, @var{ytrain},
## applies the fitted model to @var{xtest}, and returns a goodness of fit
## measure based on comparing the predicted and actual @var{ytest}.
## @code{crossval} returns an array containing the values returned by @var{f}
## for every cross-validation fold or resampling applied to the given data.
##
## @var{X} should be an @var{n} by @var{m} matrix of predictor values
##
## @var{y} should be an @var{n} by @var{1} vector of predicand values
##
## @var{params} may include parameter-value pairs as follows:
##
## @table @asis
## @item @qcode{"KFold"}
## Divide set into @var{k} equal-size subsets, using each one successively
## for validation.
##
## @item @qcode{"HoldOut"}
## Divide set into two subsets, training and validation. If the value
## @var{k} is a fraction, that is the fraction of values put in the
## validation subset (by default @var{k}=0.1); if it is a positive integer,
## that is the number of values in the validation subset.
##
## @item @qcode{"LeaveOut"}
## Leave-one-out partition (each element is placed in its own subset).
## The value is ignored.
##
## @item @qcode{"Partition"}
## The value should be a @var{cvpartition} object.
##
## @item @qcode{"Given"}
## The value should be an @var{n} by @var{1} vector specifying in which
## partition to put each element.
##
## @item @qcode{"stratify"}
## The value should be an @var{n} by @var{1} vector containing class
## designations for the elements, in which case the @qcode{"KFold"} and
## @qcode{"HoldOut"} partitionings attempt to ensure each partition
## represents the classes proportionately.
##
## @item @qcode{"mcreps"}
## The value should be a positive integer specifying the number of times
## to resample based on different partitionings. Currently only works with
## the partition type @qcode{"HoldOut"}.
##
## @end table
##
## Only one of @qcode{"KFold"}, @qcode{"HoldOut"}, @qcode{"LeaveOut"},
## @qcode{"Given"}, @qcode{"Partition"} should be specified. If none is
## specified, the default is @qcode{"KFold"} with @var{k} = 10.
##
## @seealso{cvpartition}
## @end deftypefn

## Author: Nir Krakauer

function results = crossval (f, X, y, varargin)

  [n m] = size (X);
  
  if numel(y) != n
    error('X, y sizes incompatible')
  endif
  
  #extract optional parameter-value argument pairs
  if numel(varargin) > 1
    vargs = varargin;
    nargs = numel (vargs);
    values = vargs(2:2:nargs);
    names = vargs(1:2:nargs)(1:numel(values));
    validnames = {'KFold', 'HoldOut', 'LeaveOut', 'Partition', 'Given', 'stratify', 'mcreps'};    
    for i = 1:numel(names)
       names(i) = validatestring (names(i){:}, validnames);
    end
    for i = 1:numel(validnames)
      name = validnames(i){:};
      name_pos = strmatch (name, names);
      if !isempty(name_pos)
        eval([name ' = values(name_pos){:};'])
      endif
    endfor
  endif
    
  #construct CV partition
  if exist ("Partition", "var")
    P = Partition;
  elseif exist ("Given", "var")
    P = cvpartition (Given, "Given");
  elseif exist ("KFold", "var")
    if !exist ("stratify", "var")
      stratify = n;
    endif
    P = cvpartition (stratify, "KFold", KFold);
  elseif exist ("HoldOut", "var")
    if !exist ("stratify", "var")
      stratify = n;
    endif
    P = cvpartition (stratify, "HoldOut", HoldOut);
    if !exist ("mcreps", "var") || isempty (mcreps)
      mcreps = 1;
    endif
  elseif exist ("LeaveOut", "var")
    P = cvpartition (n, "LeaveOut");
  else #KFold
    if !exist ("stratify", "var")
      stratify = n;
    endif
    P = cvpartition (stratify, "KFold");
  endif
  
  nr = get(P, "NumTestSets"); #number of test sets to do cross validation on
  nreps = 1;
  if strcmp(get(P, "Type"), 'holdout') && exist("mcreps", "var") && mcreps > 1
    nreps = mcreps;
  endif
  results = nan (nreps, nr);
  for rep = 1:nreps
    if rep > 1
      P = repartition (P);
    endif
    for i = 1:nr
      inds_train = training (P, i);  
      inds_test = test (P, i);
      result = f (X(inds_train, :), y(inds_train), X(inds_test, :), y(inds_test));
      results(rep, i) = result;
    endfor
  endfor

endfunction

%!test
%! load fisheriris.txt
%! y = fisheriris(:, 2);
%! X = [ones(size(y)) fisheriris(:, 3:5)];
%! f = @(X1, y1, X2, y2) meansq (y2 - X2*regress(y1, X1));
%! results0 = crossval (f, X, y);
%! results1 = crossval (f, X, y, 'KFold', 10);
%! folds = 5;
%! results2 = crossval (f, X, y, 'KFold', folds);
%! results3 = crossval (f, X, y, 'Partition', cvpartition (numel (y), 'KFold', folds));
%! results4 = crossval (f, X, y, 'LeaveOut', 1);
%! mcreps = 2; n_holdout = 20;
%! results5 = crossval (f, X, y, 'HoldOut', n_holdout, 'mcreps', mcreps);
%!
%! ## ensure equal representation of iris species in the training set -- tends
%! ## to slightly reduce cross-validation mean square error 
%! results6 = crossval (f, X, y, 'KFold', 5, 'stratify', fisheriris(:, 1));
%!
%! assert (results0, results1);
%! assert (results2, results3);
%! assert (size(results4), [1 numel(y)]);
%! assert (mean(results4), 4.5304, 1E-4);
%! assert (size(results5), [mcreps 1]);

