## Copyright (C) 2014 Nir Krakauer
## Copyright (C) 2025 Yassin Achengli
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
## along with this program; If not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn  {statistics} {@var{results} =} crossval (@var{f}, @var{X}, @var{y})
## @deftypefnx {statistics} {@var{results} =} crossval (@var{f}, @var{X}, @var{y}, @var{name}, @var{value})
##
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
## Optional arguments may include name-value pairs as follows:
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

function results = crossval (f, X, varargin)
  if (nargin < 2)
    print_usage
  endif

  n = size (X, 1);

  ## extract optional parameter-value argument pairs
  if (ischar (f) && nargin < 5)
    error ("categorical evaluation needs X, y inputs. Only X given")
  else
    y = cell2mat (varargin(1));
    varargin = varargin(2:end);
  endif

  if (ischar (f) && ! strcmp (f, "mse") && ! strcmp (f, "mcr"))
    error ("Bad criterion. Must be MSE or MCR")
  endif

  if (numel (varargin) > 1)
    vargs = varargin;
    nargs = numel (vargs);
    values = vargs(2:2:nargs);
    names = vargs(1:2:nargs)(1:numel(values));
    validnames = {"KFold", "HoldOut", "LeaveOut", "Partition", ...
    "Given", "stratify", "mcreps", "Predfun"};
    for i = 1:numel (names)
      names(i) = validatestring (names(i){:}, validnames);
    end
    for i = 1:numel(validnames)
      name = validnames(i){:};
      name_pos = find (strcmp (name, names));
      if (! isempty (name_pos))
        eval ([name " = values(name_pos){:};"])
      endif
    endfor
  endif

  if (ischar (f) && ! (exist ("Predfun", "var")))
    error ("crossval for error validation needs defined Predfun parameter")
  endif

  ## construct CV partition
  if exist ("Partition", "var")
    P = Partition;
  elseif exist ("Given", "var")
    P = cvpartition (Given, "Given");
  elseif exist ("KFold", "var")
    if (! exist ("stratify", "var"))
      stratify = n;
    endif
    P = cvpartition (stratify, "KFold", KFold);
  elseif (exist ("HoldOut", "var"))
    if (! exist ("stratify", "var"))
      stratify = n;
    endif
    P = cvpartition (stratify, "HoldOut", HoldOut);
    if (! exist ("mcreps", "var") || isempty (mcreps))
      mcreps = 1;
    endif
  elseif (exist ("LeaveOut", "var"))
    P = cvpartition (n, "LeaveOut");
  else #KFold
    if (! exist ("stratify", "var"))
      stratify = n;
    endif
    P = cvpartition (stratify, "KFold");
  endif

  nr = P.NumTestSets;
  nreps = 1;
  if (strcmp (P.Type, "holdout") && exist ("mcreps", "var") &&  mcreps > 1)
    nreps = mcreps;
  endif

  if (ischar (f))
    results = nan (nreps, nr);
    for rep = 1:nreps
      if (rep > 1)
        P = repartition (P);
      endif
      for idx = 1:nr
        idx_train = training (P, idx);
        idx_test = test (P, idx);
        y_fit = Predfun (X(idx_train, :), y(idx_train), X(idx_test, :));
        if (strcmp (f, "mse"))
          N = numel(y_fit) - 1;
          if (N < 1)
            N = 1;
          endif
          err = sum ((y_fit - y(idx_test)).^2) / (N - 1);
          results(rep, idx) = err;
        else # MCR 
          err = numel (y_fit(y_fit == y(idx_test))) / numel (y_fit);
          results(rep, idx) = err;
        endif
      endfor
    endfor
    results = mean (mean (results));
  else
    # Model execution
    for rep = 1:nreps
      if (rep > 1)
        P = repartition (P);
      endif
      for idx = 1:nr
        idx_train = training (P, idx);
        idx_test = test (P, idx);
        result = f (X(idx_train, :), X(idx_test, :));

        results(rep, idx) = result;
      endfor
    endfor
  endif
endfunction

%!test
%! load fisheriris
%! y = meas(:, 1);
%! X = [ones(size(y)) meas(:, 2:4)];
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
%! results6 = crossval (f, X, y, 'KFold', 5, 'stratify', grp2idx(species));
%!
%!# assert (results0, results1, 2e-15);
%!# assert (results2, results3, 5e-17);
%!# assert (size(results4), [1 numel(y)]);
%!# assert (mean(results4), 0.1018, 1e-4);
%!# assert (size(results5), [mcreps 1]);


