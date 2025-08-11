## Copyright (C) 2014 Nir Krakauer
## Copyright (C) 2025 Yassin Achengli
## Copyright (C) 2025 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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

function results = crossval (f, varargin)

  ## Parse optional Name-Value paired arguments
  optNames = {'Holdout', 'KFold', 'Leaveout', 'MCReps', ...
              'Partition', 'Stratify', 'Predfun'};
  dfValues = {[], [], [], 1, [], [], []};
  [Holdout, KFold, Leaveout, MCReps, Partition, Stratify, Predfun, args] = ...
                                   pairedArgs (optNames, dfValues, varargin(:));

  ## Check first input argument
  if (ischar (f))
    ## Check for valid criterion
    if (! ismember (f, {'mse',',mcr'}))
      error ("crossval: criterion must be 'mse' or 'mcr'.");
    endif
    ## Check for valid prediction function handle
    if (is_function_handle (Predfun))
      ## Check for valid prediction function handle
      XTrain = [1, 2; 1, 2; 1, 2; 1, 2; 1, 2];
      yTrain = [1; 2; 1; 2; 1];
      XTest = [1, 2; 1, 2];
      try
        yFit = Predfun (XTrain, yTrain, XTest);
      catch
        error ("crossval: bad prediction function handle for error evaluation.");
      end_try_catch
      if (! iscolumn (yFit) || numel (yFit) != 2)
        error (strcat ("crossval: prediction function handle must return", ...
                       " a column vector with the same rows as XTest."));
      endif
    else
      error (strcat ("crossval: prediction function handle", ...
                     " is required for error evaluation."));
    endif
    ## At least two additional input arguments (X and y) are required
    nargs = numel (args);
    if (nargs < 2)
      error ("crossval: X and Y are required for error evaluation.");
    endif
    y = args{end};
    if (nargs == 2)
      X = args{1};        # numeric matrix with single data variable
      n = size (X, 1);
    else
      X = args(1:end-1);  # cell array with multiple data variables
      n = size (args{1}, 1);
    endif

  elseif (is_function_handle (f))
    ## Check for valid function handle
    XTrain = [1, 2; 1, 2; 1, 2; 1, 2; 1, 2];
    XTest = [1, 2; 1, 2];
    try
      value = f (XTrain, XTest);
    catch
      error ("crossval: bad function handle to cross-validate.");
    end_try_catch
    if (! isrow (value))
      error ("crossval: function handle must return a scalar or a row vector.");
    endif
    ## At least one additional input argument (X) is required
    nargs = numel (args);
    if (nargs < 1)
      error ("crossval: X is required for values evaluation.");
    endif
    if (nargs == 1)
      X = args{1};        # numeric matrix with single data variable
      n = size (X, 1);
    else
      X = args(1:end-1);  # cell array with multiple data variables
      n = size (args{1}, 1);
    endif
  else
    error ("crossval: invalid first input argument.");
  endif

  ## Input validation for valid values are handled by the cvpartition class.
  ## Check for single paired argument for CV partition.
  vcpa = sum (isempty (Holdout), isempty (KFold), ...
              isempty (Leaveout), isempty (Partition));
  if (vcpa > 1)
    error (strcat ("crossval: you can only set one", ...
                   " cvpartition type in paired arguments."));
  endif
  ## Check for Partition and Stratify paired arguments
  if (! isempty (Partition) && ! isempty (Stratify))
    error ("crossval: you cannot specify both 'Partition' and 'Stratify'.");
  endif

  ## Construct the CV partition
  if (! isempty (Partition))
    P = Partition;
    if (P.IsCustom || ismember (P.Type, {'resubstitution', 'leaveout'}))
      MCReps = 1;
    endif
  elseif (! isempty (Leaveout))
    P = cvpartition (n, "LeaveOut");
    MCReps = 1;
  elseif (! isempty (Holdout))
    if (isempty (Stratify))
      P = cvpartition (n, "HoldOut", Holdout);
    else
      P = cvpartition (Stratify, "HoldOut", Holdout);
    endif
  elseif (! isempty (KFold))
    if (isempty (Stratify))
      P = cvpartition (n, "KFold", KFold);
    else
      P = cvpartition (Stratify, "KFold", KFold);
    endif
  else # KFold by default
    if (isempty (Stratify))
      P = cvpartition (n, "KFold");
    else
      P = cvpartition (Stratify, "KFold");
    endif
  endif

  ## FIX ME: this requires further work to handle multiple data variables
  ## Currently, it works only for single data variable, where X is numeric.
  if (ischar (f)) # error evaluation
    results = nan (MCReps, nr);
    for rep = 1:MCReps
      if (rep > 1)
        P = repartition (P);
      endif
      for idx = 1:nr
        idx_train = training (P, idx);
        idx_test = test (P, idx);
        y_fit = Predfun (X(idx_train, :), y(idx_train), X(idx_test, :));
        if (strcmp (f, "mse"))
          N = numel (y_fit) - 1;
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

  else  # model execution
    for rep = 1:MCReps
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


