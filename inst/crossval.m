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
## The value is ignored, but it is required.
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
    ## Check for user supplied prediction function handle
    if (! is_function_handle (Predfun))
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
      is_cellarray = false;
    else
      X = args(1:end-1);  # cell array with multiple data variables
      n = size (args{1}, 1);
      is_cellarray = true;
    endif
    ## Check for valid prediction function handle on user data
    try
      if (is_cellarray)
        yFit = Predfun (X{:}, y, X{:});
      else
        yFit = Predfun (X, y, X);
      endif
    catch
      error ("crossval: bad prediction function handle for error evaluation.");
    end_try_catch
    if (! iscolumn (yFit) || numel (yFit) != numel (y))
      error (strcat ("crossval: prediction function handle must return", ...
                     " a column vector with the same rows as XTest."));
    endif

  elseif (is_function_handle (f))
    ## At least one additional input argument (X) is required
    nargs = numel (args);
    if (nargs < 1)
      error ("crossval: X is required for values evaluation.");
    endif
    if (nargs == 1)
      X = args{1};        # numeric matrix with single data variable
      n = size (X, 1);
      is_cellarray = false;
    else
      X = args(1:end-1);  # cell array with multiple data variables
      n = size (args{1}, 1);
      is_cellarray = true;
    endif
    ## Check for valid function handle on user data
    try
      if (is_cellarray)
        value = f (X{:}, X{:});
      else
        value = f (X, X);
      endif
    catch
      error ("crossval: bad function handle to cross-validate.");
    end_try_catch
    if (isscalar (value))
      is_scalar = true;
    elseif (isrow (value))
      is_scalar = false;
    else
      error ("crossval: function handle must return a scalar or a row vector.");
    endif

  else
    error ("crossval: invalid first input argument.");
  endif

  ## Input validation for valid values are handled by the cvpartition class.
  ## Check for single paired argument for CV partition.
  vcpa = sum ([isempty(Holdout), isempty(KFold), ...
              isempty(Leaveout), isempty(Partition)]);
  if (vcpa < 3) # at least 3 must be empty
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
    nSets = P.NumTestSets;
    if (P.IsCustom || ismember (P.Type, {'resubstitution', 'leaveout'}))
      MCReps = 1;
    endif
  elseif (! isempty (Leaveout))
    P = cvpartition (n, "LeaveOut");
    nSets = P.NumTestSets;
    MCReps = 1;
  elseif (! isempty (Holdout))
    if (isempty (Stratify))
      P = cvpartition (n, "HoldOut", Holdout);
    else
      P = cvpartition (Stratify, "HoldOut", Holdout);
    endif
    nSets = P.NumTestSets;
  elseif (! isempty (KFold))
    if (isempty (Stratify))
      P = cvpartition (n, "KFold", KFold);
    else
      P = cvpartition (Stratify, "KFold", KFold);
    endif
    nSets = P.NumTestSets;
  else # KFold by default
    if (isempty (Stratify))
      P = cvpartition (n, "KFold");
    else
      P = cvpartition (Stratify, "KFold");
    endif
    nSets = P.NumTestSets;
  endif

  ## Apply cross-validation scheme
  if (ischar (f)) # error evaluation
    results = nan (MCReps, nSets);
    for rep = 1:MCReps
      if (rep > 1)
        P = repartition (P);
      endif
      for idx = 1:nSets
        idx_train = training (P, idx);
        idx_test = test (P, idx);
        if (is_cellarray)
          Xtrain = cellfun (@(x) x(idx_train, :), X, "UniformOutput", false);
          Xtest = cellfun (@(x) x(idx_test, :), X, "UniformOutput", false);
          y_fit = Predfun (Xtrain{:}, y(idx_train), Xtest{:});
        else
          y_fit = Predfun (X(idx_train, :), y(idx_train), X(idx_test, :));
        endif
        if (strcmp (f, "mse"))
          err = sum ((y_fit - y(idx_test)).^2) / numel (y_fit);
          results(rep, idx) = err;
        else # MCR
          err = sum (y_fit == y(idx_test)) / numel (y_fit);
          results(rep, idx) = err;
        endif
      endfor
    endfor
    results = mean (mean (results));

  else  # model execution
    if (is_scalar)
      results = nan (MCReps, nSets);
      for rep = 1:MCReps
        if (rep > 1)
          P = repartition (P);
        endif
        for idx = 1:nSets
          idx_train = training (P, idx);
          idx_test = test (P, idx);
          if (is_cellarray)
            Xtrain = cellfun (@(x) x(idx_train, :), X, "UniformOutput", false);
            Xtest = cellfun (@(x) x(idx_test, :), X, "UniformOutput", false);
            result = f (Xtrain{:}, Xtest{:});
          else
            result = f (X(idx_train, :), X(idx_test, :));
          endif
          results(rep, idx) = result;
        endfor
      endfor
    else  # concatenate Monte Carlo repetitions along first dimension
      results = [];
      for rep = 1:MCReps
        if (rep > 1)
          P = repartition (P);
        endif
        tmpresults = [];
        for idx = 1:nSets
          idx_train = training (P, idx);
          idx_test = test (P, idx);
          if (is_cellarray)
            Xtrain = cellfun (@(x) x(idx_train, :), X, "UniformOutput", false);
            Xtest = cellfun (@(x) x(idx_test, :), X, "UniformOutput", false);
            result = f (Xtrain{:}, Xtest{:});
          else
            result = f (X(idx_train, :), X(idx_test, :));
          endif
          tmpresults = [tmpresults; result];
        endfor
        results = [results; tmpresults];
      endfor
    endif
  endif
endfunction

%!demo
%! ## Determine the optimal number of clusters using cross-validation
%!
%! ## Declare a function to compute the sum of squared distances
%! ## between data points and a varying number of clusters.
%! function D = dist2clusters (X, Y, k)
%!   [Z, Zmu, Zstd] = zscore (X);
%!   [~, C] = kmeans (Z, k);
%!   ZY = (Y - Zmu) ./ Zstd;
%!   d = pdist2 (C, ZY, 'euclidean', 'Smallest', 1);
%!   D = sum (d .^ 2);
%! endfunction
%!
%! load fisheriris
%! for k = 1:8
%!   fcn = @(X, Y) dist2clusters (X, Y, k);
%!   distances = crossval (fcn, meas);
%!   cvdist(k) = sum (distances);
%! endfor
%!
%! plot(cvdist)
%! xlabel('Number of Clusters')
%! ylabel('CV Sum of Squared Distances')
%! xlim ([1,8]);

## Test output
%!test
%! function yfit = regf (Xtrain, ytrain, Xtest)
%!   b = regress (ytrain, Xtrain);
%!   yfit = Xtest * b;
%! endfunction
%!
%! load carsmall
%! data = [Acceleration Horsepower Weight MPG];
%! data(any(isnan(data),2),:) = [];
%!
%! y = data(:,4);
%! X = [ones(length(y),1) data(:,1:3)];
%! rand ("seed", 3);
%! cvMSE = crossval('mse',X,y,'Predfun',@regf);
%! assert (cvMSE, 18.720, 1e-3);

## Test input validation
%!error <crossval: criterion must be 'mse' or 'mcr'.> ...
%! crossval ('fe', rand (10, 1), rand (10, 1), 1);
%!error <crossval: prediction function handle is required for error evaluation.>  ...
%! crossval ('mse', rand (10, 1), rand (10, 1), 1);
%!error <crossval: X and Y are required for error evaluation.> ...
%! crossval ('mse', rand (10, 1), 'Predfun', @(x,y) x + y);
%!error <crossval: bad prediction function handle for error evaluation.> ...
%! crossval ('mse', rand (10, 3), rand (10, 1), 'Predfun', @(x,y) sum (x + y));
%!error <crossval: prediction function handle must return a column vector with the same rows as XTest.> ...
%! crossval ('mse', rand (10, 3), rand (10, 1), 'Predfun', @(x,y,z) sum (x + y));
%!error <crossval: X is required for values evaluation.> crossval (@(x) x);
%!error <crossval: bad function handle to cross-validate.> ...
%! crossval (@(x) x, rand (10, 3), rand (10, 1));
%!error <crossval: function handle must return a scalar or a row vector.> ...
%! crossval (@(x,y) [x, y], rand (10, 3), rand (10, 1));
%!error <crossval: invalid first input argument.> crossval ({1}, 1, 1);
%!error <crossval: you can only set one cvpartition type in paired arguments.> ...
%! crossval (@(x,y) sum ([x; y]), rand (10, 3), 'Holdout', 0.1, 'Leaveout', true)
%!error <crossval: you cannot specify both 'Partition' and 'Stratify'.> ...
%! crossval (@(x,y) sum ([x; y]), rand (10, 3), 'Partition', cvpartition (10, 'Leaveout'), 'Stratify', true)
