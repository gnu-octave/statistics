## Copyright (C) 2025 Andreas Bertsatos <abertsatos@biol.uoa.gr>
##
## This file is part of the statistics package for GNU Octave.
##
## This program is free software; you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation; either version 3 of the License, or (at your option) any later
## version.
##
## This program is distributed in the hope that it will be useful, but WITHOUT
## ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
## FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
## details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.

classdef cvpartition
## -*- texinfo -*-
## @deftp {statistics} cvpartition
##
## Partition data for cross-validation
##
## The @code{cvpartition} class generates a partitioning scheme on a dataset to
## facilitate cross-validation of statistical models utilizing training and
## testing subsets of the dataset.
##
## Create a @code{cvpartition} object by using the class constructor.
##
## @seealso{crossval}
## @end deftp

  properties (SetAccess = private)
    ## -*- texinfo -*-
    ## @deftp {cvpartition} {property} NumObservations
    ##
    ## Number of observations
    ##
    ## A positive integer scalar specifying the number of observations in the
    ## dataset (including any missing data, where applicable).  This property
    ## is read-only.
    ##
    ## @end deftp
    NumObservations = [];

    ## -*- texinfo -*-
    ## @deftp {cvpartition} {property} NumTestSets
    ##
    ## Number of test sets
    ##
    ## A positive integer scalar specifying the number of folds for partition
    ## types @qcode{'kfold'} and @qcode{'leaveout'}.  When partition type is
    ## @qcode{'holdout'} and @qcode{'resubstitution'}, then @qcode{NumTestSets}
    ## is 1.  This property is read-only.
    ##
    ## @end deftp
    NumTestSets     = [];

    ## -*- texinfo -*-
    ## @deftp {cvpartition} {property} TrainSize
    ##
    ## Size of each train set
    ##
    ## A positive integer scalar specifying the size of the train set for
    ## partition types @qcode{'holdout'} and @qcode{'resubstitution'} or a
    ## vector of positive integers specifying the size of each training set for
    ## partition types @qcode{'kfold'} and @qcode{'leaveout'}.  This property
    ## is read-only.
    ##
    ## @end deftp
    TrainSize       = [];

    ## -*- texinfo -*-
    ## @deftp {cvpartition} {property} TestSize
    ##
    ## Size of each test set
    ##
    ## A positive integer scalar specifying the size of the test set for
    ## partition types @qcode{'holdout'} and @qcode{'resubstitution'} or a
    ## vector of positive integers specifying the size of each testing set for
    ## partition types @qcode{'kfold'} and @qcode{'leaveout'}.  This property
    ## is read-only.
    ##
    ## @end deftp
    TestSize        = [];

    ## -*- texinfo -*-
    ## @deftp {cvpartition} {property} Type
    ##
    ## Type of validation partition
    ##
    ## A character vector specifying the type of the @qcode{cvpartition} object.
    ## It can be @qcode{kfold}, @qcode{holdout}, @qcode{leaveout}, or
    ## @qcode{resubstitution}.  This property is read-only.
    ##
    ## @end deftp
    Type            = '';

    ## -*- texinfo -*-
    ## @deftp {cvpartition} {property} IsCustom
    ##
    ## Flag for custom partition
    ##
    ## A logical scalar specifying whether the @qcode{cvpartition} object
    ## was created using custom partition partitioning (@qcode{true}) or
    ## not (@qcode{false}).  This property is read-only.
    ##
    ## @end deftp
    IsCustom        = [];

    ## -*- texinfo -*-
    ## @deftp {cvpartition} {property} IsGrouped
    ##
    ## Flag for grouped partition
    ##
    ## A logical scalar specifying whether the @qcode{cvpartition} object was
    ## created using grouping variables (@qcode{true}) or not (@qcode{false}).
    ## This property is read-only.
    ##
    ## @end deftp
    IsGrouped       = [];

    ## -*- texinfo -*-
    ## @deftp {cvpartition} {property} IsStratified
    ##
    ## Flag for stratified partition
    ##
    ## A logical scalar specifying whether the @qcode{cvpartition} object was
    ## created with a @qcode{'stratify'} value of @qcode{true}.  This property
    ## is read-only.
    ##
    ## @end deftp
    IsStratified    = [];

  endproperties

  properties (SetAccess = private, Hidden)
    missidx = [];
    indices = [];
    cvptype = '';
    classes = [];
    classID = [];
    grpvars = [];
  endproperties

  methods (Hidden)

    ## Custom display
    function display (this)
      in_name = inputname (1);
      if (! isempty (in_name))
        fprintf ('%s =\n', in_name);
      endif
      disp (this);
    endfunction

    ## Custom display
    function disp (this)
      fprintf ("\n%s\n", this.cvptype);
      ## Print selected properties
      fprintf ("%+25s: %d\n", 'NumObservations', this.NumObservations);
      fprintf ("%+25s: %d\n", 'NumTestSets', this.NumTestSets);
      vlen = numel (this.TrainSize);
      if (vlen <= 10)
        str = repmat ({"%d"}, 1, vlen);
        str = strcat ('[', strjoin (str, ' '), ']');
        str1 = sprintf (str, this.TrainSize);
        str2 = sprintf (str, this.TestSize);
      else
        str = repmat ({"%d"}, 1, 10);
        str = strcat ('[', strjoin (str, ' '), ' ... ]');
        str1 = sprintf (str, this.TrainSize(1:10));
        str2 = sprintf (str, this.TestSize(1:10));
      endif
      fprintf ("%+25s: %s\n", 'TrainSize', str1);
      fprintf ("%+25s: %s\n", 'TestSize', str2);
      fprintf ("%+25s: %d\n", 'IsCustom', this.IsCustom);
      fprintf ("%+25s: %d\n", 'IsGrouped', this.IsGrouped);
      fprintf ("%+25s: %d\n\n", 'IsStratified', this.IsStratified);
    endfunction

    ## Class specific subscripted reference
    function varargout = subsref (this, s)
      chain_s = s(2:end);
      s = s(1);
      t = "Invalid %s indexing for referencing values in a cvpartition object.";
      switch (s.type)
        case '()'
          error (t, '()');
        case '{}'
          error (t, '{}');
        case '.'
          if (! ischar (s.subs))
            error (strcat ("cvpartition.subsref: '.' indexing", ...
                           " argument must be a character vector."));
          endif
          try
            out = this.(s.subs);
          catch
            error ("cvpartition.subref: unrecognized property: '%s'", s.subs);
          end_try_catch
      endswitch
      ## Chained references
      if (! isempty (chain_s))
        out = subsref (out, chain_s);
      endif
      varargout{1} = out;
    endfunction

    ## Class specific subscripted assignment
    function this = subsasgn (this, s, val)
      if (numel (s) > 1)
        error (strcat ("cvpartition.subsasgn:", ...
                       " chained subscripts not allowed."));
      endif
      t = "Invalid %s indexing for assigning values to a cvpartition object.";
      switch s.type
        case '()'
          error (t, '()');
        case '{}'
          error (t, '{}');
        case '.'
          if (! ischar (s.subs))
            error (strcat ("cvpartition.subsasgn: '.' indexing", ...
                           " argument must be a character vector."));
          endif
          error (strcat ("cvpartition.subsasgn: unrecognized", ...
                         " or read-only property: '%s'"), s.subs);
      endswitch
    endfunction

  endmethods

  methods (Access = public)

    ## -*- texinfo -*-
    ## @deftypefn  {cvpartition} {@var{C} =} cvpartition (@var{n}, @qcode{'KFold'})
    ## @deftypefnx {cvpartition} {@var{C} =} cvpartition (@var{n}, @qcode{'KFold'}, @var{k})
    ## @deftypefnx {cvpartition} {@var{C} =} cvpartition (@var{n}, @qcode{'KFold'}, @var{k}, @qcode{'GroupingVariables'}, @var{grpvars})
    ## @deftypefnx {cvpartition} {@var{C} =} cvpartition (@var{n}, @qcode{'Holdout'})
    ## @deftypefnx {cvpartition} {@var{C} =} cvpartition (@var{n}, @qcode{'Holdout'}, @var{p})
    ## @deftypefnx {cvpartition} {@var{C} =} cvpartition (@var{n}, @qcode{'Leaveout'})
    ## @deftypefnx {cvpartition} {@var{C} =} cvpartition (@var{n}, @qcode{'Resubstitution'})
    ## @deftypefnx {cvpartition} {@var{C} =} cvpartition (@var{X}, @qcode{'KFold'})
    ## @deftypefnx {cvpartition} {@var{C} =} cvpartition (@var{X}, @qcode{'KFold'}, @var{k})
    ## @deftypefnx {cvpartition} {@var{C} =} cvpartition (@var{X}, @qcode{'KFold'}, @var{k}, @qcode{'Stratify'}, @var{opt})
    ## @deftypefnx {cvpartition} {@var{C} =} cvpartition (@var{X}, @qcode{'Holdout'})
    ## @deftypefnx {cvpartition} {@var{C} =} cvpartition (@var{X}, @qcode{'Holdout'}, @var{p})
    ## @deftypefnx {cvpartition} {@var{C} =} cvpartition (@var{X}, @qcode{'Holdout'}, @var{p}, @qcode{'Stratify'}, @var{opt})
    ## @deftypefnx {cvpartition} {@var{C} =} cvpartition (@qcode{'CustomPartition'}, @var{testSets})
    ##
    ## Repartition data for cross-validation.
    ##
    ## @code{@var{C} = cvpartition (@var{n}, @qcode{'KFold'})} creates a
    ## @qcode{cvpartition} object @var{C}, which defines a random nonstratified
    ## partition for k-fold cross-validation on @var{n} observations with each
    ## fold (subsample) having approximately the same number of observations.
    ## The default number of folds is 10 for @code{@var{n} >= 10} or equal to
    ## @var{n} otherwise.
    ##
    ## @code{@var{C} = cvpartition (@var{n}, @qcode{'KFold'}, @var{k})} also
    ## creates a nonstratified random partition for k-fold cross-validation with
    ## the number of folds defined by @var{k}, which must be a positive integer
    ## scalar smaller than the number of observations @var{n}.
    ##
    ## @code{@var{C} = cvpartition (@var{n}, @qcode{'KFold'}, @var{k},
    ## @qcode{'GroupingVariables'}, @var{grpvars})} creates a @qcode{cvpartition}
    ## object @var{C} that defines a random partition for k-fold cross-validation
    ## with each fold containing the same combination of group labels as defined
    ## by @var{grpvars}.  The grouping variables specified in @var{grpvars} can
    ## be one of the following:
    ##
    ## @itemize
    ## @item A numeric vector, logical vector, categorical vector, character
    ## array, string array, or cell array of character vectors containing one
    ## grouping variable.
    ## @item A numeric matrix or cell array containing two or more grouping
    ## variables. Each column in the matrix or array must correspond to one
    ## grouping variable.
    ## @end itemize
    ##
    ## @code{@var{C} = cvpartition (@var{n}, @qcode{'Holdout'})} creates a
    ## @qcode{cvpartition} object @var{C}, which defines a random nonstratified
    ## partition for holdout validation on @var{n} observations.  90% of the
    ## observations are assigned to the training set and the remaining 10% to
    ## the test set.
    ##
    ## @code{@var{C} = cvpartition (@var{n}, @qcode{'Holdout'}, @var{p})} also
    ## creates a nonstratified random partition for holdout validation with the
    ## percentage of training and test sets defined by @var{p}, which can be a
    ## scalar value in the range @math{(0,1)} or a positive integer scalar in
    ## the range @math{[1,@var{n})}.
    ##
    ## @code{@var{C} = cvpartition (@var{n}, @qcode{'Leaveout'})} creates a
    ## @qcode{cvpartition} object @var{C}, which defines a random partition for
    ## leave-one-out cross-validation on @var{n} observations.  This is a
    ## special case of k-fold cross-validation with the number of folds equal to
    ## the number of observations.
    ##
    ## @code{@var{C} = cvpartition (@var{n}, @qcode{'Resubstitution'})} creates
    ## a @qcode{cvpartition} object @var{C} without partitioning the data and
    ## both training and test sets containing all observations @var{n}.
    ##
    ## @code{@var{C} = cvpartition (@var{X}, @qcode{'KFold'})} creates a
    ## @qcode{cvpartition} object @var{C}, which defines a stratified random
    ## partition for k-fold cross-validation according to the class proportions
    ## in @var{Χ}.  @var{X} can be a numeric, logical, categorical, or string
    ## vector, or a character array or a cell array of character vectors.
    ## Missing values in @var{X} are discarded.  The default number of folds is
    ## 10 for @code{numel (@var{X}) >= 10} or equal to @code{numel (@var{X})}
    ## otherwise.
    ##
    ## @code{@var{C} = cvpartition (@var{X}, @qcode{'KFold'}, @var{k})} also
    ## creates a stratified random partition for k-fold cross-validation with
    ## the number of folds defined by @var{k}, which must be a positive integer
    ## scalar smaller than the number of observations in @var{X}.
    ##
    ## @code{@var{C} = cvpartition (@var{X}, @qcode{'KFold'}, @var{k},
    ## @qcode{'Stratify'}, @var{opt})} creates a random partition for k-fold
    ## cross-validation, which is stratified if @var{opt} is @qcode{true}, or
    ## nonstratified if @var{opt} is @qcode{false}.
    ##
    ## @code{@var{C} = cvpartition (@var{X}, @qcode{'Holdout'})} creates a
    ## @qcode{cvpartition} object @var{C}, which defines a stratified random
    ## partition for holdout validation while maintaining the class proportions
    ## in @var{Χ}.  90% of the observations are assigned to the training set and
    ## the remaining 10% to the test set.
    ##
    ## @code{@var{C} = cvpartition (@var{X}, @qcode{'Holdout'}, @var{p})} also
    ## creates a stratified random partition for holdout validation with the
    ## percentage of training and test sets defined by @var{p}, which can be a
    ## scalar value in the range @math{(0,1)} or a positive integer scalar in
    ## the range @math{[1,@var{n})}.
    ##
    ## @code{@var{C} = cvpartition (@var{X}, @qcode{'Holdout'}, @var{p},
    ## @qcode{'Stratify'}, @var{opt})} creates a random partition for holdout
    ## validation, which is stratified if @var{opt} is @qcode{true}, or
    ## nonstratified if @var{opt} is @qcode{false}.
    ##
    ## @code{@var{C} = cvpartition (@qcode{'CustomPartition'}, @var{testSets})}
    ## creates a custom partition according to @var{testSets}, which can be a
    ## positive integer vector, a logical vector, or a logical matrix according
    ## to the following options:
    ## @itemize
    ## @item A positive integer vector of length @math{N} with values in the
    ## range @math{[1,K]}, where @math{K < N}, will specify a K-fold
    ## cross-validation partition, in which each value indicates the test set
    ## of each observation.  Alternatively, the same vector with values in the
    ## range @math{[1,N]} will specify a leave-one-out cross-validation.
    ## @item A logical vector will specify a holdout validation, in which the
    ## @qcode{true} elements correspond to the test set and the @qcode{false}
    ## elements correspond to the training set.
    ## @item A logical matrix with @math{K} columns will specify a K-fold
    ## cross-validation partition, in which each column corresponds to a fold
    ## and each row to an observation.  Alternatively, an @math{NxN} logical
    ## matrix will specify a leave-one-out cross-validation, where @math{N} is
    ## the number of observations.  @qcode{true} elements correspond to the
    ## test set and the @qcode{false} elements correspond to the training set.
    ## @end itemize
    ##
    ## @seealso{cvpartition, summary, test, training}
    ## @end deftypefn

    function this = cvpartition (X, varargin)

      ## Check for appropriate number of input arguments
      if (nargin < 2)
        error ("cvpartition: too few input arguments.");
      endif
      if (nargin > 5)
        error ("cvpartition: too many input arguments.");
      endif

      ## Check for custom partition
      if (strcmpi (X, "CustomPartition"))
        testSets = varargin{1};
        ## Check for valid test set
        if (! (isnumeric (testSets) || islogical (testSets)))
          error ("cvpartition: TESTSETS must be numeric of logical.");
        endif
        if (isnumeric (testSets))
          if (! isvector (testSets))
            error ("cvpartition: TESTSETS must be a numeric vector.");
          endif
          [~, idx, inds] = unique (testSets);
          this.NumObservations = numel (testSets);
          this.NumTestSets = numel (idx);
          nvec = this.NumObservations * ones (1, this.NumTestSets);
          if (this.NumTestSets < this.NumObservations)
            this.indices = inds;
            for i = 1:this.NumTestSets
              this.TestSize(i) = sum (inds == i);
            endfor
            this.TrainSize = nvec - this.TestSize;
            this.Type = 'kfold';
            this.cvptype = 'K-fold cross validation partition';
          else
            this.TrainSize = nvec - 1;
            this.TestSize = nvec - this.TrainSize;
            this.Type = 'leaveout';
            this.cvptype = 'Leave-one-out cross validation partition';
          endif
        else  # logical vector of matrix
          if (! ismatrix (testSets))
            error ("cvpartition: TESTSETS must be a logical vector or matrix.");
          elseif (isvector (testSets))
            this.NumObservations = numel (testSets);
            this.NumTestSets = 1;
            this.indices = testSets;
            this.TrainSize = sum (! testSets);
            this.TestSize = sum (testSets);
            this.Type = 'holdout';
            this.cvptype = 'Hold-out cross validation partition';
          else  # logical matrix
            ## Each observation must be present in exactly one test set
            if (any (sum (testSets, 2) > 1))
              error (strcat ("cvpartition: each observation in TESTSETS", ...
                             " must be exactly one in each row."));
            endif
            [this.NumObservations, this.NumTestSets] = size (testSets);
            nvec = this.NumObservations * ones (1, this.NumTestSets);
            if (this.NumTestSets < this.NumObservations)
              this.indices = zeros (this.NumObservations, 1);
              for i = 1:this.NumTestSets
                this.TestSize(i) = sum (testSets(:,i));
                this.indices(testSets(:,i)) = i;
              endfor
              this.TrainSize = nvec - this.TestSize;
              this.Type = 'kfold';
              this.cvptype = 'K-fold cross validation partition';
            elseif (this.NumTestSets == this.NumObservations)
              this.TrainSize = nvec - 1;
              this.TestSize = nvec - this.TrainSize;
              this.Type = 'leaveout';
              this.cvptype = 'Leave-one-out cross validation partition';
            else
              error (strcat ("cvpartition: a logical matrix in TESTSETS", ...
                             " must not have more columns that rows."));
            endif
          endif
        endif
        this.IsCustom = true;
        this.IsGrouped = false;
        this.IsStratified = false;

      ## Check first input being a scalar value
      elseif (isscalar (X))
        if (! (isnumeric (X) && X > 0 && fix (X) == X))
          error ("cvpartition: X must be a scalar positive integer value.");
        endif
        ## Get number of observations and partition type
        this.NumObservations = X;
        type = varargin{1};
        this.IsCustom = false;
        this.IsStratified = false;

        ## "Resubstitution"
        if (strcmpi (type, 'resubstitution'))
          this.NumTestSets = 1;
          this.TrainSize = X;
          this.TestSize = X;
          this.Type = 'resubstitution';
          this.cvptype = 'Resubstitution (no partition of data)';
          this.IsGrouped = false;

        ## "Leaveout"
        elseif (strcmpi (type, 'leaveout'))
          this.NumTestSets = X;
          this.TrainSize = (X - 1) * ones (1, X);
          this.TestSize = ones (1, X);
          this.Type = 'leaveout';
          this.cvptype = 'Leave-one-out cross validation partition';
          this.IsGrouped = false;

        ## "Holdout"
        elseif (strcmpi (type, 'holdout'))
          if (nargin > 2)
            p = varargin{2};
            if (! isnumeric (p) || ! isscalar (p))
              error (strcat ("cvpartition: P value for 'holdout'", ...
                             " must be a numeric scalar."));
            endif
            if (! ((p > 0 && p < 1) || (p == fix (p) && p > 0 && p < X)))
              error (strcat ("cvpartition: P value for 'holdout' must be", ...
                             " a scalar in the range (0,1) or an integer", ...
                             " scalar in the range [1, N)."));
            endif
          else
            p = 0.1;
          endif
          this.NumTestSets = 1;
          if (p < 1)            # target fraction to sample
            p = round (p * X);  # number of samples
          endif
          inds = false (X, 1);
          inds(randsample (X, p)) = true;  # indices for test set
          this.indices = inds;
          this.TrainSize = sum (! inds);
          this.TestSize = sum (inds);
          this.Type = 'holdout';
          this.cvptype = 'Hold-out cross validation partition';
          this.IsGrouped = false;

        ## "KFold"
        elseif (strcmpi (type, 'kfold'))
          this.Type = 'kfold';
          if (nargin > 2)
            k = varargin{2};
            if (! isnumeric (k) || ! isscalar (k))
              error (strcat ("cvpartition: K value for 'kfold'", ...
                             " must be a numeric scalar."));
            endif
          else
            if (X < 10)
              k = X;
            else
              k = 10;
            endif
          endif
          ## No grouping variables
          if (nargin < 4)
            if (! (k == fix (k) && k > 0 && k <= X))
              error (strcat ("cvpartition: K value for 'kfold' must be", ...
                             " an integer scalar in the range [1, N]."));
            endif
            this.NumTestSets = k;
            indices = floor ((0:(X - 1))' * (k / X)) + 1;
            indices = randsample (indices, X);
            nvec = X * ones (1, k);
            for i = 1:k
              this.TestSize(i) = sum (indices == i);
            endfor
            this.indices = indices;
            this.TrainSize = nvec - this.TestSize;
            this.cvptype = 'K-fold cross validation partition';
            this.IsGrouped = false;
          else  # with grouping variables
            if (! strcmpi (varargin{3}, 'groupingvariables'))
              error (strcat ("cvpartition: invalid optional paired", ...
                             " argument for 'GroupingVariables'."));
            endif
            if (nargin < 5)
              error (strcat ("cvpartition: missing value for optional", ...
                             " paired argument 'GroupingVariables'."));
            endif
            grpvars = varargin{4};
            if (isvector (grpvars))
              ## Remove any missing values
              this.missidx = ismissing (grpvars);
              if (any (this.missidx))
                grpvars(this.missidx) = [];
                X -= sum (this.missidx);
              endif
              ## Get indices for each group
              if (isa (grpvars, 'categorical'))
                [~, idx, inds] = unique (grpvars, 'stable');
              else
                [~, idx, inds] = __unique__ (grpvars, 'stable');
              endif
            elseif (ismatrix (grpvars))
              ## Remove any missing values
              this.missidx = any (ismissing (grpvars), 2);
              if (any (this.missidx))
                grpvars(this.missidx, :) = [];
                X -= sum (this.missidx);
              endif
              ## Get indices for each group
              if (isa (grpvars, 'categorical'))
                [~, idx, inds] = unique (grpvars, 'rows', 'stable');
              else
                [~, idx, inds] = __unique__ (grpvars, 'rows', 'stable');
              endif
            else
              error (strcat ("cvpartition: invalid value for optional", ...
                             " paired argument 'GroupingVariables'."));
            endif
            if (X != numel (inds))
              error (strcat ("cvpartition: grouping variable does", ...
                             " not match the number of observations."));
            endif
            this.grpvars = grpvars;
            ## Get number of groups and group sizes
            NumGroups = numel (idx);
            for i = 1:NumGroups
              GroupSize(i) = sum (inds == i);
            endfor
            ## Compare k-fold to number of groups and reduce K accordingly
            if (k > NumGroups)
              warning (strcat ("cvpartition: number of folds K is greater", ...
                               " than the groups in 'GroupingVariables'.", ...
                               " K is set to the number of groups."));
                k = NumGroups;
            endif
            ## If k == NumGroups, then each group becomes a test in a fold.
            ## If k < NumGroups, then cluster NumGroups to k folds.
            indices = zeros (X, 1);
            if (k == NumGroups)
              for i = 1:k
                indices(inds == i) = i;
              endfor
            else
              [GroupIdx, ~, GroupSz] = multiway (GroupSize, k, 'completeKK');
              for i = 1:k
                idxGV = find (GroupIdx == i);
                vecGV = arrayfun(@(x) x == inds, idxGV, "UniformOutput", false);
                index = vecGV{1};
                if (numel (vecGV) > 1)
                  for j = 2:numel (vecGV)
                    index = index | vecGV{j};
                  endfor
                endif
                indices(index) = i;
              endfor
            endif
            ## Randomize the order of folds
            random_idx = randsample ([1:k], k);
            randomized = zeros (size (inds));
            for i = 1:k
              randomized(indices == i) = random_idx(i);
            endfor
            ## Save values to properties
            this.indices = randomized;
            this.NumTestSets = k;
            nvec = X * ones (1, k);
            for i = 1:k
              this.TestSize(i) = sum (this.indices == i);
            endfor
            this.TrainSize = nvec - this.TestSize;
            this.cvptype = 'Group K-fold cross validation partition';
            this.IsGrouped = true;
          endif

        ## Invalid paired argument
        else
          error ("cvpartition: invalid optional paired argument.");
        endif

      ## Check first input being a vector for stratification
      elseif (isvector (X))
        ## Get number of observations (including missing values)
        this.NumObservations = numel (X);

        ## Remove missing values from partitioning.
        ## Keep missing index to include them in the test indices.
        this.missidx = ismissing (X);
        X(this.missidx) = [];

        ## Get stratify option
        if (nargin < 4)
          this.IsStratified = true;
        else
          if (! strcmpi (varargin{3}, 'stratify'))
              error (strcat ("cvpartition: invalid optional paired", ...
                             " argument for stratification."));
          endif
          if (nargin < 5)
            error (strcat ("cvpartition: missing value for optional", ...
                           " paired argument 'stratify'."));
          endif
          if (! isscalar (varargin{4}) || ! islogical (varargin{4}))
            error (strcat ("cvpartition: invalid value for optional", ...
                           " paired argument 'stratify'."));
          endif
          this.IsStratified = varargin{4};
        endif

        ## Handle stratification
        if (this.IsStratified)
          [classID, idx, classes] = unique (X);
          NumClasses = numel (idx);
          for i = 1:NumClasses
            ClassSize(i) = sum (classes == i);
          endfor
          this.classes = classes;
          this.classID = classID;
        endif
        X = numel (X);

        ## Get partition type
        type = varargin{1};
        this.IsCustom = false;
        this.IsGrouped = false;

        ## "Holdout"
        if (strcmpi (type, 'holdout'))
          this.Type = 'holdout';
          if (nargin > 2)
            p = varargin{2};
            if (! isnumeric (p) || ! isscalar (p))
              error (strcat ("cvpartition: P value for 'holdout'", ...
                             " must be a numeric scalar."));
            endif
            if (! ((p > 0 && p < 1) || (p == fix (p) && p > 0 && p < X)))
              error (strcat ("cvpartition: P value for 'holdout' must be", ...
                             " a scalar in the range (0,1) or an integer", ...
                             " scalar in the range [1, N), where N is the", ...
                             " number of nonmissing observations in X."));
            endif
          else
            p = 0.1;
          endif
          this.NumTestSets = 1;
          if (this.IsStratified)
            if (p < 1)
              f = p;              # target fraction to sample
              p = round (p * X);  # number of test samples
            else
              f = p / X;
            endif
            inds = zeros (X, 1, "logical");
            k_check = 0;
            for i = 1:NumClasses
              ki = round (f * ClassSize(i));
              inds(find (classes == i)(randsample (ClassSize(i), ki))) = true;
              k_check += ki;
            endfor
            if (k_check < p)      # add random elements to test set to make it p
              inds(find (! inds)(randsample (X - k_check, p - k_check))) = true;
            elseif (k_check > p)  # remove random elements from test set
              inds(find (inds)(randsample (k_check, k_check - p))) = false;
            endif
            this.cvptype = 'Stratified hold-out cross validation partition';
          else
            if (p < 1)            # target fraction to sample
              p = round (p * X);  # number of samples
            endif
            inds = false (X, 1);
            inds(randsample (X, p)) = true;  # indices for test set
            this.cvptype = 'Hold-out cross validation partition';
          endif
          this.indices = inds;
          this.TrainSize = sum (! inds);
          this.TestSize = sum (inds);

        ## "KFold"
        elseif (strcmpi (type, 'kfold'))
          this.Type = 'kfold';
          if (nargin > 2)
            k = varargin{2};
            if (! isnumeric (k) || ! isscalar (k))
              error (strcat ("cvpartition: K value for 'kfold'", ...
                             " must be a numeric scalar."));
            endif
            if (! (k == fix (k) && k > 0 && k <= X))
              error (strcat ("cvpartition: K value for 'kfold' must be", ...
                             " an integer scalar in the range [1, N],", ...
                             " where N is the number of nonmissing", ...
                             " observations in X."));
            endif
          else
            if (X < 10)
              k = X;
            else
              k = 10;
            endif
          endif
          this.NumTestSets = k;
          if (this.IsStratified)
            inds = nan (X, 1);
            pooled_idx = false (X, 1);
            do_warn = true;
            do_ceil = false;
            for i = 1:NumClasses
              cls_size = ClassSize(i);
              cls_k_eq = fix (cls_size / k) == (cls_size / k);
              ## Check that the elements in each class exceed the number of
              ## requested folds, otherwise emit a warning and add the class
              ## elements into a pooled class
              if (cls_size < k)
                if (do_warn)
                  warning (strcat ("One or more of the unique class values", ...
                                   " in the stratification variable is not", ...
                                   " present in one or more folds."));
                  do_warn = false;
                endif
                pooled_idx = pooled_idx | classes == i;
              elseif (fix (X / k) == X / k)
                ## Make sure that when X / k = integer, all
                ## test/training sizes must be equal across all folds
                if (do_ceil && ! cls_k_eq)
                  idx = ceil ((0:(cls_size - 1))' * (k / cls_size));
                  idx(idx == 0) = max (idx);
                  do_ceil = false;
                else
                  idx = floor ((0:(cls_size - 1))' * (k / cls_size)) + 1;
                  tmp = arrayfun (@(x) numel (find (x == idx)), [1:k]);
                  if (any (diff (tmp)))
                    do_ceil = true;
                  endif
                endif
                inds(classes == i) = randsample (idx, cls_size);
              else
                ## Alternate ordering over classes so that
                ## the subsets are more nearly the same size
                if (! do_ceil || cls_k_eq)
                  idx = floor ((0:(cls_size - 1))' * (k / cls_size)) + 1;
                  if (! cls_k_eq)
                    do_ceil = true;
                  endif
                else
                  idx = floor (((cls_size - 1):-1:0)' * (k / cls_size)) + 1;
                  do_ceil = false;
                endif
                inds(classes == i) = randsample (idx, cls_size);
              endif
            endfor
            ## Stratify pooled classes (if any).  They must be distributed
            ## in a way to make the test/training sizes as equal as possible
            ## across folds.
            pooled_inds = find (pooled_idx);
            while (numel (pooled_inds) > 0)
              tmp = arrayfun (@(x) numel (find (x == inds)), [1:k]);
              [min_cls, min_idx] = min (tmp);
              [max_cls, max_idx] = max (tmp);
              if (min_cls != max_cls)
                inds(pooled_inds(1)) = min_idx;
              else
                inds(pooled_inds(1)) = randsample (k, 1);
              endif
              pooled_inds(1) = [];
            endwhile
            this.cvptype = 'Stratified K-fold cross validation partition';
          else
            inds = floor ((0:(X - 1))' * (k / X)) + 1;
            inds = randsample (inds, X);
            this.cvptype = 'K-fold cross validation partition';
          endif
          this.indices = inds;
          nvec = X * ones (1, k);
          for i = 1:k
            this.TestSize(i) = sum (inds == i);
          endfor
          this.TrainSize = nvec - this.TestSize;

        ## Invalid paired argument
        else
          error ("cvpartition: invalid optional paired argument.");
        endif

      ## Otherwise first input is invalid
      else
        error ("cvpartition: invalid first input argument.");
      endif

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {cvpartition} {@var{Cnew} =} repartition (@var{C})
    ## @deftypefnx {cvpartition} {@var{Cnew} =} repartition (@var{C}, @var{sval})
    ## @deftypefnx {cvpartition} {@var{Cnew} =} repartition (@var{C}, @qcode{'legacy'})
    ##
    ## Repartition data for cross-validation.
    ##
    ## @code{@var{Cnew} = repartition (@var{C})} creates a @qcode{cvpartition}
    ## object @var{Cnew} that defines a new random partition of the same type as
    ## the @qcode{cvpartition} @var{C}.
    ##
    ## @code{@var{Cnew} = repartition (@var{C}, @var{sval})} also uses the value
    ## of @var{sval} to set the state of the random generator used in
    ## repartitioning @var{C}.  If @var{sval} is a vector, then the random
    ## generator is set using the @qcode{"state"} keyword as in
    ## @code{rand ("state", @var{sval})}.  If @var{sval} is a scalar, then the
    ## @qcode{"seed"} keyword is used as in @code{rand ("seed", @var{sval})} to
    ## specify that old generators should be used.
    ##
    ## @code{@var{Cnew} = repartition (@var{C}, @qcode{'legacy'})} only applies
    ## to @qcode{cvpartition} objects @var{C} that use k-fold partitioning and
    ## it will repartition @var{C} in the same non-random manner that was
    ## previously used by the old-style @qcode{cvpartition} class of the
    ## statistics package.  The @qcode{'legacy'} option does not apply to
    ## stratified or grouped partitions.
    ##
    ## @seealso{cvpartition, summary, test, training}
    ## @end deftypefn

    function this = repartition (this, sval = [])

      ## Emit error for custom partitions
      if (this.IsCustom)
        error ("cvpartition.repartition: cannot repartition a custom partition.");
      endif

      ## Handle legacy code with no randomization of kfold option
      if (strcmpi (sval, "legacy"))
        if (strcmpi (this.Type, "kfold"))
          X = this.NumObservations;
          k = this.NumTestSets;
          if (! (this.IsGrouped || this.IsStratified))
            inds = floor ((0:(X - 1))' * (k / X)) + 1;
            this.indices = inds;
            nvec = X * ones (1, k);
            for i = 1:k
              this.TestSize(i) = sum (inds == i);
            endfor
            this.TrainSize = nvec - this.TestSize;
          else  # legacy option does not apply for grouped or stratified
            error (strcat ("cvpartition.repartition: 'legacy' flag does", ...
                           " not apply to stratified or grouped 'kfold'", ...
                           " partitioned objects."));
          endif
          return;
        else
          error (strcat ("cvpartition.repartition: 'legacy' flag is only", ...
                         " valid for 'kfold' partitioned objects."));
        endif
      endif

      ## Check sval
      if (! isempty (sval))
        if (! (isvector (sval) && isnumeric (sval) && isreal (sval)))
          error (strcat ("cvpartition.repartition: SVAL must be", ...
                         " a real scalar or vector."));
        endif
        if (isscalar (sval))
          rand ("sval", sval);
        else
          rand ("state", sval);
        endif
      endif

      ## Handle repartitioning of randomized holdout and kfold options
      if (strcmpi (this.Type, "holdout"))
        p = this.TestSize;
        if (this.IsStratified)
          X = sum (! this.missidx);
          inds = false (X, 1);
          NumClasses = numel (this.classID);
          classes = this.classes;
          for i = 1:NumClasses
            ClassSize(i) = sum (classes == i);
          endfor
          f = p / X;
          k_check = 0;
          for i = 1:NumClasses
            ki = round (f * ClassSize(i));
            inds(find (classes == i)(randsample (ClassSize(i), ki))) = true;
            k_check += ki;
          endfor
          if (k_check < p)      # add random elements to test set to make it p
            inds(find (! inds)(randsample (X - k_check, p - k_check))) = true;
          elseif (k_check > p)  # remove random elements from test set
            inds(find (inds)(randsample (k_check, k_check - p))) = false;
          endif
        else
          X = this.NumObservations;
          inds = false (X, 1);
          inds(randsample (X, p)) = true;  # indices for test set
        endif
        this.indices = inds;

      elseif (strcmpi (this.Type, "kfold"))
        k = this.NumTestSets;
        if (! (this.IsGrouped || this.IsStratified))
          X = this.NumObservations;
          inds = floor ((0:(X - 1))' * (k / X)) + 1;
          inds = randsample (inds, X);
          this.indices = inds;
          nvec = X * ones (1, k);
          for i = 1:k
            this.TestSize(i) = sum (inds == i);
          endfor
          this.TrainSize = nvec - this.TestSize;
        elseif (this.IsGrouped)
          ## We only need resample the order of folds in this case
          ## Randomize the order of folds
          random_idx = randsample ([1:k], k);
          randomized = zeros (size (this.indices));
          for i = 1:k
            randomized(this.indices == i) = random_idx(i);
          endfor
          ## Save values to properties
          this.indices = randomized;
          this.NumTestSets = k;
          nvec = sum (! this.missidx) * ones (1, k);
          for i = 1:k
            this.TestSize(i) = sum (this.indices == i);
          endfor
          this.TrainSize = nvec - this.TestSize;
        else  # is stratified
          X = sum (! this.missidx);
          NumClasses = numel (this.classID);
          classes = this.classes;
          for i = 1:NumClasses
            ClassSize(i) = sum (classes == i);
          endfor
          inds = nan (X, 1);
          pooled_idx = false (X, 1);
          do_warn = true;
          do_ceil = false;
          for i = 1:NumClasses
            cls_size = ClassSize(i);
            cls_k_eq = fix (cls_size / k) == (cls_size / k);
            ## Check that the elements in each class exceed the number of
            ## requested folds, otherwise emit a warning and add the class
            ## elements into a pooled class
            if (cls_size < k)
              if (do_warn)
                warning (strcat ("One or more of the unique class values", ...
                                 " in the stratification variable is not", ...
                                 " present in one or more folds."));
                do_warn = false;
              endif
              pooled_idx = pooled_idx | classes == i;
            elseif (fix (X / k) == X / k)
              ## Make sure that when X / k = integer, all
              ## test/training sizes must be equal across all folds
              if (do_ceil && ! cls_k_eq)
                idx = ceil ((0:(cls_size - 1))' * (k / cls_size));
                idx(idx == 0) = max (idx);
                do_ceil = false;
              else
                idx = floor ((0:(cls_size - 1))' * (k / cls_size)) + 1;
                tmp = arrayfun (@(x) numel (find (x == idx)), [1:k]);
                if (any (diff (tmp)))
                  do_ceil = true;
                endif
              endif
              inds(classes == i) = randsample (idx, cls_size);
            else
              ## Alternate ordering over classes so that
              ## the subsets are more nearly the same size
              if (! do_ceil || cls_k_eq)
                idx = floor ((0:(cls_size - 1))' * (k / cls_size)) + 1;
                if (! cls_k_eq)
                  do_ceil = true;
                endif
              else
                idx = floor (((cls_size - 1):-1:0)' * (k / cls_size)) + 1;
                do_ceil = false;
              endif
              inds(classes == i) = randsample (idx, cls_size);
            endif
          endfor
          ## Stratify pooled classes (if any).  They must be distributed
          ## in a way to make the test/training sizes as equal as possible
          ## across folds.
          pooled_inds = find (pooled_idx);
          while (numel (pooled_inds) > 0)
            tmp = arrayfun (@(x) numel (find (x == inds)), [1:k]);
            [min_cls, min_idx] = min (tmp);
            [max_cls, max_idx] = max (tmp);
            if (min_cls != max_cls)
              inds(pooled_inds(1)) = min_idx;
            else
              inds(pooled_inds(1)) = randsample (k, 1);
            endif
            pooled_inds(1) = [];
          endwhile
          this.indices = inds;
          nvec = X * ones (1, k);
          for i = 1:k
            this.TestSize(i) = sum (inds == i);
          endfor
          this.TrainSize = nvec - this.TestSize;
        endif
      endif

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {cvpartition} {@var{tbl} =} summary (@var{C})
    ##
    ## Summarize cross-validation partition.
    ##
    ## @code{@var{tbl} = summary (@var{C})} returns a summary table @var{tbl} of
    ## the @qcode{cvpartition} object @var{C} as long as its type is either
    ## k-fold or holdout and it is either stratified of grouped.  This function
    ## requires support for the @qcode{table} class, which is provided by the
    ## @qcode{datatypes} package.
    ##
    ## @seealso{cvpartition, repartition, test, training}
    ## @end deftypefn

    function tbl = summary (this)
      error ("cvpartition.summary: not unimplemented yet.");
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {cvpartition} {@var{idx} =} test (@var{C})
    ## @deftypefnx {cvpartition} {@var{idx} =} test (@var{C}, @var{i})
    ## @deftypefnx {cvpartition} {@var{idx} =} test (@var{C}, @qcode{"all"})
    ##
    ## Test indices for cross-validation.
    ##
    ## @code{@var{idx} = test (@var{C})} returns a logical vector @var{idx} with
    ## @qcode{true} values indicating the elements corresponding to the test
    ## set defined in the code{cvpartition} object @var{C}.  For K-fold and
    ## leave-one-out partitions, the indices corresponding to the first test set
    ## are returned.
    ##
    ## @code{@var{idx} = test (@var{C}, @var{i})} returns a logical vector or
    ## matrix with the indices of the test set indicated by @var{i}.  If @var{i}
    ## is a scalar, then @var{idx} is a logical vector with the indices of the
    ## @math{i-th} set.  If @var{i} is a vector, then @var{idx} is a logical
    ## matrix in which @code{@var{idx}(:,j)} specified the observations in the
    ## test set @code{@var{i}(j)}.  The value(s) in @var{i} must not exceed the
    ## number of tests in the @qcode{cvpartition} object @var{C}.
    ##
    ## @code{@var{idx} = test (@var{C}, @qcode{"all"})} returns a logical vector
    ## or matrix for all test sets defined in the @qcode{cvpartition} object
    ## @var{C}.  For holdout and resubstitution partition types, a vector is
    ## returned.  For K-fold and leave-one-out, a matrix is returned.
    ##
    ## @seealso{cvpartition, repartition, summary, training}
    ## @end deftypefn

    function idx = test (this, varargin)

      ## Check for sufficient input arguments
      if (nargin > 2)
        error ("cvpartition.test: too many input arguments.");
      elseif (nargin == 2)
        i = varargin{1};
        if (strcmpi (i, "all"))
          idx = logical ([]);
          switch (this.Type)
            case "kfold"
              for i = 1:this.NumTestSets
                if (this.IsStratified || this.IsGrouped)
                  cid = false (this.NumObservations, 1);
                  cid(! this.missidx) = this.indices == i;
                else
                  cid = this.indices == i;
                endif
                idx = [idx, cid];
              endfor
            case "leaveout" # no stratification
              for i = 1:this.NumTestSets
                cid = false (this.NumObservations, 1);
                cid(i) = true;
                idx = [idx, cid];
              endfor
            case "holdout"
              if (this.IsStratified)
                idx = false (this.NumObservations, 1);
                idx(! this.missidx) = this.indices;
              else
                idx = this.indices;
              endif
              idx = this.indices;
            case "resubstitution" # no stratification
              idx = true (this.NumObservations, 1);
          endswitch
          return
        elseif (isempty (i))
          i = 1;
        endif
      else
        i = 1;
      endif

      if (! (isvector (i) && isnumeric (i) &&
             all (fix (i) == i) && all (i > 0)))
        error ("cvpartition.test: set index must be a positive integer vector.");
      elseif (any (i > this.NumTestSets))
        error ("cvpartition.test: set index exceeds 'NumTestSets'.");
      endif

      switch (this.Type)
        case  "kfold"
          if (isscalar (i))
            if (this.IsStratified || this.IsGrouped)
              idx = false (this.NumObservations, 1);
              idx(! this.missidx) = this.indices == i;
            else
              idx = this.indices == i;
            endif
          else
            idx = logical ([]);
            for j = i
              if (this.IsStratified || this.IsGrouped)
                cid = false (this.NumObservations, 1);
                cid(! this.missidx) = this.indices == i;
              else
                cid = this.indices == i;
              endif
              idx = [idx, cid];
            endfor
          endif
        case "leaveout" # no stratification
          if (isscalar (i))
            idx = false (this.NumObservations, 1);
            idx(i) = true;
          else
            idx = logical ([]);
            for j = i
              new = false (this.NumObservations, 1);
              new(j) = true;
              idx = [idx, new];
            endfor
          endif
        case "holdout"
          if (this.IsStratified)
            idx = false (this.NumObservations, 1);
            idx(! this.missidx) = this.indices;
          else
            idx = this.indices;
          endif
        case "resubstitution" # no stratification
          idx = true (this.NumObservations, 1);
      endswitch

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {cvpartition} {@var{idx} =} training (@var{C})
    ## @deftypefnx {cvpartition} {@var{idx} =} training (@var{C}, @var{i})
    ## @deftypefnx {cvpartition} {@var{idx} =} training (@var{C}, @qcode{"all"})
    ##
    ## Training indices for cross-validation.
    ##
    ## @code{@var{idx} = training (@var{C})} returns a logical vector @var{idx}
    ## with @qcode{true} values indicating the elements corresponding to the
    ## training set defined in the code{cvpartition} object @var{C}.  For K-fold
    ## and leave-one-out partitions, the indices corresponding to the first
    ## training set are returned.
    ##
    ## @code{@var{idx} = training (@var{C}, @var{i})} returns a logical vector
    ## or matrix with the indices of the training set indicated by @var{i}.  If
    ## @var{i} is a scalar, then @var{idx} is a logical vector with the indices
    ## of the @math{i-th} set.  If @var{i} is a vector, then @var{idx} is a
    ## logical matrix in which @code{@var{idx}(:,j)} specified the observations
    ## in the training set @code{@var{i}(j)}.  The value(s) in @var{i} must not
    ## exceed the number of tests in the @qcode{cvpartition} object @var{C}.
    ##
    ## @code{@var{idx} = training (@var{C}, @qcode{"all"})} returns a logical
    ## vector or matrix for all training sets defined in the @qcode{cvpartition}
    ## object @var{C}.  For holdout and resubstitution partition types, a vector
    ## is returned.  For K-fold and leave-one-out, a matrix is returned.
    ##
    ## @seealso{cvpartition, repartition, summary, test}
    ## @end deftypefn

    function idx = training (this, varargin)

      ## Check for sufficient input arguments
      if (nargin > 2)
        error ("cvpartition.training: too many input arguments.");
      elseif (nargin == 2)
        i = varargin{1};
        if (strcmpi (i, "all"))
          idx = logical ([]);
          switch (this.Type)
            case "kfold"
              for i = 1:this.NumTestSets
                if (this.IsStratified || this.IsGrouped)
                  cid = false (this.NumObservations, 1);
                  cid(! this.missidx) = this.indices != i;
                else
                  cid = this.indices != i;
                endif
                idx = [idx, cid];
              endfor
            case "leaveout" # no stratification
              for i = 1:this.NumTestSets
                cid = true (this.NumObservations, 1);
                cid(i) = false;
                idx = [idx, cid];
              endfor
            case "holdout"
              if (this.IsStratified)
                idx = false (this.NumObservations, 1);
                idx(! this.missidx) = ! this.indices;
              else
                idx = ! this.indices;
              endif
            case "resubstitution" # no stratification
              idx = true (this.NumObservations, 1);
          endswitch
          return
        elseif (isempty (i))
          i = 1;
        endif
      else
        i = 1;
      endif

      if (! (isvector (i) && isnumeric (i) &&
             all (fix (i) == i) && all (i > 0)))
        error (strcat ("cvpartition.training: set index must", ...
                       " be a positive integer vector."));
      elseif (any (i > this.NumTestSets))
        error ("cvpartition.training: set index exceeds 'NumTestSets'.");
      endif

      switch (this.Type)
        case  "kfold"
          if (isscalar (i))
            if (this.IsStratified || this.IsGrouped)
              idx = false (this.NumObservations, 1);
              idx(! this.missidx) = this.indices != i;
            else
              idx = this.indices != i;
            endif
          else
            idx = logical ([]);
            for j = i
              if (this.IsStratified || this.IsGrouped)
                cid = false (this.NumObservations, 1);
                cid(! this.missidx) = this.indices != i;
              else
                cid = this.indices != i;
              endif
              idx = [idx, cid];
            endfor
          endif
        case "leaveout" # no stratification
          if (isscalar (i))
            idx = true (this.NumObservations, 1);
            idx(i) = false;
          else
            idx = logical ([]);
            for j = i
              new = true (this.NumObservations, 1);
              new(j) = false;
              idx = [idx, new];
            endfor
          endif
        case "holdout"
          if (this.IsStratified)
            idx = false (this.NumObservations, 1);
            idx(! this.missidx) = ! this.indices;
          else
            idx = ! this.indices;
          endif
        case "resubstitution" # no stratification
          idx = true (this.NumObservations, 1);
      endswitch

    endfunction

  endmethods

endclassdef

## Test output results for custom partition
%!test
%! custom = [1, 1, 1, 2, 2, 2, 1, 2, 3, 2, 3, 3, 2, 1, 3]';
%! cv = cvpartition ('CustomPartition', custom);
%! assert (cv.Type, 'kfold');
%! assert (cv.NumObservations, 15);
%! assert (cv.NumTestSets, 3);
%! assert (cv.TrainSize, [10, 9, 11]);
%! assert (cv.TestSize, [5, 6, 4]);
%! assert (cv.IsCustom, true);
%! assert (cv.IsGrouped, false);
%! assert (cv.IsStratified, false);
%! idx = training (cv, 1);
%! assert (idx, custom != 1);
%! idx = test (cv, 1);
%! assert (idx, custom == 1);
%! idx = training (cv, 2);
%! assert (idx, custom != 2);
%! idx = test (cv, 2);
%! assert (idx, custom == 2);
%! idx = training (cv, 3);
%! assert (idx, custom != 3);
%! idx = test (cv, 3);
%! assert (idx, custom == 3);
%! idx1 = training (cv, 'all');
%! idx2 = test (cv, 'all');
%! assert (idx1, ! idx2);
%!test
%! custom = logical ([1, 1, 1, 0, 0, 0, 1, 0, 1, 1])';
%! cv = cvpartition ('CustomPartition', custom);
%! assert (cv.Type, 'holdout');
%! assert (cv.NumObservations, 10);
%! assert (cv.NumTestSets, 1);
%! assert (cv.TrainSize, 4);
%! assert (cv.TestSize, 6);
%! assert (cv.IsCustom, true);
%! assert (cv.IsGrouped, false);
%! assert (cv.IsStratified, false);
%! idx = training (cv, 1);
%! assert (idx, custom != 1);
%! assert (idx, training (cv, 'all'));
%! idx = test (cv, 1);
%! assert (idx, custom == 1);
%! assert (idx, test (cv, 'all'));
%!test
%! custom = logical ([1, 0, 0; 0, 1, 0; 1, 0, 0; 0, 0, 1]);
%! cv = cvpartition ('CustomPartition', custom);
%! assert (cv.Type, 'kfold');
%! assert (cv.NumObservations, 4);
%! assert (cv.NumTestSets, 3);
%! assert (cv.TrainSize, [2, 3, 3]);
%! assert (cv.TestSize, [2, 1, 1]);
%! assert (cv.IsCustom, true);
%! assert (cv.IsGrouped, false);
%! assert (cv.IsStratified, false);
%! idx = training (cv, 1);
%! assert (idx, custom(:,1) == false);
%! idx = test (cv, 1);
%! assert (idx, custom(:,1) == true);
%! idx = training (cv, 2);
%! assert (idx, custom(:,2) == false);
%! idx = test (cv, 2);
%! assert (idx, custom(:,2) == true);
%! assert (! custom, training (cv, 'all'));
%! assert (custom, test (cv, 'all'));
%!test
%! cv = cvpartition ('CustomPartition', [1:8]);
%! assert (cv.Type, 'leaveout');
%! assert (cv.NumObservations, 8);
%! assert (cv.NumTestSets, 8);
%! assert (cv.TrainSize, [7, 7, 7, 7, 7, 7, 7, 7]);
%! assert (cv.TestSize, [1, 1, 1, 1, 1, 1, 1, 1]);
%! assert (cv.IsCustom, true);
%! assert (cv.IsGrouped, false);
%! assert (cv.IsStratified, false);
%! assert (class (training (cv, 1)), 'logical');
%! assert (sum (training (cv, 1)), 7);
%! assert (sum (training (cv, 'all')), cv.TrainSize);
%! assert (class (test (cv, 1)), 'logical');
%! assert (sum (test (cv, 1)), 1);
%! assert (sum (test (cv, 'all')), cv.TestSize);
%! assert (! training (cv, 'all'), test (cv, 'all'));
%!test
%! cv = cvpartition ('CustomPartition', logical (eye (8)));
%! assert (cv.Type, 'leaveout');
%! assert (cv.NumObservations, 8);
%! assert (cv.NumTestSets, 8);
%! assert (cv.TrainSize, [7, 7, 7, 7, 7, 7, 7, 7]);
%! assert (cv.TestSize, [1, 1, 1, 1, 1, 1, 1, 1]);
%! assert (cv.IsCustom, true);
%! assert (cv.IsGrouped, false);
%! assert (cv.IsStratified, false);
%! assert (class (training (cv, 1)), 'logical');
%! assert (sum (training (cv, 1)), 7);
%! assert (sum (training (cv, 'all')), cv.TrainSize);
%! assert (class (test (cv, 1)), 'logical');
%! assert (sum (test (cv, 1)), 1);
%! assert (sum (test (cv, 'all')), cv.TestSize);
%! assert (! training (cv, 'all'), test (cv, 'all'));

## Test output results for scalar input N
%!test
%! cv = cvpartition (10, 'resubstitution');
%! assert (cv.Type, 'resubstitution');
%! assert (cv.NumObservations, 10);
%! assert (cv.NumTestSets, 1);
%! assert (cv.TrainSize, 10);
%! assert (cv.TestSize, 10);
%! assert (cv.IsCustom, false);
%! assert (cv.IsGrouped, false);
%! assert (cv.IsStratified, false);
%! assert (class (training (cv, 1)), 'logical');
%! assert (sum (training (cv, 1)), 10);
%! assert (training (cv, 'all'), logical (ones (10, 1)));
%! assert (class (test (cv, 1)), 'logical');
%! assert (sum (test (cv, 1)), 10);
%! assert (test (cv, 'all'), logical (ones (10, 1)));
%! assert (test (cv), training (cv));
%!test
%! cv = cvpartition (10, 'leaveout');
%! assert (cv.Type, 'leaveout');
%! assert (cv.NumObservations, 10);
%! assert (cv.NumTestSets, 10);
%! assert (cv.TrainSize, ones (1, 10) * 9);
%! assert (cv.TestSize, ones (1, 10));
%! assert (cv.IsCustom, false);
%! assert (cv.IsGrouped, false);
%! assert (cv.IsStratified, false);
%! assert (class (training (cv, 1)), 'logical');
%! assert (sum (training (cv, 1)), 9);
%! assert (training (cv, 'all'), ! logical (eye (10)));
%! assert (class (test (cv, 1)), 'logical');
%! assert (sum (test (cv, 1)), 1);
%! assert (test (cv, 'all'), logical (eye (10)));
%! assert (test (cv), ! training (cv));
%! assert (test (cv, 'all'), ! training (cv, 'all'));
%!test
%! rand ('seed', 5);  # for reproducibility
%! cv = cvpartition (10, 'holdout', 0.3);
%! assert (cv.Type, 'holdout');
%! assert (cv.NumObservations, 10);
%! assert (cv.NumTestSets, 1);
%! assert (cv.TrainSize, 7);
%! assert (cv.TestSize, 3);
%! assert (cv.IsCustom, false);
%! assert (cv.IsGrouped, false);
%! assert (cv.IsStratified, false);
%! assert (class (training (cv, 1)), 'logical');
%! assert (sum (training (cv, 1)), 7);
%! assert (training (cv, 'all'), logical ([1, 0, 1, 1, 0, 1, 1, 1, 0, 1])');
%! assert (class (test (cv, 1)), 'logical');
%! assert (sum (test (cv, 1)), 3);
%! assert (test (cv, 'all'), logical ([0, 1, 0, 0, 1, 0, 0, 0, 1, 0])');
%! assert (test (cv), ! training (cv));
%! assert (test (cv, 'all'), ! training (cv, 'all'));
%!test
%! cv = cvpartition (10, 'holdout', 4);
%! assert (cv.Type, 'holdout');
%! assert (cv.NumObservations, 10);
%! assert (cv.NumTestSets, 1);
%! assert (cv.TrainSize, 6);
%! assert (cv.TestSize, 4);
%! assert (cv.IsCustom, false);
%! assert (cv.IsGrouped, false);
%! assert (cv.IsStratified, false);
%! assert (class (training (cv, 1)), 'logical');
%! assert (sum (training (cv, 1)), 6);
%! assert (class (test (cv, 1)), 'logical');
%! assert (sum (test (cv, 1)), 4);
%! assert (test (cv), ! training (cv));
%! assert (test (cv, 'all'), ! training (cv, 'all'));
%!test
%! cv = cvpartition (5, 'holdout', 4);
%! assert (cv.Type, 'holdout');
%! assert (cv.NumObservations, 5);
%! assert (cv.NumTestSets, 1);
%! assert (cv.TrainSize, 1);
%! assert (cv.TestSize, 4);
%! assert (sum (test (cv, 1)), 4);
%!test
%! cv = cvpartition (5, 'holdout', 1);
%! assert (cv.Type, 'holdout');
%! assert (cv.NumObservations, 5);
%! assert (cv.NumTestSets, 1);
%! assert (cv.TrainSize, 4);
%! assert (cv.TestSize, 1);
%! assert (sum (test (cv, 1)), 1);
%!test
%! cv = cvpartition (5, 'kfold');
%! assert (cv.Type, 'kfold');
%! assert (cv.NumObservations, 5);
%! assert (cv.NumTestSets, 5);
%!test
%! cv = cvpartition (20, 'kfold');
%! assert (cv.Type, 'kfold');
%! assert (cv.NumObservations, 20);
%! assert (cv.NumTestSets, 10);
%!test
%! cv = cvpartition (10, 'kfold', 5);
%! assert (cv.Type, 'kfold');
%! assert (cv.NumObservations, 10);
%! assert (cv.NumTestSets, 5);
%! assert (cv.TrainSize, [8, 8, 8, 8, 8]);
%! assert (cv.TestSize, [2, 2, 2, 2, 2]);
%! assert (cv.IsCustom, false);
%! assert (cv.IsGrouped, false);
%! assert (cv.IsStratified, false);
%! assert (test (cv, 1), ! training (cv, 1));
%! assert (test (cv, 'all'), ! training (cv, 'all'));
%! assert (size (test (cv, 'all')), [10, 5]);
%!test
%! grpvar = [1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 5, 5];
%! rand ('seed', 5);
%! cv = cvpartition (12, 'kfold', 5, 'GroupingVariables', grpvar);
%! assert (cv.Type, 'kfold');
%! assert (cv.NumObservations, 12);
%! assert (cv.NumTestSets, 5);
%! assert (cv.TrainSize, [10, 10, 10, 8, 10]);
%! assert (cv.TestSize, [2, 2, 2, 4, 2]);
%! assert (cv.IsCustom, false);
%! assert (cv.IsGrouped, true);
%! assert (cv.IsStratified, false);
%! assert (test (cv, 1), ! training (cv, 1));
%! assert (test (cv, 'all'), ! training (cv, 'all'));
%! assert (size (test (cv, 'all')), [12, 5]);
%! assert (sum (test (cv, 'all')), [2, 2, 2, 4, 2]);
%! assert (sum (training (cv, 'all')), [10, 10, 10, 8, 10]);
%!test
%! grpvar = [1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 3, 3];
%! rand ('seed', 5);
%! cv = cvpartition (12, 'kfold', 3, 'GroupingVariables', grpvar);
%! assert (cv.Type, 'kfold');
%! assert (cv.NumObservations, 12);
%! assert (cv.NumTestSets, 3);
%! assert (cv.TrainSize, [9, 10, 5]);
%! assert (cv.TestSize, [3, 2, 7]);
%! assert (cv.IsCustom, false);
%! assert (cv.IsGrouped, true);
%! assert (cv.IsStratified, false);
%! assert (test (cv, 1), ! training (cv, 1));
%! assert (test (cv, 'all'), ! training (cv, 'all'));
%! assert (size (test (cv, 'all')), [12, 3]);
%! assert (sum (test (cv, 'all')), [3, 2, 7]);
%! assert (sum (training (cv, 'all')), [9, 10, 5]);
%!test
%! grpvar = [1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3];
%! rand ('seed', 5);
%! cv = cvpartition (12, 'kfold', 2, 'GroupingVariables', grpvar);
%! assert (cv.Type, 'kfold');
%! assert (cv.NumObservations, 12);
%! assert (cv.NumTestSets, 2);
%! assert (cv.TrainSize, [6, 6]);
%! assert (cv.TestSize, [6, 6]);
%! assert (cv.IsCustom, false);
%! assert (cv.IsGrouped, true);
%! assert (cv.IsStratified, false);
%! assert (test (cv, 1), ! training (cv, 1));
%! assert (test (cv, 'all'), ! training (cv, 'all'));
%! assert (size (test (cv, 'all')), [12, 2]);
%! assert (sum (test (cv, 'all')), [6, 6]);
%! assert (sum (training (cv, 'all')), [6, 6]);
%!test
%! grpvar = [1, 1, 1, 2, 2, 2, 2, NaN, 2, 3, 3, 3];
%! rand ('seed', 5);
%! cv = cvpartition (12, 'kfold', 2, 'GroupingVariables', grpvar);
%! assert (cv.Type, 'kfold');
%! assert (cv.NumObservations, 12);
%! assert (cv.NumTestSets, 2);
%! assert (cv.TrainSize, [6, 5]);
%! assert (cv.TestSize, [5, 6]);
%! assert (cv.IsCustom, false);
%! assert (cv.IsGrouped, true);
%! assert (cv.IsStratified, false);
%! idx = ! isnan (grpvar);
%! assert (test (cv, 1)(idx), ! training (cv, 1)(idx));
%! assert (test (cv, 'all')(idx, :), ! training (cv, 'all')(idx, :));
%! assert (size (test (cv, 'all')), [12, 2]);
%! assert (sum (test (cv, 'all')), [5, 6]);
%! assert (sum (training (cv, 'all')), [6, 5]);
%!test
%! grpvar = [1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 3, 3];
%! rand ('seed', 5);
%! cv = cvpartition (12, 'kfold', 2, 'GroupingVariables', grpvar);
%! assert (cv.Type, 'kfold');
%! assert (cv.NumObservations, 12);
%! assert (cv.NumTestSets, 2);
%! assert (cv.TrainSize, [5, 7]);
%! assert (cv.TestSize, [7, 5]);
%! assert (cv.IsCustom, false);
%! assert (cv.IsGrouped, true);
%! assert (cv.IsStratified, false);
%! assert (test (cv, 1), ! training (cv, 1));
%! assert (test (cv, 'all'), ! training (cv, 'all'));
%! assert (size (test (cv, 'all')), [12, 2]);
%! assert (sum (test (cv, 'all')), [7, 5]);
%! assert (sum (training (cv, 'all')), [5, 7]);
%! assert (test (cv, 1)', grpvar == 2);
%! assert (test (cv, 2)', grpvar != 2);
%!test
%! grpvar = [1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3];
%! rand ('seed', 5);
%! cv = cvpartition (12, 'kfold', 2, 'GroupingVariables', grpvar);
%! assert (cv.Type, 'kfold');
%! assert (cv.NumObservations, 12);
%! assert (cv.NumTestSets, 2);
%! assert (cv.TrainSize, [7, 5]);
%! assert (cv.TestSize, [5, 7]);
%! assert (cv.IsCustom, false);
%! assert (cv.IsGrouped, true);
%! assert (cv.IsStratified, false);
%! assert (test (cv, 1), ! training (cv, 1));
%! assert (test (cv, 'all'), ! training (cv, 'all'));
%! assert (size (test (cv, 'all')), [12, 2]);
%! assert (sum (test (cv, 'all')), [5, 7]);
%! assert (sum (training (cv, 'all')), [7, 5]);
%! assert (test (cv, 1)', grpvar == 2);
%! assert (test (cv, 2)', grpvar != 2);
%!test
%! grpvar = [1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3];
%! rand ('seed', 5);
%! cv = cvpartition (12, 'kfold', 2, 'GroupingVariables', grpvar);
%! assert (cv.Type, 'kfold');
%! assert (cv.NumObservations, 12);
%! assert (cv.NumTestSets, 2);
%! assert (cv.TrainSize, [7, 5]);
%! assert (cv.TestSize, [5, 7]);
%! assert (cv.IsCustom, false);
%! assert (cv.IsGrouped, true);
%! assert (cv.IsStratified, false);
%! assert (test (cv, 1), ! training (cv, 1));
%! assert (test (cv, 'all'), ! training (cv, 'all'));
%! assert (size (test (cv, 'all')), [12, 2]);
%! assert (sum (test (cv, 'all')), [5, 7]);
%! assert (sum (training (cv, 'all')), [7, 5]);
%! assert (test (cv, 1)', grpvar == 3);
%! assert (test (cv, 2)', grpvar != 3);
%!test
%! status = warning;
%! warning ('off');
%! cv = cvpartition (5, 'kfold', 5, 'GroupingVariables', {'a';'a';'b';'b';''});
%! warning (status);
%! assert (cv.Type, 'kfold');
%! assert (cv.NumObservations, 5);
%! assert (cv.NumTestSets, 2);
%! assert (cv.TrainSize, [2, 2]);
%! assert (cv.TestSize, [2, 2]);
%! assert (cv.IsCustom, false);
%! assert (cv.IsGrouped, true);
%! assert (cv.IsStratified, false);
%! idx = ! ismissing ({'a';'a';'b';'b';''});
%! assert (test (cv, 1)(idx), ! training (cv, 1)(idx));
%! assert (test (cv, 'all')(idx,:), ! training (cv, 'all')(idx,:));
%! assert (size (test (cv, 'all')), [5, 2]);
%! assert (sum (test (cv, 'all')), [2, 2]);
%! assert (sum (test (cv, 'all'), 2), [1; 1; 1; 1; 0]);

## Test output results for vector input X
%!test
%! rand ('seed', 5);
%! cv = cvpartition ([1, 1, 1, 1, 1, 2, 2, 2, 2, 2], 'holdout', 3);
%! assert (cv.Type, 'holdout');
%! assert (cv.NumObservations, 10);
%! assert (cv.NumTestSets, 1);
%! assert (cv.TrainSize, 7);
%! assert (cv.TestSize, 3);
%! assert (cv.IsCustom, false);
%! assert (cv.IsGrouped, false);
%! assert (cv.IsStratified, true);
%! assert (test (cv, 1), ! training (cv, 1));
%! assert (test (cv), logical ([0, 0, 0, 0, 1, 0, 1, 0, 0, 1])');
%!test
%! cv = cvpartition ([1, 1, 1, 1, 1, 2, 2, 2, 2, 2], 'holdout', 4);
%! assert (cv.Type, 'holdout');
%! assert (cv.NumObservations, 10);
%! assert (cv.NumTestSets, 1);
%! assert (cv.TrainSize, 6);
%! assert (cv.TestSize, 4);
%! assert (cv.IsCustom, false);
%! assert (cv.IsGrouped, false);
%! assert (cv.IsStratified, true);
%! assert (test (cv, 1), ! training (cv, 1));
%! assert (sum (test (cv)(1:5)), 2);
%! assert (sum (test (cv)(6:10)), 2);
%!test
%! grpvar = [1, 1, 1, 1, 1, 2, 2, 2, 2, 2];
%! rand ('seed', 5);
%! cv = cvpartition (grpvar, 'holdout', 4, 'Stratify', false);
%! assert (cv.Type, 'holdout');
%! assert (cv.NumObservations, 10);
%! assert (cv.NumTestSets, 1);
%! assert (cv.TrainSize, 6);
%! assert (cv.TestSize, 4);
%! assert (cv.IsCustom, false);
%! assert (cv.IsGrouped, false);
%! assert (cv.IsStratified, false);
%! assert (test (cv, 1), ! training (cv, 1));
%! assert (sum (test (cv)(1:5)), 3);
%! assert (sum (test (cv)(6:10)), 1);
%!test
%! cv = cvpartition ([1 1 1 1 1 2 2 2 2 1], 'kfold', 2);
%! assert (cv.Type, 'kfold');
%! assert (cv.NumObservations, 10);
%! assert (cv.NumTestSets, 2);
%! assert (cv.TrainSize, [5, 5]);
%! assert (cv.TestSize, [5, 5]);
%! assert (cv.IsCustom, false);
%! assert (cv.IsGrouped, false);
%! assert (cv.IsStratified, true);
%! assert (test (cv, 1), ! training (cv, 1));
%! assert (test (cv, 'all'), ! training (cv, 'all'));
%! assert (sum (test (cv, 1)(1:5)), 3);
%! assert (sum (test (cv, 2)(1:5)), 2);
%! assert (sum (test (cv, 1)(6:10)), 2);
%! assert (sum (test (cv, 2)(6:10)), 3);
%!test
%! grpvar = [1 1 1 1 1 2 2 2 2 1];
%! rand ('seed', 5);
%! cv = cvpartition (grpvar, 'kfold', 2, 'Stratify', false);
%! assert (cv.Type, 'kfold');
%! assert (cv.NumObservations, 10);
%! assert (cv.NumTestSets, 2);
%! assert (cv.TrainSize, [5, 5]);
%! assert (cv.TestSize, [5, 5]);
%! assert (cv.IsCustom, false);
%! assert (cv.IsGrouped, false);
%! assert (cv.IsStratified, false);
%! assert (test (cv, 1), ! training (cv, 1));
%! assert (test (cv, 'all'), ! training (cv, 'all'));
%! assert (sum (test (cv, 1)(1:5)), 4);
%! assert (sum (test (cv, 2)(1:5)), 1);
%! assert (sum (test (cv, 1)(6:10)), 1);
%! assert (sum (test (cv, 2)(6:10)), 4);
%!test
%! status = warning;
%! warning ('off');
%! cv = cvpartition ({'a','a','b','b',''}, 'kfold');
%! warning (status);
%! assert (cv.Type, 'kfold');
%! assert (cv.NumObservations, 5);
%! assert (cv.NumTestSets, 4);
%! assert (cv.TrainSize, [3, 3, 3, 3]);
%! assert (cv.TestSize, [1, 1, 1, 1]);
%! assert (cv.IsCustom, false);
%! assert (cv.IsGrouped, false);
%! assert (cv.IsStratified, true);
%! idx = ! ismissing ({'a','a','b','b',''});
%! assert (test (cv, 1)(idx), ! training (cv, 1)(idx));
%! assert (test (cv, 'all')(idx,:), ! training (cv, 'all')(idx,:));
%! assert (sum (test (cv, 'all'), 2), [1; 1; 1; 1; 0]);

## Test input validation
%!error <cvpartition: too few input arguments.> cvpartition (2)
%!error <cvpartition: too many input arguments.> cvpartition (1, 2, 3, 4, 5, 6)
%!error <cvpartition: TESTSETS must be numeric of logical.> ...
%! cvpartition ("CustomPartition", 'a')
%!error <cvpartition: TESTSETS must be a numeric vector.> ...
%! cvpartition ("CustomPartition", [2, 3; 2, 3])
%!error <cvpartition: TESTSETS must be a logical vector or matrix.> ...
%! cvpartition ("CustomPartition", false (3, 3, 3))
%!error <cvpartition: each observation in TESTSETS must be exactly one in each row.> ...
%! cvpartition ("CustomPartition", [false, true; true, true; true, false])
%!error <cvpartition: a logical matrix in TESTSETS must not have more columns that rows.> ...
%! cvpartition ("CustomPartition", false (3, 5))
%!error <cvpartition: X must be a scalar positive integer value.> ...
%! cvpartition (-20, "LeaveOut")
%!error <cvpartition: X must be a scalar positive integer value.> ...
%! cvpartition (20.5, "LeaveOut")
%!error <cvpartition: P value for 'holdout' must be a numeric scalar.> ...
%! cvpartition (20, "HoldOut", [0.2, 0.3])
%!error <cvpartition: P value for 'holdout' must be a numeric scalar.> ...
%! cvpartition (20, "HoldOut", 'a')
%!error <cvpartition: P value for 'holdout' must be a scalar in the range \(0,1\) or an integer scalar in the range \[1, N\).> ...
%! cvpartition (20, "HoldOut", 0)
%!error <cvpartition: P value for 'holdout' must be a scalar in the range \(0,1\) or an integer scalar in the range \[1, N\).> ...
%! cvpartition (20, "HoldOut", -0.1)
%!error <cvpartition: P value for 'holdout' must be a scalar in the range \(0,1\) or an integer scalar in the range \[1, N\).> ...
%! cvpartition (20, "HoldOut", 21)
%!error <cvpartition: K value for 'kfold' must be a numeric scalar.> ...
%! cvpartition (20, "kfold", [2, 3])
%!error <cvpartition: K value for 'kfold' must be a numeric scalar.> ...
%! cvpartition (20, "kfold", 'a')
%!error <cvpartition: K value for 'kfold' must be an integer scalar in the range \[1, N\].> ...
%! cvpartition (20, "kfold", 2.5)
%!error <cvpartition: K value for 'kfold' must be an integer scalar in the range \[1, N\].> ...
%! cvpartition (20, "kfold", 21)
%!error <cvpartition: invalid optional paired argument for 'GroupingVariables'.> ...
%! cvpartition (10, "kfold", 3, "Group")
%!error <cvpartition: missing value for optional paired argument 'GroupingVariables'.> ...
%! cvpartition (10, "kfold", 3, "GroupingVariables")
%!error <cvpartition: invalid value for optional paired argument 'GroupingVariables'.> ...
%! cvpartition (10, "kfold", 3, "GroupingVariables", ones (3, 3, 3))
%!error <cvpartition: grouping variable does not match the number of observations.> ...
%! cvpartition (10, "kfold", 3, "GroupingVariables", {'a', 'a', 'a', 'b', 'b'})
%!warning <cvpartition: number of folds K is greater than the groups in 'GroupingVariables'. K is set to the number of groups.> ...
%! cvpartition (5, "kfold", 3, "GroupingVariables", {'a', 'a', 'a', 'b', 'b'});
%!error <cvpartition: invalid optional paired argument.> ...
%! cvpartition (20, "some")
%!error <cvpartition: invalid optional paired argument for stratification.> ...
%! cvpartition ([1, 1, 1, 2, 2], "kfold", 2, "strat")
%!error <cvpartition: missing value for optional paired argument 'stratify'.> ...
%! cvpartition ([1, 1, 1, 2, 2], "kfold", 2, "stratify")
%!error <cvpartition: invalid value for optional paired argument 'stratify'.> ...
%! cvpartition ([1, 1, 1, 2, 2], "kfold", 2, "stratify", [true, true])
%!error <cvpartition: invalid value for optional paired argument 'stratify'.> ...
%! cvpartition ([1, 1, 1, 2, 2], "kfold", 2, "stratify", 'no')
%!error <cvpartition: P value for 'holdout' must be a numeric scalar.> ...
%! cvpartition ([1, 1, 1, 2, 2], "holdout", 'a')
%!error <cvpartition: P value for 'holdout' must be a numeric scalar.> ...
%! cvpartition ([1, 1, 1, 2, 2], "holdout", 'a', "stratify", true)
%!error <cvpartition: P value for 'holdout' must be a numeric scalar.> ...
%! cvpartition ([1, 1, 1, 2, 2], "holdout", [0.2, 0.3])
%!error <cvpartition: P value for 'holdout' must be a numeric scalar.> ...
%! cvpartition ([1, 1, 1, 2, 2], "holdout", [0.2, 0.3], "stratify", true)
%!error <cvpartition: P value for 'holdout' must be a scalar in the range \(0,1\) or an integer scalar in the range \[1, N\), where N is the number of nonmissing observations in X.> ...
%! cvpartition ([1, 1, 1, 2, 2], "holdout", 0)
%!error <cvpartition: P value for 'holdout' must be a scalar in the range \(0,1\) or an integer scalar in the range \[1, N\), where N is the number of nonmissing observations in X.> ...
%! cvpartition ([1, 1, 1, 2, 2], "holdout", 0, "stratify", true)
%!error <cvpartition: P value for 'holdout' must be a scalar in the range \(0,1\) or an integer scalar in the range \[1, N\), where N is the number of nonmissing observations in X.> ...
%! cvpartition ([1, 1, 1, 2, 2], "holdout", -0.1)
%!error <cvpartition: P value for 'holdout' must be a scalar in the range \(0,1\) or an integer scalar in the range \[1, N\), where N is the number of nonmissing observations in X.> ...
%! cvpartition ([1, 1, 1, 2, 2], "holdout", -0.1, "stratify", true)
%!error <cvpartition: P value for 'holdout' must be a scalar in the range \(0,1\) or an integer scalar in the range \[1, N\), where N is the number of nonmissing observations in X.> ...
%! cvpartition ([1, 1, 1, 2, 2], "holdout", 1.2)
%!error <cvpartition: P value for 'holdout' must be a scalar in the range \(0,1\) or an integer scalar in the range \[1, N\), where N is the number of nonmissing observations in X.> ...
%! cvpartition ([1, 1, 1, 2, 2], "holdout", 1.2, "stratify", false)
%!error <cvpartition: P value for 'holdout' must be a scalar in the range \(0,1\) or an integer scalar in the range \[1, N\), where N is the number of nonmissing observations in X.> ...
%! cvpartition ([1, 1, 1, 2, 2], "holdout", 6)
%!error <cvpartition: P value for 'holdout' must be a scalar in the range \(0,1\) or an integer scalar in the range \[1, N\), where N is the number of nonmissing observations in X.> ...
%! cvpartition ([1, 1, 1, 2, 2], "holdout", 6, "stratify", false)
%!error <cvpartition: K value for 'kfold' must be a numeric scalar.> ...
%! cvpartition ([1, 1, 1, 2, 2], "kfold", 'a')
%!error <cvpartition: K value for 'kfold' must be a numeric scalar.> ...
%! cvpartition ([1, 1, 1, 2, 2], "kfold", 'a', "stratify", true)
%!error <cvpartition: K value for 'kfold' must be a numeric scalar.> ...
%! cvpartition ([1, 1, 1, 2, 2], "kfold", [2, 3])
%!error <cvpartition: K value for 'kfold' must be a numeric scalar.> ...
%! cvpartition ([1, 1, 1, 2, 2], "kfold", [2, 3], "stratify", false)
%!error <cvpartition: K value for 'kfold' must be an integer scalar in the range \[1, N\], where N is the number of nonmissing observations in X.> ...
%! cvpartition ([1, 1, 1, 2, 2], "kfold", 0)
%!error <cvpartition: K value for 'kfold' must be an integer scalar in the range \[1, N\], where N is the number of nonmissing observations in X.> ...
%! cvpartition ([1, 1, 1, 2, 2], "kfold", 0, "stratify", true)
%!error <cvpartition: K value for 'kfold' must be an integer scalar in the range \[1, N\], where N is the number of nonmissing observations in X.> ...
%! cvpartition ([1, 1, 1, 2, 2], "kfold", 1.5)
%!error <cvpartition: K value for 'kfold' must be an integer scalar in the range \[1, N\], where N is the number of nonmissing observations in X.> ...
%! cvpartition ([1, 1, 1, 2, 2], "kfold", 1.5, "stratify", true)
%!error <cvpartition: K value for 'kfold' must be an integer scalar in the range \[1, N\], where N is the number of nonmissing observations in X.> ...
%! cvpartition ([1, 1, 1, 2, 2], "kfold", 6)
%!error <cvpartition: K value for 'kfold' must be an integer scalar in the range \[1, N\], where N is the number of nonmissing observations in X.> ...
%! cvpartition ([1, 1, 1, 2, 2], "kfold", 6, "stratify", true)
%!error <cvpartition: invalid optional paired argument.> ...
%! cvpartition ([1, 1, 1, 2, 2], "leaveout")
%!error <cvpartition: invalid optional paired argument.> ...
%! cvpartition ([1, 1, 1, 2, 2], "resubstitution")
%!error <cvpartition: invalid optional paired argument.> ...
%! cvpartition ([1, 1, 1, 2, 2], "some")
%!error <cvpartition: invalid first input argument.> ...
%! cvpartition ({1, 1; 2, 2}, "kfold")

%!error <cvpartition.repartition: cannot repartition a custom partition.> ...
%! repartition (cvpartition ('CustomPartition', [1,1,2,2,3,3]))
%!error <cvpartition.repartition: 'legacy' flag does not apply to stratified or grouped 'kfold' partitioned objects.> ...
%! repartition (cvpartition ([1 1 1 1 1 2 2 2 2 1], 'kfold', 2, 'Stratify', true), 'legacy')
%!error <cvpartition.repartition: 'legacy' flag is only valid for 'kfold' partitioned objects.> ...
%! repartition (cvpartition (20, 'Leaveout', 0.2), 'legacy')
%!error <cvpartition.repartition: SVAL must be a real scalar or vector.> ...
%! repartition (cvpartition (20, 'Leaveout', 0.2), 'asd')
%!error <cvpartition.repartition: SVAL must be a real scalar or vector.> ...
%! repartition (cvpartition (20, 'Leaveout', 0.2), 2+i)
%!error <cvpartition.repartition: SVAL must be a real scalar or vector.> ...
%! repartition (cvpartition (20, 'KFold', 5), [34, 56; 2, 3])

%!error <cvpartition.test: too many input arguments.> ...
%! test (cvpartition (20, "kfold"), 2, 3)
%!error <cvpartition.test: set index must be a positive integer vector.> ...
%! test (cvpartition (20, "kfold"), 0)
%!error <cvpartition.test: set index must be a positive integer vector.> ...
%! test (cvpartition (20, "kfold"), 1.5)
%!error <cvpartition.test: set index must be a positive integer vector.> ...
%! test (cvpartition (20, "kfold"), [1, 1.5])
%!error <cvpartition.test: set index must be a positive integer vector.> ...
%! test (cvpartition (20, "kfold"), [2, 3; 2, 3])
%!error <cvpartition.test: set index exceeds 'NumTestSets'.> ...
%! test (cvpartition (20, "kfold"), 21)
%!error <cvpartition.test: set index exceeds 'NumTestSets'.> ...
%! test (cvpartition (20, "kfold"), [18, 21])

%!error <cvpartition.training: too many input arguments.> ...
%! training (cvpartition (20, "kfold"), 2, 3)
%!error <cvpartition.training: set index must be a positive integer vector.> ...
%! training (cvpartition (20, "kfold"), 0)
%!error <cvpartition.training: set index must be a positive integer vector.> ...
%! training (cvpartition (20, "kfold"), 1.5)
%!error <cvpartition.training: set index must be a positive integer vector.> ...
%! training (cvpartition (20, "kfold"), [1, 1.5])
%!error <cvpartition.training: set index must be a positive integer vector.> ...
%! training (cvpartition (20, "kfold"), [2, 3; 2, 3])
%!error <cvpartition.training: set index exceeds 'NumTestSets'.> ...
%! training (cvpartition (20, "kfold"), 21)
%!error <cvpartition.training: set index exceeds 'NumTestSets'.> ...
%! training (cvpartition (20, "kfold"), [18, 21])
