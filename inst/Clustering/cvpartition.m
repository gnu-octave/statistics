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
## @deftp {Class} cvpartition
##
## Partition data for cross-validation
##
## The @code{cvpartition} class generates a partitioning scheme on a dataset to
## facilitate cross-validation of statistical models utilizing training and
## testing subsets of the dataset.
##
## @seealso{crossval}
## @end deftp

  properties (SetAccess = private)
    ## -*- texinfo -*-
    ## @deftp {cvpartition} {property} NumObservations
    ##
    ## A positive integer scalar specifying the number of observations in the
    ## dataset (including any missing data, where applicable).  This property
    ## is read-only.
    ##
    ## @end deftp
    NumObservations = [];     # Number of observations

    ## @deftp {cvpartition} {property} NumTestSets
    ##
    ## A positive integer scalar specifying the number of folds for partition
    ## types @qcode{'kfold'} and @qcode{'leaveout'}.  When partition type is
    ## @qcode{'holdout'} and @qcode{'resubstitution'}, then @qcode{NumTestSets}
    ## is 1.  This property is read-only.
    ##
    ## @end deftp
    NumTestSets     = [];     # Number of test sets

    ## @deftp {cvpartition} {property} TrainSize
    ##
    ## A positive integer scalar specifying the size of the train set for
    ## partition types @qcode{'holdout'} and @qcode{'resubstitution'} or a
    ## vector of positive integers specifying the size of each training set for
    ## partition types @qcode{'kfold'} and @qcode{'leaveout'}.  This property
    ## is read-only.
    ##
    ## @end deftp
    TrainSize       = [];     # Size of each train set

    ## @deftp {cvpartition} {property} TestSize
    ##
    ## A positive integer scalar specifying the size of the test set for
    ## partition types @qcode{'holdout'} and @qcode{'resubstitution'} or a
    ## vector of positive integers specifying the size of each testing set for
    ## partition types @qcode{'kfold'} and @qcode{'leaveout'}.  This property
    ## is read-only.
    ##
    ## @end deftp
    TestSize        = [];     # Size of each test set

    ## @deftp {cvpartition} {property} Type
    ##
    ## A character vector specifying the type of the @qcode{cvpartition} object.
    ## It can be @qcode{kfold}, @qcode{holdout}, @qcode{leaveout}, or
    ## @qcode{resubstitution}.  This property is read-only.
    ##
    ## @end deftp
    Type            = '';     # Type of validation partition

    ## @deftp {cvpartition} {property} IsCustom
    ##
    ## A logical scalar specifying whether the @qcode{cvpartition} object
    ## was created using custom partition partitioning (@qcode{true}) or
    ## not (@qcode{false}).  This property is read-only.
    ##
    ## @end deftp
    IsCustom        = [];     # Indicator of custom partition

    ## @deftp {cvpartition} {property} IsGrouped
    ##
    ## A logical scalar specifying whether the @qcode{cvpartition} object was
    ## created using grouping variables (@qcode{true}) or not (@qcode{false}).
    ## This property is read-only.
    ##
    ## @end deftp
    IsGrouped       = [];     # Indicator of grouped partition

    ## @deftp {cvpartition} {property} IsStratified
    ##
    ## A logical scalar specifying whether the @qcode{cvpartition} object was
    ## created with a @qcode{'stratifyOption'} value of @qcode{true}.
    ## This property is read-only.
    ##
    ## @end deftp
    IsStratified    = [];     # Indicator of stratified partition

  endproperties

  properties (SetAccess = private, Hidden)
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
            error ("cvpartition.subref: unrecongized property: '%s'", s.subs);
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
          error (strcat ("cvpartition.subsasgn: unrecongized", ...
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
    ## @code{@var{C} =} cvpartition (@qcode{'CustomPartition'}, @var{testSets})}
    ## creates a custom partition according to @var{testSets}, which can be a
    ## positive integer vector, a logical vector, or a logical matrix according
    ## to the following options:
    ## @itemize
    ## @item A positive integer vector of length @math{n} with values in the
    ## range @math{[1,k]}, where @math{k < n}, will specify a k-fold
    ## cross-validation partition, in which each value indicates the test set
    ## of each observation.  Alternatively, the same vector with values in the
    ## range @math{[1,n]} will specify a leave-one-out cross-validation.
    ## @item A logical vector will specify a holdout validation, in which the
    ## @qcode{true} elements correspond to the test set and the @qcode{false}
    ## elements correspond to the traning set.
    ## @item A logical matrix with @math{k} columns will specify a k-fold
    ## cross-validation partition, in which each collumn corresponds to a fold
    ## and each row to an observation.  Alternatively, an @math{nxn} logical
    ## matrix will specify a leave-one-out cross-validation, where @math{n} is
    ## the number of observations.  @qcode{true} elements correspond to the
    ## test set and the @qcode{false} elements correspond to the traning set.
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
          error ("cvpartition: testSets must be numeric of logical.");
        endif
        if (isnumeric (testSets))
          if (! isvector (testSets))
            error ("cvpartition: testSets must be a numeric vector.");
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
            error ("cvpartition: testSets must be a logical vector or matrix.");
          elseif (isvector (testSets))
            this.NumObservations = numel (testSets);
            this.NumTestSets = 1;
            this.indices = testSets;
            this.TrainSize = sum (! testSets);
            this.TestSize = sum (testSets);
            this.Type = 'holdout';
            this.cvptype = 'Hold-out cross validation partition';
          else  # logical matrix
            [~, idx, inds] = unique (testSets, 'rows');
            this.NumObservations = size (testSets, 1);
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
        ## Get partition type
        type = varargin{1};
        this.IsCustom = false;
        this.IsStratified = false;

        ## "Resubstitution"
        if (strcmpi (type, 'resubstitution'))
          this.NumObservations = X;
          this.NumTestSets = 1;
          this.TrainSize = X;
          this.TestSize = X;
          this.Type = 'resubstitution';
          this.cvptype = 'Resubstitution (no partition of data)';
          this.IsGrouped = false;

        ## "Leaveout"
        elseif (strcmpi (type, 'leaveout'))
          this.NumObservations = X;
          this.NumTestSets = 1;
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
              error (strcat ("cvpartition: p value for 'holdout'", ...
                             " must be a numeric scalar."));
            endif
            if (! ((p > 0 && p < 1) || (p == fix (p) && p > 0 && p < X)))
              error (strcat ("cvpartition: p value for 'holdout' must be", ...
                             " a scalar in the range (0,1) or an integer", ...
                             " scalar in the range [1, n)."));
            endif
          else
            p = 0.1;
          endif
          this.NumObservations = X;
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
              error (strcat ("cvpartition: k value for 'kfold'", ...
                             " must be a numeric scalar."));
            endif
            if (! (k == fix (k) && k > 0 && k < X))
              error (strcat ("cvpartition: k value for 'kfold' must be", ...
                             " an integer scalar in the range [1, n)."));
            endif
          else
            if (X <= 10)
              k = X - 1;
            else
              k = 10;
            endif
          endif
          ## No grouping variables
          if (nargin < 4)
            this.NumObservations = X;
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
              if (isa (grpvars, 'categorical'))
                [~, idx, inds] = unique (grpvars, 'stable');
              else
                [~, idx, inds] = __unique__ (grpvars, 'stable');
              endif
            elseif (ismatrix (grpvars))
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
            ## Each k-fold attemps to spit the groups to equal sizes in such a
            ## way so that eash test set contains unique groups that are not
            ## present in the corresponding training set but also not shared
            ## with other test sets.
            GroupIdx = multiway (GroupSize, k);
            GroupIdx = randsample (GroupIdx, k);
            indices = zeros (X, 1);
            for i = 1:k
              idxGV = inds(idx(GroupIdx == i));
              vecGV = arrayfun(@(x) x == inds, idxGV, "UniformOutput", false);
              index = vecGV{1};
              if (numel (vecGV) > 1)
                for j = 2:numel (vecGV)
                  index = index | vecGV{j};
                endfor
              endif
              indices(index) = i;
            endfor
            this.indices = indices;
            this.NumObservations = X;
            this.NumTestSets = k;
            nvec = X * ones (1, k);
            for i = 1:k
              this.TestSize(i) = sum (this.indices == i);
            endfor
            this.TrainSize = nvec - this.TestSize;
            this.cvptype = 'Group k-fold cross validation partition';
            this.IsGrouped = true;
          endif

        ## Invalid paired argument
        else
          error ("cvpartition: invalid optional paired argument.");
        endif

      ## Check first input being a vector for stratification
      elseif (isvector (X))
        ## Remove missing values
        X(ismissing (X)) = [];
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
        else
          this.cvptype = 'Hold-out cross validation partition';
        endif

        ## Get number of observations
        X = numel (X);
        this.NumObservations = X;

        ## Get partition type
        type = varargin{1};
        this.IsCustom = false;
        this.IsGrouped = false;

        ## "Holdout"
        if (strcmpi (type, 'holdout'))
          if (nargin > 2)
            p = varargin{2};
            if (! isnumeric (p) || ! isscalar (p))
              error (strcat ("cvpartition: p value for 'holdout'", ...
                             " must be a numeric scalar."));
            endif
            if (! ((p > 0 && p < 1) || (p == fix (p) && p > 0 && p < X)))
              error (strcat ("cvpartition: p value for 'holdout' must be", ...
                             " a scalar in the range (0,1) or an integer", ...
                             " scalar in the range [1, n)."));
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
          else
            if (p < 1)            # target fraction to sample
              p = round (p * X);  # number of samples
            endif
            inds = false (X, 1);
            inds(randsample (X, p)) = true;  # indices for test set
          endif
          this.indices = inds;
          this.TrainSize = sum (! inds);
          this.TestSize = sum (inds);
          this.Type = 'holdout';
          this.cvptype = 'Stratified hold-out cross validation partition';

        ## "KFold"
        elseif (strcmpi (type, 'kfold'))
          if (nargin > 2)
            k = varargin{2};
            if (! isnumeric (k) || ! isscalar (k))
              error (strcat ("cvpartition: k value for 'kfold'", ...
                             " must be a numeric scalar."));
            endif
            if (! (k == fix (k) && k > 0 && k < X))
              error (strcat ("cvpartition: k value for 'kfold' must be", ...
                             " an integer scalar in the range [1, n)."));
            endif
          else
            if (X <= 10)
              k = X - 1;
            else
              k = 10;
            endif
          endif
          this.NumTestSets = k;
          if (this.IsStratified)
            inds = nan (X, 1);
            for i = 1:NumClasses
              ## Alternate ordering over classes so that
              ## the subsets are more nearly the same size
              if (mod (i, 2))
                idx = floor ((0:(ClassSize(i)-1))' * (k / ClassSize(i))) + 1;
                inds(classes == i) = randsample (idx, ClassSize(i));
              else
                idx = floor (((ClassSize(i)-1):-1:0)' * (k / ClassSize(i))) + 1;
                inds(classes == i) = randsample (idx, ClassSize(i));
              endif
            endfor
          else
            inds = floor ((0:(X - 1))' * (k / X)) + 1;
            inds = randsample (inds, X);
          endif
          this.indices = inds;
          nvec = X * ones (1, k);
          for i = 1:k
            this.TestSize(i) = sum (inds == i);
          endfor
          this.TrainSize = nvec - this.TestSize;
          this.Type = 'kfold';
          this.cvptype = 'Stratified k-fold cross validation partition';

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
    ## @deftypefn {cvpartition} {@var{Cnew} =} repartition (@var{C})
    ## @deftypefn {cvpartition} {@var{Cnew} =} repartition (@var{C}, @var{sval})
    ## @deftypefn {cvpartition} {@var{Cnew} =} repartition (@var{C}, @qcode{'legacy'})
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
    ## statistics package.
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
          elseif (this.IsGrouped)
            grpvars = this.grpvars;
            if (isvector (grpvars))
              if (isa (grpvars, 'categorical'))
                [~, idx, inds] = unique (grpvars, 'stable');
              else
                [~, idx, inds] = __unique__ (grpvars, 'stable');
              endif
            else
              if (isa (grpvars, 'categorical'))
                [~, idx, inds] = unique (grpvars, 'rows', 'stable');
              else
                [~, idx, inds] = __unique__ (grpvars, 'rows', 'stable');
              endif
            endif
            ## Get number of groups and group sizes
            NumGroups = numel (idx);
            for i = 1:NumGroups
              GroupSize(i) = sum (inds == i);
            endfor
            ## Each k-fold attemps to spit the groups to equal sizes in such a
            ## way so that eash test set contains unique groups that are not
            ## present in the corresponding training set but also not shared
            ## with other test sets.
            GroupIdx = multiway (GroupSize, k);
            indices = zeros (X, 1);
            for i = 1:k
              idxGV = inds(idx(GroupIdx == i));
              vecGV = arrayfun(@(x) x == inds, idxGV, "UniformOutput", false);
              index = vecGV{1};
              if (numel (vecGV) > 1)
                for j = 2:numel (vecGV)
                  index = index | vecGV{j};
                endfor
              endif
              indices(index) = i;
            endfor
            this.indices = indices;
            nvec = X * ones (1, k);
            for i = 1:k
              this.TestSize(i) = sum (this.indices == i);
            endfor
            this.TrainSize = nvec - this.TestSize;
          else  # is stratified
            NumClasses = numel (this.classID);
            classes = this.classes;
            for i = 1:NumClasses
              ClassSize(i) = sum (classes == i);
            endfor
            inds = nan (X, 1);
            for i = 1:NumClasses
              ## Alternate ordering over classes so that
              ## the subsets are more nearly the same size
              if (mod (i, 2))
                inds(classes == i) = floor ((0:(ClassSize(i)-1))' * ...
                                            (k / ClassSize(i))) + 1;
              else
                inds(classes == i) = floor (((ClassSize(i)-1):-1:0)' * ...
                                            (k / ClassSize(i))) + 1;
              endif
            endfor
            this.indices = inds;
            nvec = X * ones (1, k);
            for i = 1:k
              this.TestSize(i) = sum (inds == i);
            endfor
            this.TrainSize = nvec - this.TestSize;
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
        X = this.NumObservations;
        p = this.TestSize;
        inds = false (X, 1);
        if (this.IsStratified)
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
          inds(randsample (X, p)) = true;  # indices for test set
        endif
        this.indices = inds;
      elseif (strcmpi (this.Type, "kfold"))
        X = this.NumObservations;
        k = this.NumTestSets;
        if (! (this.IsGrouped || this.IsStratified))
          inds = floor ((0:(X - 1))' * (k / X)) + 1;
          inds = randsample (inds, X);
          this.indices = inds;
          nvec = X * ones (1, k);
          for i = 1:k
            this.TestSize(i) = sum (inds == i);
          endfor
          this.TrainSize = nvec - this.TestSize;
        elseif (this.IsGrouped)
          grpvars = this.grpvars;
          if (isvector (grpvars))
            if (isa (grpvars, 'categorical'))
              [~, idx, inds] = unique (grpvars, 'stable');
            else
              [~, idx, inds] = __unique__ (grpvars, 'stable');
            endif
          else
            if (isa (grpvars, 'categorical'))
              [~, idx, inds] = unique (grpvars, 'rows', 'stable');
            else
              [~, idx, inds] = __unique__ (grpvars, 'rows', 'stable');
            endif
          endif
          ## Get number of groups and group sizes
          NumGroups = numel (idx);
          for i = 1:NumGroups
            GroupSize(i) = sum (inds == i);
          endfor
          ## Each k-fold attemps to spit the groups to equal sizes in such a
          ## way so that eash test set contains unique groups that are not
          ## present in the corresponding training set but also not shared
          ## with other test sets.
          GroupIdx = multiway (GroupSize, k);
          GroupIdx = randsample (GroupIdx, k);
          indices = zeros (X, 1);
          for i = 1:k
            idxGV = inds(idx(GroupIdx == i));
            vecGV = arrayfun(@(x) x == inds, idxGV, "UniformOutput", false);
            index = vecGV{1};
            if (numel (vecGV) > 1)
              for j = 2:numel (vecGV)
                index = index | vecGV{j};
              endfor
            endif
            indices(index) = i;
          endfor
          this.indices = indices;
          nvec = X * ones (1, k);
          for i = 1:k
            this.TestSize(i) = sum (this.indices == i);
          endfor
          this.TrainSize = nvec - this.TestSize;
        else  # is stratified
          NumClasses = numel (this.classID);
          classes = this.classes;
          for i = 1:NumClasses
            ClassSize(i) = sum (classes == i);
          endfor
          inds = nan (X, 1);
          for i = 1:NumClasses
            ## Alternate ordering over classes so that
            ## the subsets are more nearly the same size
            if (mod (i, 2))
              idx = floor ((0:(ClassSize(i)-1))' * (k / ClassSize(i))) + 1;
              inds(classes == i) = randsample (idx, ClassSize(i));
            else
              idx = floor (((ClassSize(i)-1):-1:0)' * (k / ClassSize(i))) + 1;
              inds(classes == i) = randsample (idx, ClassSize(i));
            endif
          endfor
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
    ## set defined in the code{cvpartition} object @var{C}.  For k-fold and
    ## leave-one-out partitions, the indices corresponding to the first test set
    ## are returned.
    ##
    ## @code{@var{idx} = test (@var{C}, @var{i})} returns a logical vector or
    ## matrix with the indices of the test set indicated by @var{i}.  If @var{i}
    ## is a scalar, then @var{idx} is a logical vector with the indices of the
    ## @math{i-th} set.  If @var{i} is a vector, then @var{idx} is a logical
    ## matrix in which @code{@var{idx}(:,j)} specified the observations in the
    ## test set @code{@var{i}(j)}.  The value(s) in @var{i} must not excced the
    ## number of tests in the @qcode{cvpartition} object @var{C}.
    ##
    ## @code{@var{idx} = test (@var{C}, @qcode{"all"})} returns a logical vector
    ## or matrix for all test sets defined in the @qcode{cvpartition} object
    ## @var{C}.  For holdout and resubstitution partition types, a vector is
    ## returned.  For k-fold and leave-one-out, a matrix is returned.
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
          idx = [];
          switch (this.Type)
            case "kfold"
              for i = 1:this.NumTestSets
                cid = this.indices == i;
                idx = [idx, cid];
              endfor
            case "leaveout"
              for i = 1:this.NumTestSets
                cid = false (this.NumObservations, 1);
                cid(i) = true;
                idx = [idx, cid];
              endfor
            case "holdout"
              idx = this.indices;
            case "resubstitution"
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
            idx = this.indices == i;
          else
            idx = [];
            for j = i
              new = this.indices == j;
              idx = [idx, new];
            endfor
          endif
        case "leaveout"
          if (isscalar (i))
            idx = false (this.NumObservations, 1);
            idx(i) = true;
          else
            idx = [];
            for j = i
              new = false (this.NumObservations, 1);
              new(j) = true;
              idx = [idx, new];
            endfor
          endif
        case "holdout"
          idx = this.indices;
        case "resubstitution"
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
    ## training set defined in the code{cvpartition} object @var{C}.  For k-fold
    ## and leave-one-out partitions, the indices corresponding to the first
    ## training set are returned.
    ##
    ## @code{@var{idx} = training (@var{C}, @var{i})} returns a logical vector
    ## or matrix with the indices of the training set indicated by @var{i}.  If
    ## @var{i} is a scalar, then @var{idx} is a logical vector with the indices
    ## of the @math{i-th} set.  If @var{i} is a vector, then @var{idx} is a
    ## logical matrix in which @code{@var{idx}(:,j)} specified the observations
    ## in the training set @code{@var{i}(j)}.  The value(s) in @var{i} must not
    ## excced the number of tests in the @qcode{cvpartition} object @var{C}.
    ##
    ## @code{@var{idx} = training (@var{C}, @qcode{"all"})} returns a logical
    ## vector or matrix for all training sets defined in the @qcode{cvpartition}
    ## object @var{C}.  For holdout and resubstitution partition types, a vector
    ## is returned.  For k-fold and leave-one-out, a matrix is returned.
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
          idx = [];
          switch (this.Type)
            case "kfold"
              for i = 1:this.NumTestSets
                cid = this.indices != i;
                idx = [idx, cid];
              endfor
            case "leaveout"
              for i = 1:this.NumTestSets
                cid = true (this.NumObservations, 1);
                cid(i) = false;
                idx = [idx, cid];
              endfor
            case "holdout"
              idx = ! this.indices;
            case "resubstitution"
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
            idx = this.indices != i;
          else
            idx = [];
            for j = i
              new = this.indices != j;
              idx = [idx, new];
            endfor
          endif
        case "leaveout"
          if (isscalar (i))
            idx = true (this.NumObservations, 1);
            idx(i) = false;
          else
            idx = [];
            for j = i
              new = true (this.NumObservations, 1);
              new(j) = false;
              idx = [idx, new];
            endfor
          endif
        case "holdout"
          idx = ! this.indices;
        case "resubstitution"
          idx = true (this.NumObservations, 1);
      endswitch

    endfunction

  endmethods

endclassdef

## Test Input Validation
%!error <cvpartition: too few input arguments.> cvpartition (2)
%!error <cvpartition: too many input arguments.> cvpartition (1, 2, 3, 4, 5, 6)
%!error <cvpartition: testSets must be numeric of logical.> ...
%! cvpartition ("CustomPartition", 'a')
%!error <cvpartition: testSets must be a numeric vector.> ...
%! cvpartition ("CustomPartition", [2, 3; 2, 3])
%!error <cvpartition: testSets must be a logical vector or matrix.> ...
%! cvpartition ("CustomPartition", false (3, 3, 3))
%!error <cvpartition: X must be a scalar positive integer value.> ...
%! cvpartition (-20, "LeaveOut")
%!error <cvpartition: X must be a scalar positive integer value.> ...
%! cvpartition (20.5, "LeaveOut")
%!error <cvpartition: p value for 'holdout' must be a numeric scalar.> ...
%! cvpartition (20, "HoldOut", [0.2, 0.3])
%!error <cvpartition: p value for 'holdout' must be a numeric scalar.> ...
%! cvpartition (20, "HoldOut", 'a')
%!error <cvpartition: p value for 'holdout' must be a scalar in the range \(0,1\) or an integer scalar in the range \[1, n\).> ...
%! cvpartition (20, "HoldOut", 0)
%!error <cvpartition: p value for 'holdout' must be a scalar in the range \(0,1\) or an integer scalar in the range \[1, n\).> ...
%! cvpartition (20, "HoldOut", -0.1)
%!error <cvpartition: p value for 'holdout' must be a scalar in the range \(0,1\) or an integer scalar in the range \[1, n\).> ...
%! cvpartition (20, "HoldOut", 21)
%!error <cvpartition: k value for 'kfold' must be a numeric scalar.> ...
%! cvpartition (20, "kfold", [2, 3])
%!error <cvpartition: k value for 'kfold' must be a numeric scalar.> ...
%! cvpartition (20, "kfold", 'a')
%!error <cvpartition: k value for 'kfold' must be an integer scalar in the range \[1, n\).> ...
%! cvpartition (20, "kfold", 2.5)
%!error <cvpartition: k value for 'kfold' must be an integer scalar in the range \[1, n\).> ...
%! cvpartition (20, "kfold", 21)
%!error <cvpartition: invalid optional paired argument for 'GroupingVariables'.> ...
%! cvpartition (10, "kfold", 3, "Group")
%!error <cvpartition: missing value for optional paired argument 'GroupingVariables'.> ...
%! cvpartition (10, "kfold", 3, "GroupingVariables")
%!error <cvpartition: invalid value for optional paired argument 'GroupingVariables'.> ...
%! cvpartition (10, "kfold", 3, "GroupingVariables", ones (3, 3, 3))
%!error <cvpartition: grouping variable does not match the number of observations.> ...
%! cvpartition (10, "kfold", 3, "GroupingVariables", {'a', 'a', 'a', 'b', 'b'})
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
%!error <cvpartition: p value for 'holdout' must be a numeric scalar.> ...
%! cvpartition ([1, 1, 1, 2, 2], "holdout", 'a')
%!error <cvpartition: p value for 'holdout' must be a numeric scalar.> ...
%! cvpartition ([1, 1, 1, 2, 2], "holdout", 'a', "stratify", true)
%!error <cvpartition: p value for 'holdout' must be a numeric scalar.> ...
%! cvpartition ([1, 1, 1, 2, 2], "holdout", [0.2, 0.3])
%!error <cvpartition: p value for 'holdout' must be a numeric scalar.> ...
%! cvpartition ([1, 1, 1, 2, 2], "holdout", [0.2, 0.3], "stratify", true)
%!error <cvpartition: p value for 'holdout' must be a scalar in the range \(0,1\) or an integer scalar in the range \[1, n\).> ...
%! cvpartition ([1, 1, 1, 2, 2], "holdout", 0)
%!error <cvpartition: p value for 'holdout' must be a scalar in the range \(0,1\) or an integer scalar in the range \[1, n\).> ...
%! cvpartition ([1, 1, 1, 2, 2], "holdout", 0, "stratify", true)
%!error <cvpartition: p value for 'holdout' must be a scalar in the range \(0,1\) or an integer scalar in the range \[1, n\).> ...
%! cvpartition ([1, 1, 1, 2, 2], "holdout", -0.1)
%!error <cvpartition: p value for 'holdout' must be a scalar in the range \(0,1\) or an integer scalar in the range \[1, n\).> ...
%! cvpartition ([1, 1, 1, 2, 2], "holdout", -0.1, "stratify", true)
%!error <cvpartition: p value for 'holdout' must be a scalar in the range \(0,1\) or an integer scalar in the range \[1, n\).> ...
%! cvpartition ([1, 1, 1, 2, 2], "holdout", 1.2)
%!error <cvpartition: p value for 'holdout' must be a scalar in the range \(0,1\) or an integer scalar in the range \[1, n\).> ...
%! cvpartition ([1, 1, 1, 2, 2], "holdout", 1.2, "stratify", false)
%!error <cvpartition: p value for 'holdout' must be a scalar in the range \(0,1\) or an integer scalar in the range \[1, n\).> ...
%! cvpartition ([1, 1, 1, 2, 2], "holdout", 6)
%!error <cvpartition: p value for 'holdout' must be a scalar in the range \(0,1\) or an integer scalar in the range \[1, n\).> ...
%! cvpartition ([1, 1, 1, 2, 2], "holdout", 6, "stratify", false)
%!error <cvpartition: k value for 'kfold' must be a numeric scalar.> ...
%! cvpartition ([1, 1, 1, 2, 2], "kfold", 'a')
%!error <cvpartition: k value for 'kfold' must be a numeric scalar.> ...
%! cvpartition ([1, 1, 1, 2, 2], "kfold", 'a', "stratify", true)
%!error <cvpartition: k value for 'kfold' must be a numeric scalar.> ...
%! cvpartition ([1, 1, 1, 2, 2], "kfold", [2, 3])
%!error <cvpartition: k value for 'kfold' must be a numeric scalar.> ...
%! cvpartition ([1, 1, 1, 2, 2], "kfold", [2, 3], "stratify", false)
%!error <cvpartition: k value for 'kfold' must be an integer scalar in the range \[1, n\).> ...\
%! cvpartition ([1, 1, 1, 2, 2], "kfold", 0)
%!error <cvpartition: k value for 'kfold' must be an integer scalar in the range \[1, n\).> ...\
%! cvpartition ([1, 1, 1, 2, 2], "kfold", 0, "stratify", true)
%!error <cvpartition: k value for 'kfold' must be an integer scalar in the range \[1, n\).> ...\
%! cvpartition ([1, 1, 1, 2, 2], "kfold", 1.5)
%!error <cvpartition: k value for 'kfold' must be an integer scalar in the range \[1, n\).> ...\
%! cvpartition ([1, 1, 1, 2, 2], "kfold", 1.5, "stratify", true)
%!error <cvpartition: k value for 'kfold' must be an integer scalar in the range \[1, n\).> ...\
%! cvpartition ([1, 1, 1, 2, 2], "kfold", 5)
%!error <cvpartition: k value for 'kfold' must be an integer scalar in the range \[1, n\).> ...\
%! cvpartition ([1, 1, 1, 2, 2], "kfold", 5, "stratify", true)
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
