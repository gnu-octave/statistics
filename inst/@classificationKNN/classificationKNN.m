## Copyright (C) 2023 Mohammed Azmat Khan <azmat.dev0@gmail.com>
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

classdef classificationKNN
## -*- texinfo -*-
## @deftypefn  {statistics} {@var{P} =} classificationKNN
## @deftypefnx {statistics} {@var{P} =} classificationKNN (@var{X}, @var{Y})
## @deftypefnx {statistics} {@var{P} =} classificationKNN (@var{X}, @var{Y}, @var{name}, @var{value})
##
## Create a kNN algorithm, classificationKNN Object using
## @var{X}, @var{Y} and other additional Name-Value pairs.
## Object of class classificationKNN can be used to store the training data
## and values can be altered via Object to predict labels and other data
## using built-in functions for new observations. It is recommended to use
## fitcknn to create an object of classificationKNN
##
## @code{@var{P} = classificationKNN (@var{X}, @var{Y})} returns an empty Object
## of class classificationKNN with empty properties.
##
## @code{@var{P} = classificationKNN (@var{X}, @var{Y})} returns an Object
## of class classificationKNN, with @var{X} as predictor data and @var{Y}
## with the class labels of observations in @var{X}.
##
## @code{@var{P} = classificationKNN (@var{X}, @var{Y}, @var{name}, @var{value})}
## returns an Object of class classificationKNN, with @var{X} as predictor data
## and @var{Y} with the class labels of observations in @var{X} with additional
## properties specified in @qcode{Name-Value} pairs.
##
## An object of @qcode{classificationKNN} contains the following properties :
##
## @multitable @columnfractions 0.05 0.2 0.75
## @headitem @tab @var{property} @tab @var{Description}
##
## @item @tab @qcode{"K"} @tab is the number of nearest neighbors to be found
## in the kNN search.  It must be a positive integer value and by default it is
## 1. @var{"K"} can be changed by using @code{@var{obj.K} = 10}.
##
## @item @tab @qcode{"weights"} @tab is a @math{Nx1} numeric non-negative matrix of
## the observational weights, each row in @var{weights} corresponds to the row
## in @var{Y} and indicates the relative importance or weight to be considered
## in calculating the Nearest-neighbour, negative values are removed before
## calculations if weights are specified. this propery cannot be changed via
## object. default value @qcode{@var{weight} = ones(rows(Y),1)}.
##
## @item @tab @qcode{"cost"} @tab is a @math{NxR} numeric matrix containing
## misclassification cost for the corresponding instances in @var{X} where
## @var{R} is the number of unique categories in @var{Y}. If an instance is
## correctly classified into its category the cost is calculated to be 1, If
## not then 0. cost matrix can be altered use @code{@var{obj.cost} = somecost}.
## default value @qcode{@var{cost} = ones(rows(X),numel(unique(Y)))}.
##
## @item @tab @qcode{"exponenet"} @tab is the Minkowski distance exponent and
## it must be a positive scalar.  This argument is only valid when the selected
## distance metric is @qcode{"minkowski"}.  By default it is 2.
##
## @item @tab @qcode{"scale"} @tab is the scale parameter for the standardized
## Euclidean distance and it must be a nonnegative numeric vector of equal
## length to the number of columns in @var{X}.  This argument is only valid when
## the selected distance metric is @qcode{"seuclidean"}, in which case each
## coordinate of @var{X} is scaled by the corresponding element of
## @qcode{"scale"}, as is each query point in @var{Y}.  By default, the scale
## parameter is the standard deviation of each coordinate in @var{X}.
##
## @item @tab @qcode{"Prior"} @tab is the matrix of prior probablities for each
## unique class in @var{Y}. @var{Prior} can be changed via object using
## @code{@var{obj.Prior} = somePrior}
##
## @item @tab @qcode{"cov"} @tab is the covariance matrix for computing the
## mahalanobis distance and it must be a positive definite matrix matching the
## the number of columns in @var{X}.  This argument is only valid when the
## selected distance metric is @qcode{"mahalanobis"}.
##
## @item @tab @qcode{"BucketSize"} @tab is the maximum number of data points in
## the leaf node of the Kd-tree and it must be a positive integer.  This
## argument is only valid when the selected search method is @qcode{"kdtree"}.
##
## @item @tab @qcode{"Distance"} @tab is the distance metric used by
## @code{knnsearch} as specified below:
## @end multitable
##
## @multitable @columnfractions 0.1 0.25 0.65
## @item @tab @qcode{"euclidean"} @tab Euclidean distance.
## @item @tab @qcode{"seuclidean"} @tab standardized Euclidean distance.  Each
## coordinate difference between the rows in @var{X} and the query matrix
## @var{Y} is scaled by dividing by the corresponding element of the standard
## deviation computed from @var{X}.  To specify a different scaling, use the
## @qcode{"scale"} name-value argument.
## @item @tab @qcode{"cityblock"} @tab City block distance.
## @item @tab @qcode{"chebychev"} @tab Chebychev distance (maximum coordinate
## difference).
## @item @tab @qcode{"minkowski"} @tab Minkowski distance.  The default exponent
## is 2.  To specify a different exponent, use the @qcode{"exponent"} name-value
## argument.
## @item @tab @qcode{"mahalanobis"} @tab Mahalanobis distance, computed using a
## positive definite covariance matrix.  To change the value of the covariance
## matrix, use the @qcode{"cov"} name-value argument.
## @item @tab @qcode{"cosine"} @tab Cosine distance.
## @item @tab @qcode{"correlation"} @tab One minus the sample linear correlation
## between observations (treated as sequences of values).
## @item @tab @qcode{"spearman"} @tab One minus the sample Spearman's rank
## correlation between observations (treated as sequences of values).
## @item @tab @qcode{"hamming"} @tab Hamming distance, which is the percentage
## of coordinates that differ.
## @item @tab @qcode{"jaccard"} @tab One minus the Jaccard coefficient, which is
## the percentage of nonzero coordinates that differ.
## @end multitable
##
## @multitable @columnfractions 0.05 0.2 0.75
## @item @tab @qcode{"NSMethod"} @tab is the nearest neighbor search method used
## by @code{knnsearch} as specified below.
## @end multitable
##
## @multitable @columnfractions 0.1 0.25 0.65
## @item @tab @qcode{"kdtree"} @tab Creates and uses a Kd-tree to find nearest
## neighbors.  @qcode{"kdtree"} is the default value when the number of columns
## in @var{X} is less than or equal to 10, @var{X} is not sparse, and the
## distance metric is @qcode{"euclidean"}, @qcode{"cityblock"},
## @qcode{"chebychev"}, or @qcode{"minkowski"}.  Otherwise, the default value is
## @qcode{"exhaustive"}.  This argument is only valid when the distance metric
## is one of the four aforementioned metrics.
## @item @tab @qcode{"exhaustive"} @tab Uses the exhaustive search algorithm by
## computing the distance values from all the points in @var{X} to each point in
## @var{Y}.
## @end multitable
##
## @multitable @columnfractions 0.05 0.2 0.75
## @item @tab @qcode{"IncludeTies"} @tab is a boolean flag to indicate if the
## returned values should contain the indices that have same distance as the
## @math{K^th} neighbor.  When @qcode{false}, @code{knnsearch} chooses the
## observation with the smallest index among the observations that have the same
## distance from a query point.  When @qcode{true}, @code{knnsearch} includes
## all nearest neighbors whose distances are equal to the @math{K^th} smallest
## distance in the output arguments. To specify @math{K}, use the @qcode{"K"}
## name-value pair argument.
##
## @item @tab @qcode{"breakties"} @tab is the method to break ties in kNN
## algorithm while deciding the label for an observation.
##
## @item @tab @qcode{"classnames"} @tab is a vector of unique classes in
## training data @var{Y}.If not given as name-value pair, the function
## calculates and stores the unique values in @var{Y}. @var{"classnames"} cannot
## be changed via @var{obj}.
##
## @item @tab @qcode{"NosClasses"} @tab is a numeric value that stores the
## number of unique classes in @var{Y}.If not given as name-value pair, the
## function calculates and stores the number of unique values in @var{Y}.
## @var{"NosClasses"} cannot be changed via @var{obj}.
##
## @item @tab @qcode{"Xclass"} @tab must be a @math{MxP} numeric matrix of
## query/new points that are to be classified into the clas labels.
## @var{Xclass} must have same numbers of columns as @var{X}. @var{Xclass}
## can be changed via object @var{obj} using @code{@var{obj.Xclass} = newXclass}.
## to predict labels for new set of values while keeping the rest of the model
## same with all its properties.
##
## @item @tab @qcode{"NN"} @tab stores the Nearest-Neighbours indices for the
## query points in @var{Xclass}. this property cannot be modified.
##
## @item @tab @qcode{"Score"} @tab store the predicted labels for the query
## points in @var{Xclass}.this property cannot be modified.
##
## @end multitable
##
## Methods of class classificationKNN are :
##
##
## @end deftypefn

  properties (Access = public)

    ## properties that can be set via object
    X = [];
    Y = [];
    Xclass = [];

    prior;
    cov;
    weights;
    exponent;
    distance;
    NSmethod;
    standerdize;

    ## additional arguments
   endproperties
   properties (Acesss = protected)
    k;
    NN;
    cost;                   #cost ,alterable
    bucketsize;
    Includeties;
    Breakties;
    NumObsv;
    Scale;                   # read only

    ## variable dependent properties
    classNames;
    NosClasses;
    Score;
    label;

  endproperties

  methods (Access = public)
    function this = classificationKNN (X, Y, varargin) ##obj constructor
      if (nargin < 2 && nargin != 0) ## return empty object when no arguments
        error ("@classificationKNN: too few arguments.");
      endif

      if (nargin >= 2) ## add upper condtion for nargin

        if (rows(X) != rows(Y))
          error ("@classificationKNN: number of rows in X and Y must be equal.");
        endif

        if ( !isnumeric(X) || !isfinite(X))
          error ("@classificationKNN: Invalid values in X.");
        else
          this.X = X;
          NumObsv = size (X,1);
        endif

        ## assign Y
        this.Y = Y;

        ## adding default values
            k = 1;
            cost  = [];                   #cost,alterable
            Scale = [];                   # read only
            prior = [];
            cov   = [];
            weights  = [];
            exponent = 2;
            bucketsize  = 50;
            distance    = "euclidean";
            NSmethod    = "exhaustive";
            standerdize = false;
            Includeties = false;
            Breakties   = "equal";
            NosClasses  = [];
            classNames  = [];
            Xclass      = [];

        while (numel (varargin) > 0)
          switch (tolower (varargin{1}))
          case "k"
            k = varargin{2};
          case "weights"
            weights = varargin{2};
          case "cost"
            cost = varargin{2};
          case "scale"
            Scale = varargin{2};
          case "prior"
            weights = varargin{2};
          case "cov"
            cov = varargin{2};
          case "exponent"
            exponent = varargin{2};
          case "bucketsize"
            bucketsize = varargin{2};
          case "distance"
            distance = varargin{2};
          case "includeties"
            Includeties = varargin{2};
          case "nsmethod"
            NSmethod = varargin{2};
          case "standerdize"
            standerdize = varargin{2};
          case "includeties"
            IncludeTies = varargin{2};
          case "breakties"
            breakties = varargin{2};
          case "classnames"
            classNames = varargin{2};
          case "nosclasses"
            NosClasses = varargin{2};
          case "xclass"
            Xclass = varargin{2};
          otherwise
            error ("@classificationKNN: Invalid NAME in optional pairs of arguments.");
          endswitch
          varargin(1:2) = [];
        endwhile


        ## assign properties related to Y

        if (isempty (classNames))
          this.classNames = unique (Y);
        endif
        if (isempty (NosClasses))
          this.NosClasses = numel (this.classNames);
        endif

        ##------checking optional parameters------##

        ## only checking if the arguments are provided and storing if they
        ## pass the check,
        ## If not provided the default values will be stored.

        ## check if Xclass
        if ( !isempty (Xclass))
          if (columns(X) != columns(Xclass))
            error ("@classificationKNN: number of columns in Xclass must be equal to X.");
          else
            this.Xclass = Xclass;
          endif
        endif

        ## check k
        if ( !isempty (k))
          if ( !isscalar(k) || !isnumeric(k) || k < 1 || k != round (k))
            error ("@classificationKNN: Invalid value of k.");
          else
            this.k = k;
          endif
        endif

        ## check weights
        if ( !isempty(weights))
          if ( !isvector(weights) || rows(X) != length(weights) || sum(weights < 0) != 0)
            error("@classificationKNN: Weights has invalid observarions.");
          else
            this.weights = weights;
          endif
        endif
        ## check minkowski distance exponent
        if ( !isempty (exponent))
          if ( !isscalar(exponent) || !isnumeric(exponent) || exponent < 0)
            error ("@classificationKNN: Invalid value of Minkowski Exponent.");
          else
            this.exponent = exponent;
          endif
        endif

        ## check cost
        if ( !isempty(cost))
          if ( !ismatrix (cost) || !isnumeric(cost) || columns (cost)!= numel(unique(Y)) || rows(cost) != rows(X) )
            error ("@classificationKNN: Invalid value in cost or the size of cost.");
          else
            this.cost = cost;
          endif
        endif
        ## check scale
        if ( !isempty(Scale))
          if ( !ismatrix (Scale) || any(Scale < 0) || numel(Scale) != rows(X))
            error ("@classificationKNN: Invalid value in Scale or the size of scale.");
          else
            this.Scale = Scale;
          endif
        endif
        ## check cov
        if ( !isempty(cov))
          if ( !strcmp(distance,"mahalanobis") || !ismatrix(cov) || !isnumeric(cov))
            error ("@classificationKNN: Invalid value in cov, cov can only be given for mahalanobis distance.");
          else
            this.cov = cov;
          endif
        endif
        ## check bucketsize
        if ( !isempty (bucketsize))
          if ( !isscalar(bucketsize) || bucketsize < 0)
            error ("@classificationKNN: Invalid value of bucketsize.");
          else
            this.bucketsize = bucketsize;
          endif
        endif

        ## standerdize
        if ( !isempty (standerdize))
          if ( !islogical(standerdize) && standerdize != 1 && standerdize != 0)
            error ("@classificationKNN: value of standerdize must be boolean.");
          else
            this.standerdize = standerdize;
          endif
        endif
        ## includeties
        if ( !isempty (Includeties))
          if ( !islogical(Includeties) && Includeties != 1 && Includeties != 0)
            error ("@classificationKNN: value of Includeties must be boolean");
          else
            this.Includeties = Includeties;
          endif
        endif

        if (!isempty (distance))
          this.distance = distance;
        endif

        if (!isempty (NSmethod))
          this.NSmethod = NSmethod;
        endif


        ## ------ param check and assign end ------ ##

      endif

    endfunction ## constructor


    ## predict function to predict the labels for Xclass
    function [label, Score, Cost] = predict (this)

      ## check Xclass
      if ( isempty (this.Xclass))
        error ("@classificationKNN.predict: Xclass is empty");
      elseif ( columns(this.X) != columns(this.Xclass))
        error ("@classificationKNN.predict: Xclass must have columns equal to X");
      endif

      ## check cost
      if (isempty (this.cost))
        ## if empty assign all cost = 1
        this.cost = ones (rows (this.X), this.NosClasses);
      endif

      if (isempty (this.X))
        ## no data in X
        label     = repmat(classNames(1,:),0,1);
        posterior = NaN(0,classNos);
        cost      = NaN(0,classNos);
      else
        ## calculate the NNs using knnsearch
        [idx, dist] = knnsearch (this.X, this.Xclass, "k", this.k, ...
                     "NSMethod", this.NSmethod, "Distance", this.distance, ...
                     "P", this.exponent, "Scale", this.Scale, "cov", this.cov, ...
                     "bucketsize", this.bucketsize, "sortindices", true, ...
                     "includeties",this.Includeties);

        [label, score, cost_temp] = predictlabel (this, idx);

        cost  = this.cost(rows(cost_temp),columns(cost_temp)) .* cost_temp;

        ## store predcited in the object variables

        this.NN    = idx;
        this.label = label;
        this.Score = score;
        this.cost  = cost;

      endif

    endfunction

    ## function predict label
    function [labels, score, cost_temp] = predictlabel (this, idx)
      ## assign intial values
      freq = [];
      wsum  = sum (this.weights);
      score = [];
      labels = [];
      cost_temp = [];

        for i = 1:rows (idx)
          if (!isempty (this.weights))
            ## weighted kNN
            for id = 1:numel (this.classNames)
              freq = [freq; (sum (strcmpi (this.classNames(id,1), this.Y(idx(i,:))) .* this.weights)) / wsum;]
              score_temp = (freq ./ this.k)';
            endfor
          else
            ## Non-weighted kNN
            for id = 1:size (this.classNames,1) #u{iu(:),2}
              freq(id,1) = (sum (strcmpi (this.classNames(id,1), this.Y(idx(i,:)))));
            endfor
            ## score calculation
            score_temp = (freq ./ this.k)';
            cost_temp  = [cost_temp; ones(1,this.NosClasses) - score_temp];
            score = [score; score_temp];
          endif

          [val, iu] = max (freq);

          ## set label for the index idx
          labels = [labels; this.classNames(iu,1)];

        endfor

    endfunction

  endmethods


endclassdef


%!demo
%! ##create a classificationKNN object using fisheriris dataset
%! load fisheriris
%! x = meas;
%! y = species;
%! a = classificationKNN (x, y, 'k', 5)

## Test constructor
%!test
%! a = classificationKNN ();
%! assert (class (a), "classificationKNN");
%! assert ({a.Breakties, a.Includeties, a.NN, a.NSmethod}, {[],[],[],[]})
%! assert ({a.NosClasses, a.NumObsv, a.Scale, a.Score, a.cost}, {[],[],[],[],[]})
%! assert ({a.X, a.Xclass, a.Y, a.bucketsize, a.classNames}, {[],[],[],[],[]})
%! assert ({a.cov, a.distance, a.exponent, a.k, a.label}, {[],[],[],[],[]})
%! assert ({a.prior, a.standerdize, a.weights}, {[],[],[]})

%!test
%! warning("off");
%! x = [1,2,3;4,5,6;7,8,9;3,2,1];
%! y = ['a';'a';'b';'b'];
%! a = classificationKNN (x, y);
%! assert (class (a), "classificationKNN");
%! assert ({a.X, a.Y, a.k},{[1,2,3;4,5,6;7,8,9;3,2,1], ['a';'a';'b';'b'], 1})
%! assert ({a.NSmethod, a.distance, a.bucketsize},{"exhaustive","euclidean",50})

%!test
%! warning("off");
%! x = [1,2,3;4,5,6;7,8,9;3,2,1];
%! y = ['a';'a';'b';'b'];
%! k = 10;
%! a = classificationKNN (x, y, "K" ,k);
%! assert (class (a), "classificationKNN");
%! assert ({a.X, a.Y, a.k},{[1,2,3;4,5,6;7,8,9;3,2,1], ['a';'a';'b';'b'], 10})
%! assert ({a.NSmethod, a.distance, a.bucketsize},{"exhaustive","euclidean",50})

%!test
%! warning("off");
%! x = [1,2,3;4,5,6;7,8,9;3,2,1];
%! y = ['a';'a';'b';'b'];
%! weights = ones (4,1);
%! a = classificationKNN (x, y, "weights" , weights);
%! assert (class (a), "classificationKNN");
%! assert ({a.X, a.Y, a.k},{[1,2,3;4,5,6;7,8,9;3,2,1], ['a';'a';'b';'b'], 1})
%! assert (a.weights, ones (4,1))
%! assert ({a.NSmethod, a.distance, a.bucketsize},{"exhaustive","euclidean",50})

%!test
%! warning("off");
%! x = [1,2,3;4,5,6;7,8,9;3,2,1];
%! y = ['a';'a';'b';'b'];
%! a = classificationKNN (x, y, "exponent" , 10);
%! assert (class (a), "classificationKNN");
%! assert ({a.X, a.Y, a.k},{[1,2,3;4,5,6;7,8,9;3,2,1], ['a';'a';'b';'b'], 1})
%! assert (a.exponent, 10)
%! assert ({a.NSmethod, a.distance, a.bucketsize},{"exhaustive","euclidean",50})

%!test
%! warning("off");
%! x = [1,2,3;4,5,6;7,8,9;3,2,1];
%! y = ['a';'a';'b';'b'];
%! cov = rand (4,1);
%! a = classificationKNN (x, y, "cov" , cov, 'distance', 'mahalanobis');
%! assert (class (a), "classificationKNN");
%! assert ({a.X, a.Y, a.k},{[1,2,3;4,5,6;7,8,9;3,2,1], ['a';'a';'b';'b'], 1})
%! assert (a.cov, cov)
%! assert ({a.NSmethod, a.distance, a.bucketsize},{"exhaustive","mahalanobis",50})

%!test
%! warning("off");
%! x = [1,2,3;4,5,6;7,8,9;3,2,1];
%! y = ['a';'a';'b';'b'];
%! a = classificationKNN (x, y, "bucketsize" , 20, 'distance', 'mahalanobis');
%! assert (class (a), "classificationKNN");
%! assert ({a.X, a.Y, a.k},{[1,2,3;4,5,6;7,8,9;3,2,1], ['a';'a';'b';'b'], 1})
%! assert ({a.NSmethod,a.distance,a.bucketsize},{"exhaustive","mahalanobis",20})

%!test
%! warning("off");
%! x = [1,2,3;4,5,6;7,8,9;3,2,1];
%! y = ['a';'a';'b';'b'];
%! a = classificationKNN (x, y, 'standerdize', true);
%! assert (class (a), "classificationKNN");
%! assert ({a.X, a.Y, a.k},{[1,2,3;4,5,6;7,8,9;3,2,1], ['a';'a';'b';'b'], 1})
%! assert (a.standerdize, true);
%! assert ({a.NSmethod,a.distance,a.bucketsize},{"exhaustive","euclidean",50})

%!test
%! warning("off");
%! x = [1,2,3;4,5,6;7,8,9;3,2,1];
%! y = ['a';'a';'b';'b'];
%! a = classificationKNN (x, y, 'includeties', true);
%! assert (class (a), "classificationKNN");
%! assert ({a.X, a.Y, a.k},{[1,2,3;4,5,6;7,8,9;3,2,1], ['a';'a';'b';'b'], 1})
%! assert (a.Includeties, true);
%! assert ({a.NSmethod,a.distance,a.bucketsize},{"exhaustive","euclidean",50})

%!test
%! warning("off");
%! x = [1,2,3;4,5,6;7,8,9;3,2,1];
%! y = ['a';'a';'b';'b'];
%! cost = ones (4,2);
%! a = classificationKNN (x, y, 'cost', cost );
%! assert (class (a), "classificationKNN")
%! assert ({a.X, a.Y, a.k},{[1,2,3;4,5,6;7,8,9;3,2,1], ['a';'a';'b';'b'], 1})
%! assert (a.cost, [1,1;1,1;1,1;1,1])
%! assert ({a.NSmethod,a.distance,a.bucketsize},{"exhaustive","euclidean",50})

%!test
%! warning("off");
%! x = [1,2,3;4,5,6;7,8,9;3,2,1];
%! y = ['a';'a';'b';'b'];
%! scale = [1,2,3,4];
%! a = classificationKNN (x, y, 'scale', scale );
%! assert (class (a), "classificationKNN")
%! assert ({a.X, a.Y, a.k},{[1,2,3;4,5,6;7,8,9;3,2,1], ['a';'a';'b';'b'], 1})
%! assert (a.Scale, [1,2,3,4])
%! assert ({a.NSmethod,a.distance,a.bucketsize},{"exhaustive","euclidean",50})


## Test input validation
%!error<@classificationKNN: too few arguments.>classificationKNN(ones(4,1))
%!error<@classificationKNN: number of rows in X and Y must be equal.> ...
%! classificationKNN(ones (4,2), ones (1,4))
%!error<@classificationKNN: Invalid values in X.> ...
%! classificationKNN([1;2;3;'a';4], ones (5,1))
%!error<@classificationKNN: Invalid NAME in optional pairs of arguments.> ...
%! classificationKNN(ones (4,2), ones (4,1), "some",'some')
%!error<@classificationKNN: number of columns in Xclass must be equal to X.> ...
%! classificationKNN(ones (4,2),ones (4,1), "Xclass", ones (1,4))
%!error<@classificationKNN: Invalid value of k.> ...
%! classificationKNN(ones (4,2),ones (4,1), "K", NaN)
%!error<@classificationKNN: Invalid value of k.> ...
%! classificationKNN(ones (4,2),ones (4,1), "K", -5)
%!error<@classificationKNN: Invalid value of k.> ...
%! classificationKNN(ones (4,2),ones (4,1), "K", 3.14)
%!error<@classificationKNN: Weights has invalid observarions.> ...
%! classificationKNN(ones(4,2),ones(4,1), "weights", ones(2,2))
%!error<@classificationKNN: Weights has invalid observarions.> ...
%! classificationKNN(ones(4,2),ones(4,1), "weights", [1;2;3;'a';4])
%!error<@classificationKNN: Invalid value of Minkowski Exponent.> ...
%! classificationKNN(ones(4,2),ones(4,1), "exponent", -5)
%!error<@classificationKNN: Invalid value in cost or the size of cost.> ...
%! classificationKNN(ones(4,2),ones(4,1), "cost", ones (2,2))
%!error<@classificationKNN: Invalid value in cost or the size of cost.> ...
%! classificationKNN(ones(4,2),ones(4,1), "cost", [1;2;3;'a';4])
%!error<@classificationKNN: Invalid value in Scale or the size of scale.> ...
%! classificationKNN(ones(4,2),ones(4,1), "Scale", [1;2;3;'a';4])
%!error<@classificationKNN: Invalid value in Scale or the size of scale.> ...
%! classificationKNN(ones(4,2),ones(4,1), "Scale", zeros (3,3))
%!error<@classificationKNN: Invalid value in cov, cov can only be given for mahalanobis distance.> ...
%! classificationKNN(ones(4,2),ones(4,1), "cov", ones (4,1))
%!error<@classificationKNN: Invalid value of bucketsize.> ...
%! classificationKNN(ones(4,2),ones(4,1), "bucketsize", -5)
%!error<@classificationKNN: value of standerdize must be boolean.> ...
%! classificationKNN(ones(4,2),ones(4,1), "standerdize", 5)
%!error<@classificationKNN: number of columns in Xclass must be equal to X.> ...
%! classificationKNN (ones (4,3), ones (4,1), 'Xclass', ones(4,1))
