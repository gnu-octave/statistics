## Copyright (C) 2013-2019 Fernando Damian Nieuwveldt <fdnieuwveldt@gmail.com>
## Copyright (C) 2021 Stefano Guidoni <ilguido@users.sf.net>
##
## This file is part of the statistics package for GNU Octave.
##
## This program is free software; you can redistribute it and/or
## modify it under the terms of the GNU General Public License
## as published by the Free Software Foundation; either version 3
## of the License, or (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn  {statistics} @var{coeff} = pca(@var{X})
## @deftypefnx {statistics} @var{coeff} = pca(@var{X}, Name, Value)
## @deftypefnx {statistics} [@var{coeff}, @var{score}, @var{latent}] = pca(@dots{})
## @deftypefnx {statistics} [@var{coeff}, @var{score}, @var{latent}, @var{tsquared}] = pca(@dots{})
## @deftypefnx {statistics} [@var{coeff}, @var{score}, @var{latent}, @var{tsquared}, @var{explained}, @var{mu}] = pca(@dots{})
##
## Performs a principal component analysis on a data matrix X.
##
## A principal component analysis of a data matrix of @code{n} observations in a
## @code{p}-dimensional space returns a @code{p}-by-@code{p} transformation
## matrix, to perform a change of basis on the data. The first component of the
## new basis is the direction that maximizes the variance of the projected data.
##
## Input argument:
## @itemize @bullet
## @item
## @var{x} : a @code{n}-by-@code{p} data matrix
## @end itemize
##
## Pair arguments:
## @itemize @bullet
## @item
## @code{Algorithm} : the algorithm to use, it can be either @code{eig},
## for eigenvalue decomposition, or @code{svd} (default), for singular value
## decomposition
## @item
## @code{Centered} : boolean indicator for centering the observation data, it is
## @code{true} by default
## @item
## @code{Economy} : boolean indicator for the economy size output, it is
## @code{true} by default; @code{pca} returns only the elements of @var{latent}
## that are not necessarily zero, and the corresponding columns of @var{coeff}
## and @var{score}, that is, when @code{n <= p}, only the first @code{n - 1}
## @item
## @code{NumComponents} : the number of components @code{k} to return, if
## @code{k < p}, then only the first @code{k} columns of @var{coeff}
## and @var{score} are returned
## @item
## @code{Rows} : action to take with missing values, it can be either
## @code{complete} (default), missing values are removed before computation,
## @code{pairwise} (only with algorithm @code{eig}), the covariance of rows with
## missing data is computed using the available data, but the covariance matrix
## could be not positive definite, which triggers the termination of @code{pca},
## @code{complete}, missing values are not allowed, @code{pca} terminates with
## an error if there are any
## @item
## @code{Weights} : observation weights, it is a vector of positive values of
## length @code{n}
## @item
## @code{VariableWeights} : variable weights, it can be either a vector of
## positive values of length @code{p} or the string @code{variance} to use the
## sample variance as weights
## @end itemize
##
## Return values:
## @itemize @bullet
## @item
## @var{coeff} : the principal component coefficients, a @code{p}-by-@code{p}
## transformation matrix
## @item
## @var{score} : the principal component scores, the representation of @var{x}
## in the principal component space
## @item
## @var{latent} : the principal component variances, i.e., the eigenvalues of
## the covariance matrix of @var{x}
## @item
## @var{tsquared} : Hotelling's T-squared Statistic for each observation in
## @var{x}
## @item
## @var{explained} : the percentage of the variance explained by each principal
## component
## @item
## @var{mu} : the estimated mean of each variable of @var{x}, it is zero if the
## data are not centered
## @end itemize
##
## Matlab compatibility note:  the alternating least square method 'als' and
## associated options 'Coeff0', 'Score0', and 'Options' are not yet implemented
##
## @subheading References
##
## @enumerate
## @item
## Jolliffe, I. T., Principal Component Analysis, 2nd Edition, Springer, 2002
##
## @end enumerate
## @end deftypefn

## BUGS:
## - tsquared with weights (xtest below)

##FIXME
## -- can change isnan to ismissing once the latter is implemented in Octave
## -- change mystd to std and remove the helper function mystd once weighting is available in the Octave function

function [coeff, score, latent, tsquared, explained, mu] = pca (X, varargin)

  if (nargin < 1)
     print_usage ();
  endif

  [nobs, nvars] = size (X);

  ## default options
  optAlgorithmS = "svd";
  optCenteredB = true;
  optEconomyB = true;
  optNumComponentsI = nvars;
  optWeights = [];
  optVariableWeights = [];
  optRowsB = false;
  TF = [];

  ## parse parameters
  pair_index = 1;
  while (pair_index <= (nargin - 1))
    switch (lower (varargin{pair_index}))
      ## decomposition algorithm: singular value decomposition, eigenvalue
      ## decomposition or (currently unavailable) alternating least square
      case "algorithm"
        optAlgorithmS = varargin{pair_index + 1};
        switch (optAlgorithmS)
          case {"svd", "eig"}
            ;
          case "als"
            error ("pca: alternating least square algorithm not implemented");
          otherwise
            error ("pca: invalid algorithm %s", optAlgorithmS);
        endswitch
      ## centering of the columns, around the mean
      case "centered"
        if (isbool (varargin{pair_index + 1}))
          optCenteredB = varargin{pair_index + 1};
        else
          error ("pca: 'centered' requires a boolean value");
        endif
      ## limit the size of the output to the degrees of freedom, when a smaller
      ## number than the number of variables
      case "economy"
        if (isbool (varargin{pair_index + 1}))
          optEconomyB = varargin{pair_index + 1};
        else
          error ("pca: 'economy' requires a boolean value");
        endif
      ## choose the number of components to show
      case "numcomponents"
        optNumComponentsI = varargin{pair_index + 1};
        if ((! isscalar (optNumComponentsI)) ||
            (! isnumeric (optNumComponentsI)) ||
            optNumComponentsI != floor (optNumComponentsI) ||
            optNumComponentsI <= 0 ||
            optNumComponentsI > nvars)
          error (["pca: the number of components must be a positive integer"...
                  "number smaller or equal to the number of variables"]);
        endif
      ## observation weights: some observations can be more accurate than others
      case "weights"
        optWeights = varargin{pair_index + 1};
        if ((! isvector (optWeights)) ||
            length (optWeights) != nobs ||
            length (find (optWeights < 0)) > 0)
          error ("pca: weights must be a numerical array of positive numbers");
        endif

        if (rows (optWeights) == 1 )
          optWeights = transpose (optWeights);
        endif
      ## variable weights: weights used for the variables
      case "variableweights"
        optVariableWeights = varargin{pair_index + 1};
        if (ischar (optVariableWeights) &&
            strcmpi (optVariableWeights, "variance"))
          optVariableWeights = "variance"; # take care of this later
        elseif ((! isvector (optVariableWeights)) ||
            length (optVariableWeights) != nvars ||
            (! isnumeric (optVariableWeights)) ||
            length (find (optVariableWeights < 0)) > 0)
          error (["pca: variable weights must be a numerical array of "...
                  "positive numbers or the string 'variance'"]);
        else
          optVariableWeights = 1 ./ sqrt (optVariableWeights);

          ## it is used as a row vector
          if (columns (optVariableWeights) == 1 )
            optVariableWeights = transpose (optVariableWeights);
          endif
        endif
      ## rows: policy for missing values
      case "rows"
        switch (varargin{pair_index + 1})
          case "complete"
            optRowsB = false;
          case "pairwise"
            optRowsB = true;
          case "all"
            if (any (isnan (X)))
              error (["pca: when all rows are requested the dataset cannot"...
                      " include NaN values"]);
            endif
          otherwise
            error ("pca: %s is an invalid value for rows", ...
                   varargin{pair_index + 1});
        endswitch
      case {"coeff0", "score0", "options"}
        error ("pca: parameter %s is only valid with the 'als' method, which is not yet implemented", varargin{pair_index});
      otherwise
        error ("pca: unknown property %s", varargin{pair_index});
    endswitch

    pair_index += 2;
  endwhile

  ## Preparing the dataset according to the chosen policy for missing values
  if (optRowsB)
    if (! strcmp (optAlgorithmS, "eig"))
      optAlgorithmS = "eig";
      warning (["pca: setting algorithm to 'eig' because 'rows' option is "...
                "set to 'pairwise'"]);
    endif

    TF = isnan (X);
    missingRows = zeros (nobs, 1);
    nmissing = 0;
  else
    ## "complete": remove all the rows with missing values
    TF = isnan (X);
    missingRows = any (TF, 2);
    nmissing = sum (missingRows);
  endif

  ## indices of the available rows
  ridcs = find (missingRows == 0);

  ## Center the columns to mean zero if requested
  if (optCenteredB)
    if (isempty (optWeights) && nmissing == 0 && ! optRowsB)
      ## no weights and no missing values
      mu = mean (X);
    elseif (nmissing == 0 && ! optRowsB)
      ## weighted observations: some observations are more valuable, i.e. they
      ## can be trusted more
      mu = sum (optWeights .* X) ./ sum (optWeights);
    else
      ## missing values: the mean is computed column by column
      mu = zeros (1, nvars);

      if (isempty (optWeights))
        for iter = 1 : nvars
          mu(iter) = mean (X(find (TF(:, iter) == 0), iter));
        endfor
      else
        ## weighted mean with missing data
        for iter = 1 : nvars
          mu(iter) =  sum (X(find (TF(:, iter) == 0), iter) .* ...
                          optWeights(find (TF(:, iter) == 0))) ./ ...
                      sum (optWeights(find (TF(:, iter) == 0)));
        endfor
      endif
    endif

    Xc = X - mu;
  else
    Xc = X;

    ## The mean of the variables of the original dataset:
    ## return zero if the dataset is not centered
    mu = zeros (1, nvars);
  endif

  ## Change the columns according to the variable weights
  if (! isempty (optVariableWeights))
    if (ischar (optVariableWeights))
      if (isempty (optWeights))
        sqrtBias = 1; # see below
        optVariableWeights = std (X);
      else
        ## unbiased variance estimation: the bias when using reliability weights
        ## is 1 - var(weights) / std(weigths)^2
        sqrtBias = sqrt (1 - (sumsq (optWeights) / sum (optWeights) ^ 2));
        optVariableWeights = mystd (X, optWeights) / sqrtBias;
      endif
    endif
    Xc = Xc ./ optVariableWeights;
  endif

  ## Compute the observation weight matrix
  if (isempty (optWeights))
    Wd = eye (nobs - nmissing);
  else
    Wd = diag (optWeights) ./ sum (optWeights);
  endif

  ## Compute the coefficients
  switch (optAlgorithmS)
    case "svd"
      ## Check if there are more variables than observations
      if (nvars <= nobs)
        [U, S, coeff] = svd (sqrt (Wd) * Xc(ridcs,:), "econ");
      else
        ## Calculate the svd on the transpose matrix, much faster
        if (optEconomyB)
          [coeff, S, V] = svd (Xc(ridcs,:)' * sqrt (Wd), "econ");
        else
          [coeff, S, V] = svd (Xc(ridcs,:)' * sqrt (Wd));
        endif
      endif
    case "eig"
      ## this method requires the computation of the sample covariance matrix
      if (optRowsB)
        ## pairwise:
        ## in this case the degrees of freedom for each element of the matrix
        ## are equal to the number of valid rows for the couple of columns
        ## used to compute the element
        Xpairwise = Xc;
        Xpairwise(find (isnan (Xc))) = 0;

        Ndegrees = (nobs - 1) * ones (nvars, nvars);
        for i_iter = 1 : nvars
          for j_iter = i_iter : nvars
            Ndegrees(i_iter, j_iter) = Ndegrees(i_iter, j_iter) - ...
                                       sum (any (TF(:,[i_iter j_iter]), 2));
            Ndegrees(j_iter, i_iter) = Ndegrees(i_iter, j_iter);
          endfor
        endfor

        Mcov = Xpairwise' * Wd * Xpairwise ./ Ndegrees;
      else
        ## the degrees of freedom are not really important here
        ndegrees = nobs - nmissing - 1;
        Mcov = Xc(ridcs, :)' * Wd * Xc(ridcs, :) / ndegrees;
      endif

      [coeff, S] = eigs (Mcov, nvars);
  endswitch

  ## Change the coefficients according to the variable weights
  if (! isempty (optVariableWeights))
    coeff = coeff .* transpose (optVariableWeights);
  endif

  ## MATLAB compatibility: the sign convention is that the
  ## greatest absolute value for each column is positive
  switchSignV = find (max (coeff) < abs (min (coeff)));
  if (! isempty (switchSignV))
    coeff(:, switchSignV) = -1 * coeff(:, switchSignV);
  endif

  ## Compute the scores
  if (nargout > 1)
    ## This is for the score when using variable weights, it is not really
    ## a new definition of Xc
    if (! isempty (optVariableWeights))
      Xc = Xc ./ optVariableWeights;
    endif

    ## Get the Scores
    score = Xc(ridcs,:) * coeff;

    ## Get the rank of the score matrix
    r = rank (score);

    ## If there is missing data, put it back
    ## FIXME: this needs tests
    if (nmissing)
      scoretmp = zeros (nobs, nvars);
      scoretmp(find (missingRows == 0), :) = score;
      scoretmp(find (missingRows), :) = NaN;
      score = scoretmp;
    endif

    ## Only use the first r columns, pad rest with zeros if economy != true
    score = score(:, 1:r) ;

    if (! optEconomyB)
      score = [score, (zeros (nobs , nvars-r))];
    else
      coeff = coeff(: , 1:r);
    endif
  endif

  ## Compute the variances
  if (nargout > 2)
    ## degrees of freedom: n - 1 for centered data
    if (optCenteredB)
      dof = size (Xc(ridcs,:), 1) - 1;
    else
      dof = size (Xc(ridcs,:), 1);
    endif

    ## This is the same as the eigenvalues of the covariance matrix of X
    if (strcmp (optAlgorithmS, "eig"))
      latent = diag (S, 0);
    else
      latent  = (diag (S'*S) / dof)(1:r);
    endif

    ## If observation weights were used, we need to scale back these values
    if (! isempty (optWeights))
      latent = latent .* sum (optWeights(ridcs));
    endif

    if (! optEconomyB)
      latent= [latent; (zeros (nvars - r, 1))];
    endif
  endif

  ## Compute the Hotelling T-square statistics
  ## MATLAB compatibility: when using weighted observations the T-square
  ## statistics differ by some rounding error
  if (nargout > 3)
    ## Calculate the Hotelling T-Square statistic for the observations
    ## formally: tsquared = sumsq (zscore (score(:, 1:r)),2);
    if (! isempty (optWeights))
      ## probably splitting the weights, using the square roots, is not the
      ## best solution, numerically
      weightedScore = score .* sqrt (optWeights);
      tsquared = mahal (weightedScore(ridcs, 1:r), weightedScore(ridcs, 1:r))...
                 ./ optWeights;
    else
      tsquared = mahal (score(ridcs, 1:r), score(ridcs, 1:r));
    endif
  endif

  ## Compute the variance explained by each principal component
  if (nargout > 4)
    explained = 100 * latent / sum (latent);
  endif

  ## When a number of components is chosen, the coefficients and score matrix
  ## only show that number of columns
  if (optNumComponentsI != nvars)
    coeff = coeff(:, 1:optNumComponentsI);
    score = score(:, 1:optNumComponentsI);
  endif
endfunction

#return the weighted standard deviation
function retval = mystd (x, w)
  (dim = find (size(x) != 1, 1)) || (dim = 1);
  den = sum (w);
  mu = sum (w .* x, dim) ./ sum (w);
  retval = sum (w .* ((x - mu) .^ 2), dim) / den;
  retval = sqrt (retval);
endfunction

%!shared COEFF,SCORE,latent,tsquare,m,x,R,V,lambda,i,S,F

#NIST Engineering Statistics Handbook example (6.5.5.2)
%!test
%! x=[7 4 3
%!    4 1 8
%!    6 3 5
%!    8 6 1
%!    8 5 7
%!    7 2 9
%!    5 3 3
%!    9 5 8
%!    7 4 5
%!    8 2 2];
%! R = corrcoef (x);
%! [V, lambda] = eig (R);
%! [~, i] = sort(diag(lambda), "descend"); #arrange largest PC first
%! S = V(:, i) * diag(sqrt(diag(lambda)(i)));
%!assert(diag(S(:, 1:2)*S(:, 1:2)'), [0.8662; 0.8420; 0.9876], 1E-4); #contribution of first 2 PCs to each original variable
%! B = V(:, i) * diag( 1./ sqrt(diag(lambda)(i)));
%! F = zscore(x)*B;
%! [COEFF,SCORE,latent,tsquare] = pca(zscore(x, 1));
%!assert(tsquare,sumsq(F, 2),1E4*eps);

%!test
%! x=[1,2,3;2,1,3]';
%! [COEFF,SCORE,latent,tsquare] = pca(x, "Economy", false);
%! m=[sqrt(2),sqrt(2);sqrt(2),-sqrt(2);-2*sqrt(2),0]/2;
%! m(:,1) = m(:,1)*sign(COEFF(1,1));
%! m(:,2) = m(:,2)*sign(COEFF(1,2));

%!assert(COEFF,m(1:2,:),10*eps);
%!assert(SCORE,-m,10*eps);
%!assert(latent,[1.5;.5],10*eps);
%!assert(tsquare,[4;4;4]/3,10*eps);

#test with observation weights (using Matlab's results for this case as a reference)
%! [COEFF,SCORE,latent,tsquare] = pca(x, "Economy", false, "weights", [1 2 1], "variableweights", "variance");
%!assert(COEFF, [0.632455532033676 -0.632455532033676; 0.741619848709566 0.741619848709566], 10*eps);
%!assert(SCORE, [-0.622019449426284 0.959119380657905; -0.505649896847432 -0.505649896847431; 1.633319243121148 0.052180413036957], 10*eps);
%!assert(latent, [1.783001790889027; 0.716998209110974], 10*eps);
%!xtest assert(tsquare, [1.5; 0.5; 1.5], 10*eps);  #currently, [4; 2; 4]/3 is actually returned; see comments above

%!test
%! x=x';
%! [COEFF,SCORE,latent,tsquare] = pca(x, "Economy", false);
%! m=[sqrt(2),sqrt(2),0;-sqrt(2),sqrt(2),0;0,0,2]/2;
%! m(:,1) = m(:,1)*sign(COEFF(1,1));
%! m(:,2) = m(:,2)*sign(COEFF(1,2));
%! m(:,3) = m(:,3)*sign(COEFF(3,3));

%!assert(COEFF,m,10*eps);
%!assert(SCORE(:,1),-m(1:2,1),10*eps);
%!assert(SCORE(:,2:3),zeros(2),10*eps);
%!assert(latent,[1;0;0],10*eps);
%!assert(tsquare,[0.5;0.5],10*eps)

%!test
%! [COEFF,SCORE,latent,tsquare] = pca(x);

%!assert(COEFF,m(:, 1),10*eps);
%!assert(SCORE,-m(1:2,1),10*eps);
%!assert(latent,[1],10*eps);
%!assert(tsquare,[0.5;0.5],10*eps)

%!error <invalid algorithm> pca([1 2; 3 4], "Algorithm", "xxx")
%!error <'centered' requires a boolean value> pca([1 2; 3 4], "Centered", "xxx")
%!error <must be a positive integer> pca([1 2; 3 4], "NumComponents", -4)
%!error <invalid value for rows> pca([1 2; 3 4], "Rows", 1)
%!error <weights must be> pca([1 2; 3 4], "Weights", [1 2 3])
%!error <weights must be> pca([1 2; 3 4], "Weights", [-1 2])
%!error <variable weights must be> pca([1 2; 3 4], "VariableWeights", [-1 2])
%!error <variable weights must be> pca([1 2; 3 4], "VariableWeights", "xxx")
%!error <unknown property> pca([1 2; 3 4], "XXX", 1)
