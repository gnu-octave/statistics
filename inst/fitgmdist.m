## Copyright (C) 2018 John Donoghue <john.donoghue@ieee.org>
## Copyright (C) 2015 Lachlan Andrew <lachlanbis@gmail.com>
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

## -*- texinfo -*-
## @deftypefn {Function File} {@var{GMdist} =} fitgmdist (@var{data}, @var{k}, @var{param1}, @var{value1}, @dots{})
## Fit a Gaussian mixture model with @var{k} components to @var{data}.
## Each row of @var{data} is a data sample.  Each column is a variable.
##
## Optional parameters are:
##    @itemize
##    @item 'start':  initialization conditions.  Possible values are:
##       @itemize
##       @item 'randSample' (default) takes means uniformly from rows of data
##       @item 'plus'       use k-means++ to initialize means
##       @item 'cluster'    Performs an initial clustering with 10% of the data
##       @item vector       A vector whose length is the number of rows in data,
##                    and whose values are 1 to k specify the components
##                    each row is initially allocated to.  The mean, variance
##                    and weight of each component is calculated from that
##       @item structure    with elements  mu,  Sigma  ComponentProportion
##       @end itemize
##       For 'randSample', 'plus' and 'cluster', the initial variance of each
##       component is the variance of the entire data sample.
##
##    @item 'Replicates'    Number of random restarts to perform
##
##    @item 'RegularizationValue'
##    @item 'Regularize'  A small number added to the diagonal entries
##                    of the covariance to prevent singular covariances
##
##    @item 'SharedCovariance'
##    @item 'SharedCov' (logical) True if all components must share the
##                    same variance, to reduce the number of free parameters
##
##    @item 'CovarianceType'
##    @item 'CovType'       (string). Possible values are:
##       @itemize
##       @item 'full'       (default)  Allow arbitrary covariance matrices
##       @item 'diagonal'   Force covariances to be diagonal, to reduce the
##                    number of free parameters.
##       @end itemize
##
##    @item 'Option'        A structure with all of the following fields:
##       @itemize
##       @item 'MaxIter'    Maximum number of EM iterations (default 100)
##       @item 'TolFun'     Threshold increase in likelihood to terminate EM
##                          (default 1e-6)
##       @item 'Display'
##          @itemize
##          @item 'off' (default): display nothing
##          @item 'final': display the number of iterations and likelihood
##                         once execution completes
##          @item 'iter':  display the above after each iteration
##          @end itemize
##       @end itemize
##    @item 'Weight'  A column vector or n-by-2 matrix. The first column
##                    consists of non-negative weights given to the
##                    samples.
##                    If these are all integers, this is equivalent
##                    to specifying @var{weight}(i) copies of row i of
##                    @var{data}, but potentially faster.
##
##                    If a row of @var{data} is used to represent samples
##                    that are similar but not identical, then the second
##                    column of @var{weight} indicates the variance of
##                    those original samples.  Specifically, in the EM
##                    algorithm, the contribution of row i towards the
##                    variance is set to at least @var{weight}(i,2), to
##                    prevent spurious components with zero variance.
##    @end itemize
##    
## @seealso{gmdistribution, kmeans}
## @end deftypefn

function obj = fitgmdist(data, k, varargin)

  if nargin < 2 || mod (nargin, 2) == 1
    print_usage;
  endif

  [~, prop] = parseparams (varargin);

  ## defaults for options
  diagonalCovar  = false;      # "full".  (true is "diagonal")
  sharedCovar    = false;
  start          = "randSample";
  replicates     = 1;
  option.MaxIter= 100;
  option.TolFun = 1e-6;
  option.Display= "off";       # "off" (1 is "final", 2 is "iter")
  Regularizer    = 0;
  weights        = [];         # Each row i counts as "weights(i,1)" rows


  # Remove rows containing NaN / NA
  data = data(!any (isnan (data), 2),:);

  # used for getting the number of samples
  nRows = rows (data);
  nCols = columns (data);

  # Parse options
  while (!isempty (prop))
    try
      switch (lower (prop{1}))
        case {"sharedcovariance",...
              "sharedcov"},          sharedCovar    = prop{2};
        case {"covariancetype",...
             "covartype"},           diagonalCovar  = prop{2};
        case {"regularizationvalue",...
              "regularize"},         Regularizer    = prop{2};
        case "replicates",           replicates     = prop{2};
        case "start",                start          = prop{2};
        case "weights",              weights        = prop{2};

        case "option"
            option.MaxIter = prop{2}.MaxIter;
            option.TolFun  = prop{2}.TolFun;
            option.Display = prop{2}.Display;

        otherwise
          error ("fitgmdist: Unknown option %s", prop{1});
      endswitch
    catch ME
      if (length (prop) < 2)
        error ("fitgmdist: Option '%s' has no argument", prop{1});
      else
        rethrow (ME)
      endif
    end_try_catch
    prop = prop(3:end);
  endwhile

  # Process options


  # check for the "replicates" property
  try
    if isempty (1:replicates)
      error ("fitgmdist: replicates must be positive");
    endif
  catch
    error ("fitgmdist: invalid number of replicates");
  end_try_catch

  # check for the "option" property
  MaxIter = option.MaxIter;
  TolFun  = option.TolFun;
  switch (lower (option.Display))
    case "off", Display = 0;
    case "final", Display = 1;
    case "iter", Display = 2;
    case "notify", Display = 0;
    otherwise, error ("fitgmdist: Unknown Display option %s", option.Display);
  endswitch

  try
    p = ones(1, k) / k;         # Default is uniform component proportions
  catch ME
    if (!isscalar (k) || !isnumeric (k))
      error ("fitgmdist: The second argument must be a numeric scalar");
    else
      rethrow (ME)
    endif
  end_try_catch
  # check for the "start" property
  if (ischar (start))
    start = lower (start);
    switch (start)
      case {"randsample", "plus", "cluster", "randsamplep", "plusp", "clusterp"}
      otherwise
        error ("fitgmdist: Unknown Start value %s\n", start);
    endswitch
    component_order_free = true;
  else
    component_order_free = false;
    if (!ismatrix (start) || !isnumeric (start))
      try
        mu    = start.mu;
        Sigma = start.Sigma;
        if (isfield (start, 'ComponentProprition'))
          p = start.ComponentProportion(:)';
        end
        if (any (size (data, 2) ~= [size(mu,2), size(Sigma)]) || ...
            any (k ~= [size(mu,1), size(p,2)]))
          error ('fitgmdist: Start parameter has mismatched dimensions');
        endif
      catch
        error ("fitgmdist: invalid start parameter");
      end_try_catch
      if (isfield (start, "ComponentProprition"))
        p = start.ComponentProportion(:)';
      endif
    else
      validIndices = 0;
      mu = zeros (k, nRows);
      Sigma = zeros (nRows, nRows, k);
      for i = 1:k
        idx = (start == i);
        validIndices = validIndices + sum (idx);
        mu(i,:) = mean (data(idx,:));
        Sigma(:,:,i) = cov (data(idx,:)) + Regularizer*eye (nCols);
      endfor
      if (validIndices < nRows)
        error ("fitgmdist: Start is numeric, but is not integers between 1 and k");
      endif
    endif
    start = [];       # so that variance isn't recalculated later
    replicates = 1;   # Will be the same each time anyway
  endif

  # check for the "SharedCovariance" property
  if (!islogical (sharedCovar))
    error ("fitgmdist: SharedCoveriance must be logical true or false");
  endif

  # check for the "CovarianceType" property
  if (!islogical (diagonalCovar))
    try
      if (strcmpi (diagonalCovar, "diagonal"))
        diagonalCovar = true;
      elseif (strcmpi (diagonalCovar, "full"))
        diagonalCovar = false;
      else
        error ("fitgmdist: CovarianceType must be Full or Diagonal");
      endif
    catch
        error ("fitgmdist: CovarianceType must be 'Full' or 'Diagonal'");
    end_try_catch
  endif

  # check for the "Regularizer" property
  try
    if (Regularizer < 0)
      error ("fitgmdist: Regularizer must be non-negative");
    endif
  catch ME
    if (!isscalar (Regularizer) || !isnumeric (Regularizer))
      error ("fitgmdist: Regularizer must be a numeric scalar");
    else
      rethrow (ME)
    endif
  end_try_catch

  # check for the "Weights" property and the matrix
  try
    if (!isempty (weights))
      if (columns (weights) > 2 || any (weights(:) < 0))
        error ("fitgmdist: weights must be a nonnegative numeric dx1 or dx2 matrix");
      endif
      if (rows (weights) != nRows)
        error ("fitgmdist: number of weights %d must match number of samples %d",...
               rows (weights), nRows)
      endif
      non_zero = (weights(:,1) > 0);
      weights = weights(non_zero,:);
      data    = data   (non_zero,:);

      nRows = rows (data);
      raw_samples = sum (weights(:,1));
    else
      raw_samples = nRows;
    endif

  # Validate the matrix
    if (!isreal (data(k,1)))
      error ("fitgmdist: first input argument must be a DxN real data matrix");
    endif
  catch ME
    if (!isnumeric (data) || !ismatrix (data) || !isreal (data))
      error ("fitgmdist: first input argument must be a DxN real data matrix");
    elseif (k > nRows || k < 0)
      if (exists("non_zero", "var") && k <= length(non_zero))
        error ("fitgmdist: The number of non-zero weights (%d) must be at least the number of components (%d)", nRows, k);
      else
        error ("fitgmdist: The number of components (%d) must be a positive number less than the number of data rows (%d)", k, nRows);
      endif
    elseif (!ismatrix (weights) || !isnumeric (weights))
      error ("fitgmdist: weights must be a nonnegative numeric dx1 or dx2 matrix");
    else
      rethrow (ME)
    endif
  end_try_catch
  #if k == 1
  #  replicates = 1;
  #endif


  # Done processing options
  #######################################

  # used to hold the probability of each class, given each data vector
  try
    p_x_l = zeros (nRows, k);     # probability of observation x given class l

    best = -realmax;
    best_params = [];
    diag_slice = 1:(nCols+1):(nCols)^2;

    # Create index slices to calculate symmetric completion of upper triangular Mx
    lower_half = zeros(nCols*(nCols-1)/2,1);
    upper_half = zeros(nCols*(nCols-1)/2,1);
    i = 1;
    for rw = 1:nCols
      for cl = rw+1:nCols
        upper_half(i) = sub2ind([nCols, nCols], rw, cl);
        lower_half(i) = sub2ind([nCols, nCols], cl, rw);
        i = i + 1;
      endfor
    endfor

    for rep = 1:replicates
      if (!isempty (start))
        # Initialize the means
        switch (start)
          case {"randsample"}
            if (isempty (weights))
              idx = randperm (nRows, k);
            else
              idx = randsample (nRows, k, false, weights);
            endif
            mu = data(idx, :);
          case {"plus"}               # k-means++, by Arthur and Vassilios
            mu(1,:) = data(randi (nRows),:);
            d = inf (nRows, 1);       # Distance to nearest centroid so far
            for i = 2:k
              d = min (d, sum (bsxfun (@minus, data, mu(i-1, :)).^2, 2));
                      # pick next sample with prob. prop to dist.*weights
              if (isempty (weights))
                cs = cumsum (d);
              else
                cs = cumsum (d .* weights(:,1));
              endif
              mu(i,:) = data(find (cs > rand * cs(end), 1), :);
            endfor
          case {"cluster"}
            subsamp = max (k, ceil (nRows/10));
            if (isempty (weights))
              idx = randperm (nRows, subsamp);
            else
              idx = randsample (nRows, subsamp, false, weights);
            endif
            [~, mu] = kmeans (data(idx), k, "start", "sample");
        endswitch

        # Initialize the variance, unless set explicitly
        #
          Sigma = var (data) + Regularizer;
          if (!diagonalCovar)
            Sigma = diag (Sigma);
          endif
          if (!sharedCovar)
            Sigma = repmat (Sigma, [1, 1, k]);
          endif
      endif

      # Run the algorithm
      iter = 1;

      log_likeli = -inf;
      incr = 1;

      while (incr > TolFun && iter <= MaxIter)
        iter = iter + 1;
        #######################################
        # "E step"
        # Calculate probability of class  l  given observations
        for i = 1:k
          if (sharedCovar)
            sig = Sigma;
          else
            sig = Sigma(:,:,i);
          endif
          if (diagonalCovar)
            sig = diag(sig);
          endif
          try
            p_x_l (:, i) = mvnpdf (data, mu(i, :), sig);
          catch ME
            if (strfind (ME.message, "positive definite"))
              error ("fitgmdist: Covariance is not positive definite.  Increase RegularizationValue");
            else
              rethrow (ME)
            endif
          end_try_catch
        endfor
            # Bayes' rule
        p_x_l = bsxfun (@times, p_x_l, p);                # weight by priors
        p_l_x = bsxfun (@rdivide, p_x_l, sum (p_x_l, 2)); # Normalize

        #######################################
        # "M step"
        # Calculate new parameters
        if (!isempty (weights))
          p_l_x = bsxfun (@times, p_l_x, weights(:,1));
        endif

        sum_p_l_x = sum(p_l_x);   # row vec of \sum_{data} p(class|data,params)

        p  = sum_p_l_x / raw_samples;                       # new proportions
        mu = bsxfun(@rdivide, p_l_x' * data, sum_p_l_x');   # new means
        if (sharedCovar)
          sumSigma = zeros (size (Sigma(:,:,1))); # diagonalCovar gives size
        endif
        for i = 1:k
          # Sigma
          deviation = bsxfun(@minus, data, mu(i,:));
          lhs = bsxfun(@times, p_l_x(:,i), deviation);

          # Calculate covariance
          # Iterate either over elements of the covariance matrix, since
          # there should be fewer of those than rows of data.
          for rw = 1:nCols
            for cl = rw:nCols
              sig(rw,cl) = lhs(:,rw)' * deviation(:,cl);
            endfor
          endfor
          sig(lower_half) = sig(upper_half);

          sig = sig/sum_p_l_x(i) + Regularizer*eye (nCols);

          if (columns (weights) > 1)   # don't give "singleton" clusters low var
            sig(diag_slice) = max (sig(diag_slice), weights(i,2));
          endif

          if (diagonalCovar)
            sig = diag(sig)';
          endif

          if (sharedCovar)
            sumSigma = sumSigma + sig * p(i);         # Heuristic. Should it use
          else                                        # old p?  Something else?
            Sigma(:,:,i) = sig;
          endif
        endfor
        if (sharedCovar)
          Sigma = sumSigma;
        endif
        #######################################

        # calculate the new (and relative change in) log-likelihood
        if (isempty (weights))
          new_log_likeli = sum (log (sum (p_x_l, 2)));
        else
          new_log_likeli = sum (weights(:,1) .* log (sum (p_x_l, 2)));
        endif
        incr  = (new_log_likeli - log_likeli)/max(1,abs(new_log_likeli));
        if Display == 2
          fprintf("iter %d  log-likelihood %g\n", iter-1, new_log_likeli);
          %disp(mu);
        endif
        log_likeli = new_log_likeli;
      endwhile
      if (log_likeli > best)
        best = log_likeli;
        best_params.mu    = mu;
        best_params.Sigma = Sigma;
        best_params.p     = p;
      endif
    endfor
  catch ME
    try
      if (1 < MaxIter), end
    catch
      error ("fitgmdist: invalid MaxIter");
    end_try_catch
    rethrow (ME)
  end_try_catch

  # List components in descending order of proportion,
  # unless the order was implicitly specified by "start"
  if (component_order_free)
    [~, idx] = sort (-best_params.p);
    best_params.p  = best_params.p (idx);
    best_params.mu = best_params.mu(idx,:);
    if (!sharedCovar)
      best_params.Sigma = best_params.Sigma(:,:,idx);
    endif
  endif

  # Calculate number of parameters
  if (diagonalCovar)
    params = nCols;
  else
    params = nCols * (nCols+1) / 2;
  endif
  params = params*size (Sigma, 3) + 2*rows (mu) - 1;

  extra.NegativeLogLikelihood = -best;
  extra.AIC = -2*(best - params);
  extra.BIC = -2*best + params * log (raw_samples);
  extra.Converged = (incr <= TolFun);
  extra.NumIterations = iter-1;
  extra.RegularizationValue = Regularizer;

  # This works in Octave, but not in Matlab
  #obj = gmdistribution (best_params.mu, best_params.Sigma, best_params.p', extra);
  obj = gmdistribution (best_params.mu, best_params.Sigma, best_params.p');
  if (Display == 1)
    fprintf ("  %d iterations   log-likelihood = %g\n", ...
              -extra.NegativeLogLikelihood, extra.NumIterations);
  endif
endfunction

%!xdemo <50286>
%! ## Generate a two-cluster problem
%! C1 = randn (100, 2) + 1;
%! C2 = randn (100, 2) - 1;
%! data = [C1; C2];
%!
%! ## Perform clustering
%! GMModel = fitgmdist (data, 2);
%!
%! ## Plot the result
%! figure
%! [heights, bins] = hist3([C1; C2]);
%! [xx, yy] = meshgrid(bins{1}, bins{2});
%! bbins = [xx(:), yy(:)];
%! contour (reshape (GMModel.pdf (bbins), heights));
%! hold on
%! plot (centers (:, 1), centers (:, 2), "kv", "markersize", 10);
%! hold off
