## Copyright (C) 2023 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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

## -*- texinfo -*-
## @deftypefn  {statistics} {@var{TF} =} isoutlier (@var{x})
## @deftypefnx {statistics} {@var{TF} =} isoutlier (@var{x}, @var{method})
## @deftypefnx {statistics} {@var{TF} =} isoutlier (@var{x}, @qcode{"percentiles"}, @var{threshold})
## @deftypefnx {statistics} {@var{TF} =} isoutlier (@var{x}, @var{movmethod}, @var{window})
## @deftypefnx {statistics} {@var{TF} =} isoutlier (@dots{}, @var{dim})
## @deftypefnx {statistics} {@var{TF} =} isoutlier (@dots{}, @var{Name}, @var{Value})
## @deftypefnx {statistics} {[@var{TF}, @var{L}, @var{U}, @var{C}] =} isoutlier (@dots{})
##
## Find outliers in data
##
## @code{isoutlier (@var{x})} returns a logical array whose elements are true
## when an outlier is detected in the corresponding element of @var{x}.
## @code{isoutlier} treats NaNs as missing values and removes them.
##
## @itemize
## @item
## If @var{x} is a matrix, then @code{isoutlier} operates on each column of
## @var{x} separately.
## @item
## If @var{x} is a multidimensional array, then @code{isoutlier} operates along
## the first dimension of @var{x} whose size does not equal 1.
## @end itemize
##
## By default, an outlier is a value that is more than three scaled median
## absolute deviations (MAD) from the median.  The scaled median is defined as
## @code{c*median(abs(A-median(A)))}, where @code{c=-1/(sqrt(2)*erfcinv(3/2))}.
##
## @code{isoutlier (@var{x}, @var{method})} specifies a method for detecting
## outliers.  The following methods are available:
##
## @multitable @columnfractions 0.13 0.02 0.8
## @headitem Method @tab @tab Description
## @item @qcode{"median"} @tab @tab Outliers are defined as elements more than
## three scaled MAD from the median.
## @item @qcode{"mean"} @tab @tab Outliers are defined as elements more than
## three standard deviations from the mean.
## @item @qcode{"quartiles"} @tab @tab Outliers are defined as elements more
## than 1.5 interquartile ranges above the upper quartile (75 percent) or below
## the lower quartile (25 percent).  This method is useful when the data in
## @var{x} is not normally distributed.
## @item @qcode{"grubbs"} @tab @tab Outliers are detected using Grubbs’ test for
## outliers, which removes one outlier per iteration based on hypothesis
## testing.  This method assumes that the data in @var{x} is normally
## distributed.
## @item @qcode{"gesd"} @tab @tab Outliers are detected using the generalized
## extreme Studentized deviate test for outliers.  This iterative method is
## similar to @qcode{"grubbs"}, but can perform better when there are multiple
## outliers masking each other.
## @end multitable
##
## @code{isoutlier (@var{x}, @qcode{"percentiles"}, @var{threshold})} detects
## outliers based on a percentile thresholds,  specified as a two-element row
## vector whose elements are in the interval @math{[0, 100]}.  The first element
## indicates the lower percentile threshold, and the second element indicates
## the upper percentile threshold. The first element of threshold must be less
## than the second element.
##
## @code{isoutlier (@var{x}, @var{movmethod}, @var{window})} specifies a moving
## method for detecting outliers.  The following methods are available:
##
## @multitable @columnfractions 0.13 0.02 0.8
## @headitem Method @tab @tab Description
## @item @qcode{"movmedian"} @tab @tab Outliers are defined as elements more
## than three local scaled MAD from the local median over a window length
## specified by @var{window}.
## @item @qcode{"movmean"} @tab @tab Outliers are defined as elements more than
## three local standard deviations from the from the local mean over a window
## length specified by @var{window}.
## @end multitable
##
## @var{window} must be a positive integer scalar or a two-element vector of
## positive integers.  When @var{window} is a scalar,  if it is an odd number,
## the window is centered about the current element and contains
## @qcode{@var{window} - 1} neighboring elements.  If even, then the window is
## centered about the current and previous elements.  When @var{window} is a
## two-element vector of positive integers @math{[nb, na]}, the window contains
## the current element, @math{nb} elements before the current element, and
## @math{na} elements after the current element.  When @qcode{"SamplePoints"}
## are also specified, @var{window} can take any real positive values (either as
## a scalar or a two-element vector) and in this case, the windows are computed
## relative to the sample points.
##
## @var{dim} specifies the operating dimension and it must be a positive integer
## scalar.  If not specified, then, by default, @code{isoutlier} operates along
## the first non-singleton dimension of @var{x}.
##
## The following optional parameters can be specified as @var{Name}/@var{Value}
## paired arguments.
##
## @itemize
## @item @qcode{"SamplePoints"} can be specified as a vector of sample points
## with equal length as the operating dimension.  The sample points represent
## the x-axis location of the data and must be sorted and contain unique
## elements.  Sample points do not need to be uniformly sampled.  By default,
## the vector is @qcode{[1, 2, 3, @dots{}, @var{n}]}, where
## @qcode{@var{n} = size (@var{x}, @var{dim})}.  You can use unequally spaced
## @qcode{"SamplePoints"} to define a variable-length window for one of the
## moving methods available.
##
## @item @qcode{"ThresholdFactor"} can be specified as a nonnegative scalar.
## For methods @qcode{"median"} and @qcode{"movmedian"}, the detection threshold
## factor replaces the number of scaled MAD, which is 3 by default.  For methods
## @qcode{"mean"} and @qcode{"movmean"}, the detection threshold factor replaces
## the number of standard deviations, which is 3 by default.  For methods
## @qcode{"grubbs"} and @qcode{"gesd"}, the detection threshold factor ranges
## from 0 to 1, specifying the critical @math{alpha}-value of the respective
## test, and it is 0.05 by default.  For the @qcode{"quartiles"} method, the
## detection threshold factor replaces the number of interquartile ranges, which
## is 1.5 by default.  @qcode{"ThresholdFactor"} is not supported for the
## @qcode{"quartiles"} method.
##
## @item @qcode{"MaxNumOutliers"} is only relevant to the @qcode{"gesd"} method
## and it must be a positive integer scalar specifying the maximum number of
## outliers returned by the @qcode{"gesd"} method.  By default, it is the
## integer nearest to the 10% of the number of elements along the operating
## dimension in @var{x}.  The @qcode{"gesd"} method assumes the nonoutlier input
## data is sampled from an approximate normal distribution.  When the data is
## not sampled in this way, the number of returned outliers might exceed the
## @qcode{MaxNumOutliers} value.
## @end itemize
##
## @code{[@var{TF}, @var{L}, @var{U}, @var{C}] = isoutlier (@dots{})} returns
## up to 4 output arguments as described below.
##
## @itemize
## @item @var{TF} is the outlier indicator with the same size a @var{x}.
##
## @item @var{L} is the lower threshold used by the outlier detection method.
## If @var{method} is used for outlier detection, then @var{L} has the same size
## as @var{x} in all dimensions except for the operating dimension where the
## length is 1.  If @var{movmethod} is used, then @var{L} has the same size as
## @var{x}.
##
## @item @var{U} is the upper threshold used by the outlier detection method.
## If @var{method} is used for outlier detection, then @var{U} has the same size
## as @var{x} in all dimensions except for the operating dimension where the
## length is 1.  If @var{movmethod} is used, then @var{U} has the same size as
## @var{x}.
##
## @item @var{C} is the center value used by the outlier detection method.
## If @var{method} is used for outlier detection, then @var{C} has the same size
## as @var{x} in all dimensions except for the operating dimension where the
## length is 1.  If @var{movmethod} is used, then @var{C} has the same size as
## @var{x}.  For @qcode{"median"}, @qcode{"movmedian"}, @qcode{"mean"}, and
## @qcode{"movmean"} methods, @var{C} is computed by taking into acount the
## outlier values.  For @qcode{"grubbs"} and @qcode{"gesd"} methods, @var{C} is
## computed by excluding the outliers.  For the @qcode{"percentiles"} method,
## @var{C} is the average between @var{U} and @var{L} thresholds.
## @end itemize
##
## @seealso{filloutliers, rmoutliers, ismissing}
## @end deftypefn

function [TF, L, U, C] = isoutlier (x, varargin)

  ## Check for valid input data
  if (nargin < 1)
    print_usage;
  endif

  ## Handle case if X is a scalar
  if (isscalar (x))
    TF = false;
    L = x;
    U = x;
    C = x;
    return
  endif

  ## Add defaults
  dim = [];
  method = "median";
  window = [];
  SamplePoints = [];
  ThresholdFactor = 3;
  MaxNumOutliers = [];

  ## MATLAB's constant for scaled Median Absolute Deviation
  ## c = -1 / (sqrt (2) * erfcinv (3/2))
  c = 1.482602218505602;

  ## Parse exrta arguments
  while (numel (varargin) > 0)
    if (ischar (varargin{1}))
      switch (lower (varargin{1}))
        case "median"
          method = "median";
          ThresholdFactor = 3;
          varargin(1) = [];

        case "mean"
          method = "mean";
          ThresholdFactor = 3;
          varargin(1) = [];

        case "quartiles"
          method = "quartiles";
          ThresholdFactor = 1.5;
          varargin(1) = [];

        case "grubbs"
          method = "grubbs";
          ThresholdFactor = 0.05;
          varargin(1) = [];

        case "gesd"
          method = "gesd";
          ThresholdFactor = 0.05;
          MaxNumOutliers = [];
          varargin(1) = [];

        case "movmedian"
          method = "movmedian";
          window = varargin{2};
          if (! isnumeric (window) || numel (window) < 1  ||
                numel (window) > 2 || any (window <= 0))
            error (strcat ("isoutlier: WINDOW must be a positive scalar", ...
                           " or a two-element vector of positive values"));
          endif
          varargin([1:2]) = [];

        case "movmean"
          method = "movmean";
          window = varargin{2};
          if (! isnumeric (window) || numel (window) < 1  ||
                numel (window) > 2 || any (window <= 0))
            error (strcat ("isoutlier: WINDOW must be a positive scalar", ...
                           " or a two-element vector of positive values"));
          endif
          varargin([1:2]) = [];

        case "percentiles"
          method = "percentiles";
          threshold = varargin{2};
          if (! isnumeric (threshold) || ! (numel (threshold) == 2))
            error (strcat ("isoutlier: THRESHOLD must be a two-element", ...
                           " vector whose elements are in the interval", ...
                           " [0, 100]."));
          endif
          if (! (threshold(1) < threshold(2)) || threshold(1) < 0 ||
              threshold(2) > 100)
            error (strcat ("isoutlier: THRESHOLD must be a two-element", ...
                           " vector whose elements are in the interval", ...
                           " [0, 100]."));
          endif
          varargin([1:2]) = [];

        case "samplepoints"
          SamplePoints = varargin{2};
          if (! isvector (SamplePoints) || isscalar (SamplePoints))
            error ("isoutlier: sample points must be a vector.");
          endif
          if (numel (unique (SamplePoints)) != numel (SamplePoints))
            error ("isoutlier: sample points must be unique.");
          endif
          if (any (sort (SamplePoints) != SamplePoints))
            error ("isoutlier: sample points must be sorted.");
          endif
          varargin([1:2]) = [];

        case "thresholdfactor"
          ThresholdFactor = varargin{2};
          if (! isscalar (ThresholdFactor) || ThresholdFactor <= 0)
            error ("isoutlier: threshold factor must be a nonnegative scalar.");
          endif
          varargin([1:2]) = [];

        case "maxnumoutliers"
          MaxNumOutliers = varargin{2};
          if (! isscalar (MaxNumOutliers) || MaxNumOutliers <= 0 ||
              ! (fix (MaxNumOutliers) == MaxNumOutliers))
            error (strcat ("isoutlier: maximum outlier count must be a", ...
                           " positive integer scalar."));
          endif
          varargin([1:2]) = [];

        otherwise
          error ("isoutlier: invalid input argument.");
      endswitch

    elseif (isnumeric (varargin{1}))
      dim = varargin{1};
      if (! fix (dim) == dim || dim < 1 || ! isscalar (dim) ||
          ! isscalar (varargin{1}))
        error ("isoutlier: DIM must be a positive integer scalar.");
      endif
      varargin(1) = [];
    else
      error ("isoutlier: invalid input argument.");
    endif
  endwhile

  ## Find 1st operating dimension (if empty)
  if (isempty (dim))
    szx = size (x);
    (dim = find (szx != 1, 1)) || (dim = 1);
  endif

  ## Check for valid WINDOW unless Sample Points are given
  if (isempty (SamplePoints) && ! isempty (window))
    if (! all (fix (window) == window))
      error (strcat ("isoutlier: WINDOW must be a positive integer", ...
                     " scalar or a two-element vector of positive", ...
                     " integers, unless SamplePoints are defined."));
    endif
  endif

  ## Check for valid value of ThresholdFactor for 'grubbs' and 'geds' methods
  if (any (strcmpi (method, {"grubbs", "gesd"})) && ThresholdFactor > 1)
    error (strcat ("isoutlier: threshold factor must must be in [0 1]", ...
                   " range for 'grubbs' and 'gesd' methods."));
  endif

  ## Switch methods
  switch method
    case "median"
      [L, U, C] = median_method (x, dim, ThresholdFactor, c);
      TF = x < L | x > U;
    case "mean"
      [L, U, C] = mean_method (x, dim, ThresholdFactor);
      TF = x < L | x > U;
    case "quartiles"
      [L, U, C] = quartiles_method (x, dim, ThresholdFactor);
      TF = x < L | x > U;
    case "grubbs"
      [TF, L, U, C] = grubbs_method (x, dim, ThresholdFactor);
    case "gesd"
      [L, U, C] = gesd_method (x, dim, ThresholdFactor, MaxNumOutliers);
      TF = x < L | x > U;
    case "movmedian"
      sp = SamplePoints;
      [L, U, C] = movmedian_method (x, dim, ThresholdFactor, c, window, sp);
      TF = x < L | x > U;
    case "movmean"
      sp = SamplePoints;
      [L, U, C] = movmean_method (x, dim, ThresholdFactor, window, sp);
      TF = x < L | x > U;
    case "percentiles"
      [L, U, C] = percentiles_method (x, dim, threshold);
      TF = x < L | x > U;
  endswitch

endfunction

## Find lower and upper outlier thresholds with median method
function [L, U, C] = median_method (x, dim, ThresholdFactor, c)
  C = median (x, dim, "omitnan");
  sMAD = c * mad (x, 1, dim);
  L = C - ThresholdFactor * sMAD;
  U = C + ThresholdFactor * sMAD;
endfunction

## Find lower and upper outlier thresholds with mean method
function [L, U, M] = mean_method (x, dim, ThresholdFactor)
  M = mean (x, dim, "omitnan");
  S = std (x, [], dim, "omitnan");
  L = M - ThresholdFactor * S;
  U = M + ThresholdFactor * S;
endfunction

## Find lower and upper outlier thresholds with quartiles method
function [L, U, C] = quartiles_method (x, dim, ThresholdFactor)
  Q = quantile (x, dim);
  C = Q(3);
  L = Q(2) - (Q(4) - Q(2)) * ThresholdFactor;
  U = Q(4) + (Q(4) - Q(2)) * ThresholdFactor;
endfunction

## Find lower and upper outlier thresholds with grubbs method
function [TF, L, U, C] = grubbs_method (x, dim, ThresholdFactor)
  ## Move the desired dim to be the 1st dimension (rows)
  szx   = size (x);                       # size of dimensions
  N     = szx(dim);                       # elements in operating dimension
  nd    = length (szx);                   # number of dimensions
  dperm = [dim, 1:(dim-1), (dim+1):nd];   # permutation of dimensions
  x     = permute (x, dperm);             # permute dims to first dimension
  ncols = prod (szx(dperm(2:end)));       # rest of dimensions as single column
  x     = reshape (x, N, ncols);          # reshape input
  ## Create return matrices
  L = zeros ([1, szx(dperm(2:end))]);
  U = L;
  C = L;
  TF = false (size (x));
  ## Apply processing to each column
  for i = 1:ncols
    tmp_x = x(:,i);
    TFvec = [(i-1)*size(x,1)+1:i*size(x,1)];
    TFvec(isnan (tmp_x)) = [];
    tmp_x(isnan (tmp_x)) = [];
    ## Search for outliers (one at a time)
    while (true)
      ## Get descriptive statistics
      n = length (tmp_x);
      C(i) = mean (tmp_x);
      S = std (tmp_x);
      ## Locate maximum deviation from mean
      dif_x = abs (tmp_x - C(i));
      max_x = max (dif_x);
      loc_x = find (dif_x == max_x, 1);
      ## Calculate Grubbs's critical value
      t_crit = tinv (ThresholdFactor / (2 * n), n - 2);
      G_crit = ((n - 1) / sqrt (n)) * abs (t_crit) / sqrt (n - 2 + t_crit ^ 2);
      ## Check hypothesis
      if (max_x / S > G_crit)
        tmp_x(loc_x) = [];
        TF(TFvec(loc_x)) = true;
        TFvec(loc_x) = [];
      else
        break;
      endif
    endwhile
    L(i) = C(i) - S * G_crit;
    U(i) = C(i) + S * G_crit;
  endfor
  ## Restore shape
  TF = ipermute (TF, dperm);
  L  = ipermute (L, dperm);
  U  = ipermute (U, dperm);
  C  = ipermute (C, dperm);
endfunction

## Find lower and upper outlier thresholds with gesd method
function [L, U, C] = gesd_method (x, dim, ThresholdFactor, MaxNumOutliers)
  ## Add default value in MaxNumOutliers (if empty)
  szx   = size (x);
  N     = szx(dim);
  if (isempty (MaxNumOutliers))
    MaxNumOutliers = ceil (N * 0.1);
  endif
  ## Move the desired dim to be the 1st dimension (rows)
  nd    = length (szx);                   # number of dimensions
  dperm = [dim, 1:(dim-1), (dim+1):nd];   # permutation of dimensions
  x     = permute (x, dperm);             # permute dims to first dimension
  ncols = prod (szx(dperm(2:end)));       # rest of dimensions as single column
  x     = reshape (x, N, ncols);          # reshape input
  ## Create return matrices
  L = zeros ([1, szx(dperm(2:end))]);
  U = L;
  C = L;
  TF = false (size (x));
  ## Apply processing to each column
  for i = 1:ncols
    tmp_x = x(:,i);
    vec_x = [(i-1)*size(x,1)+1:i*size(x,1)];
    vec_x(isnan (tmp_x)) = [];
    tmp_x(isnan (tmp_x)) = [];
    n = length (tmp_x);
    if (n > 1)
      mean_x = zeros (MaxNumOutliers,1);
      S = zeros (MaxNumOutliers,1);
      lambda = zeros (MaxNumOutliers,1);
      R =  zeros (MaxNumOutliers,1);
      Ridx = zeros (MaxNumOutliers,1);
      ## Search for given outliers
      for j = 1:MaxNumOutliers
        ## Get descriptive statistics
        mean_x(j) = mean (tmp_x);
        S(j) = std (tmp_x);
        ## Locate maximum deviation from mean
        dif_x = abs (tmp_x - mean_x(j));
        max_x = max (dif_x);
        loc_x = find (dif_x == max_x, 1);
        ## Calculate R
        R(j) = max_x / S(j);
        tmp_x(loc_x) = [];
        Ridx(j) = vec_x(loc_x);
        vec_x(loc_x) = [];
        ## Calculate lambda
        pp = 1 - ThresholdFactor / (2 * (n - j + 1));
        t = tinv (pp, n - j - 1);
        lambda(j) = (n - j) * t / sqrt ((n - j - 1 + t .^ 2) * (n - j + 1));
      endfor
      ## Find largest index
      idx = find (R > lambda, 1, "last");
      if (isempty (idx))
          TFidx = 1;
      else
          TFidx = min (idx + 1, MaxNumOutliers);
      endif
      L(i) = mean_x(TFidx) - S(TFidx) * lambda(TFidx);
      U(i) = mean_x(TFidx) + S(TFidx) * lambda(TFidx);
      C(i) = mean_x(TFidx);
    endif
  endfor
  ## Restore shape
  L  = ipermute (L, dperm);
  U  = ipermute (U, dperm);
  C  = ipermute (C, dperm);
endfunction

## Find lower and upper outlier thresholds with movmedian method
function [L, U, C] = movmedian_method (x, dim, ThresholdFactor, c, window, sp);
  szx = size (x);
  N   = szx(dim);
  ## Constrain window to the element in the operating dimension
  if (numel (window) == 1 && window > N)
    window = N;
  elseif (numel (window) == 2 && sum (window) > N)
    window = N;
  endif
  if (isempty (sp))
    FCN = @(x) median (x, "omitnan");
    C = movfun (FCN, x, window, "dim", dim);
    FCN = @(x) mad (x, 1);
    MAD = movfun (FCN, x, window, "dim", dim);
  else
    ## Check that sample points(sp) have the N elements
    if (numel (sp) != N)
      error (strcat ("isoutlier: sample points must have the same size", ...
                     " as the operating dimension."));
    endif
    ## Move the desired dim to be the 1st dimension (rows)
    nd    = length (szx);                 # number of dimensions
    dperm = [dim, 1:(dim-1), (dim+1):nd]; # permutation of dimensions
    x     = permute (x, dperm);           # permute dims to first dimension
    ncols = prod (szx(dperm(2:end)));     # rest of dimensions as single column
    x     = reshape (x, N, ncols);        # reshape input
    ## Find beg+end from window
    if (numel (window) == 2)
      w_lo = window(1);
      w_hi = window(2);
    else
      if (mod (window, 2) == 1)
        w_lo = w_hi = (window - 1) / 2;
      else
        w_lo = window / 2;
        w_hi = w_lo - 1;
      endif
    endif
    ## Create return matrices
    C = zeros (size (x));
    MAD = C;
    for i = 1:ncols
      tmp_x = x(:,i);
      for j = 1:N
        cp = sp - sp(j);
        nb = length (cp(cp < 0 & cp >= -w_lo));
        na = length (cp(cp > 0 & cp <= w_hi));
        sp_ind = [j-nb:j+na];
        C(j,i) = median (tmp_x(sp_ind), "omitnan");
        MAD(j,i) = mad (tmp_x(sp_ind), 1);
      endfor
    endfor
    ## Restore shape
    C = ipermute (C, dperm);
    MAD = ipermute (MAD, dperm);
  endif
  ## Compute scaled MAD
  sMAD = c * MAD;
  L = C - ThresholdFactor * sMAD;
  U = C + ThresholdFactor * sMAD;
endfunction

## Find lower and upper outlier thresholds with movmean method
function [L, U, M] = movmean_method (x, dim, ThresholdFactor, window, sp);
  ## Constrain window to the element in the operating dimension
  szx = size (x);
  N = szx(dim);
  if (numel (window) == 1 && window > N)
    window = N;
  elseif (numel (window) == 2 && sum (window) > N)
    window = N;
  endif
  if (isempty (sp))
    FCN = @(x) mean (x, "omitnan");
    M = movfun (FCN, x, window, "dim", dim);
    FCN = @(x) std (x, [], "omitnan");
    S = movfun (FCN, x, window, "dim", dim);
  else
    ## Check that sample points(sp) have the N elements
    if (numel (sp) != N)
      error (strcat ("isoutlier: sample points must have the same size", ...
                     " as the operating dimension."));
    endif
    ## Move the desired dim to be the 1st dimension (rows)
    nd    = length (szx);                 # number of dimensions
    dperm = [dim, 1:(dim-1), (dim+1):nd]; # permutation of dimensions
    x     = permute (x, dperm);           # permute dims to first dimension
    ncols = prod (szx(dperm(2:end)));     # rest of dimensions as single column
    x     = reshape (x, N, ncols);        # reshape input
    ## Find beg+end from window
    if (numel (window) == 2)
      w_lo = window(1);
      w_hi = window(2);
    else
      if (mod (window, 2) == 1)
        w_lo = w_hi = (window - 1) / 2;
      else
        w_lo = window / 2;
        w_hi = w_lo - 1;
      endif
    endif
    ## Create return matrices
    M = zeros (size (x));
    S = M;
    for i = 1:ncols
      tmp_x = x(:,i);
      for j = 1:N
        cp = sp - sp(j);
        nb = length (cp(cp < 0 & cp >= -w_lo));
        na = length (cp(cp > 0 & cp <= w_hi));
        sp_ind = [j-nb:j+na];
        M(j,i) = mean (tmp_x(sp_ind), "omitnan");
        S(j,i) = std (tmp_x(sp_ind), [], "omitnan");
      endfor
    endfor
    ## Restore shape
    M = ipermute (M, dperm);
    S = ipermute (S, dperm);
  endif
  L = M - ThresholdFactor * S;
  U = M + ThresholdFactor * S;
endfunction

## Find lower and upper outlier thresholds with percentiles method
function [L, U, C] = percentiles_method (x, dim, threshold)
  P = [threshold(1)/100, threshold(2)/100];
  Q = quantile (x, P, dim);
  L = Q(1);
  U = Q(2);
  C = (L + U) / 2;
endfunction


%!demo
%! A = [57 59 60 100 59 58 57 58 300 61 62 60 62 58 57];
%! TF = isoutlier (A, "mean")

%!demo
%! ## Use a moving detection method to detect local outliers in a sine wave
%!
%! x = -2*pi:0.1:2*pi;
%! A = sin(x);
%! A(47) = 0;
%! time = datenum (2023,1,1,0,0,0) + (1/24)*[0:length(x)-1] - 730485;
%! TF = isoutlier (A, "movmedian", 5*(1/24), "SamplePoints", time);
%! plot (time, A)
%! hold on
%! plot (time(TF), A(TF), "x")
%! datetick ('x', 20, 'keepticks')
%! legend ("Original Data", "Outlier Data")

%!demo
%! ## Locate an outlier in a vector of data and visualize the outlier
%!
%! x = 1:10;
%! A = [60 59 49 49 58 100 61 57 48 58];
%! [TF, L, U, C] = isoutlier (A);
%! plot (x, A);
%! hold on
%! plot (x(TF), A(TF), "x");
%! xlim ([1,10]);
%! line ([1,10], [L, L], "Linestyle", ":");
%! text (1.1, L-2, "Lower Threshold");
%! line ([1,10], [U, U], "Linestyle", ":");
%! text (1.1, U-2, "Upper Threshold");
%! line ([1,10], [C, C], "Linestyle", ":");
%! text (1.1, C-3, "Center Value");
%! legend ("Original Data", "Outlier Data");


## Output validation tests (checked against MATLAB)
%!test
%! A = [57 59 60 100 59 58 57 58 300 61 62 60 62 58 57];
%! assert (isoutlier (A, "mean"), logical([zeros(1,8) 1 zeros(1,6)]))
%! assert (isoutlier (A, "median"), ...
%! logical([zeros(1,3) 1 zeros(1,4) 1 zeros(1,6)]))

%!test
%! A = [57 59 60 100 59 58 57 58 300 61 62 60 62 58 57];
%! [TF, L, U, C] = isoutlier (A, "mean");
%! assert (L, -109.2459044922864, 1e-12)
%! assert (U, 264.9792378256198, 1e-12)
%! assert (C, 77.8666666666666, 1e-12)

%!test
%! A = [57 59 60 100 59 58 57 58 300 61 62 60 62 58 57];
%! [TF, L, U, C] = isoutlier (A, "median");
%! assert (L, 50.104386688966386, 1e-12)
%! assert (U, 67.895613311033610, 1e-12)
%! assert (C, 59)

%!test
%! A = magic(5) + diag(200*ones(1,5));
%! T = logical (eye (5));
%! assert (isoutlier (A, 2), T)

%!test
%! A = [57 59 60 100 59 58 57 58 300 61 62 60 62 58 57];
%! [TF, L, U, C] = isoutlier (A, "movmedian", 5);
%! l = [54.5522, 52.8283, 54.5522, 54.5522, 54.5522, 53.5522, 53.5522, ...
%!      53.5522, 47.6566, 56.5522, 57.5522, 56.5522, 51.1044, 52.3283, 53.5522];
%! u = [63.4478, 66.1717, 63.4478, 63.4478, 63.4478, 62.4478, 62.4478, ...
%!      62.4478, 74.3434, 65.4478, 66.4478, 65.4478, 68.8956, 65.6717, 62.4478];
%! c = [59, 59.5, 59, 59, 59, 58, 58, 58, 61, 61, 62, 61, 60, 59, 58];
%! assert (L, l, 1e-4)
%! assert (U, u, 1e-4)
%! assert (C, c)

%!test
%! A = [57 59 60 100 59 58 57 58 300 61 62 60 62 58 57];
%! [TF, L, U, C] = isoutlier (A, "movmedian", 5, "SamplePoints", [1:15]);
%! l = [54.5522, 52.8283, 54.5522, 54.5522, 54.5522, 53.5522, 53.5522, ...
%!      53.5522, 47.6566, 56.5522, 57.5522, 56.5522, 51.1044, 52.3283, 53.5522];
%! u = [63.4478, 66.1717, 63.4478, 63.4478, 63.4478, 62.4478, 62.4478, ...
%!      62.4478, 74.3434, 65.4478, 66.4478, 65.4478, 68.8956, 65.6717, 62.4478];
%! c = [59, 59.5, 59, 59, 59, 58, 58, 58, 61, 61, 62, 61, 60, 59, 58];
%! assert (L, l, 1e-4)
%! assert (U, u, 1e-4)
%! assert (C, c)

%!test
%! A = [57 59 60 100 59 58 57 58 300 61 62 60 62 58 57];
%! [TF, L, U, C] = isoutlier (A, "movmean", 5);
%! l = [54.0841,  6.8872, 11.5608, 12.1518, 11.0210, 10.0112, -218.2840, ...
%!      -217.2375, -215.1239, -213.4890, -211.3264, 55.5800, 52.9589, ...
%!      52.5979, 51.0627];
%! u = [63.2492, 131.1128, 122.4392, 122.2482, 122.5790, 122.7888, 431.0840, ...
%!      430.8375, 430.3239, 429.8890, 429.3264, 65.6200, 66.6411, 65.9021, ...
%!      66.9373];
%! c = [58.6667, 69, 67, 67.2, 66.8, 66.4, 106.4, 106.8, 107.6, 108.2, 109, ...
%!      60.6, 59.8, 59.25, 59];
%! assert (L, l, 1e-4)
%! assert (U, u, 1e-4)
%! assert (C, c, 1e-4)

%!test
%! A = [57 59 60 100 59 58 57 58 300 61 62 60 62 58 57];
%! [TF, L, U, C] = isoutlier (A, "movmean", 5, "SamplePoints", [1:15]);
%! l = [54.0841,  6.8872, 11.5608, 12.1518, 11.0210, 10.0112, -218.2840, ...
%!      -217.2375, -215.1239, -213.4890, -211.3264, 55.5800, 52.9589, ...
%!      52.5979, 51.0627];
%! u = [63.2492, 131.1128, 122.4392, 122.2482, 122.5790, 122.7888, 431.0840, ...
%!      430.8375, 430.3239, 429.8890, 429.3264, 65.6200, 66.6411, 65.9021, ...
%!      66.9373];
%! c = [58.6667, 69, 67, 67.2, 66.8, 66.4, 106.4, 106.8, 107.6, 108.2, 109, ...
%!      60.6, 59.8, 59.25, 59];
%! assert (L, l, 1e-4)
%! assert (U, u, 1e-4)
%! assert (C, c, 1e-4)

%!test
%! A = [57 59 60 100 59 58 57 58 300 61 62 60 62 58 57];
%! [TF, L, U, C] = isoutlier (A, "gesd");
%! assert (TF, logical ([0 0 0 1 0 0 0 0 1 0 0 0 0 0 0]))
%! assert (L, 34.235977035439944, 1e-12)
%! assert (U, 89.764022964560060, 1e-12)
%! assert (C, 62)

%!test
%! A = [57 59 60 100 59 58 57 58 300 61 62 60 62 58 57];
%! [TF, L, U, C] = isoutlier (A, "gesd", "ThresholdFactor", 0.01);
%! assert (TF, logical ([0 0 0 1 0 0 0 0 1 0 0 0 0 0 0]))
%! assert (L, 31.489256770616173, 1e-12)
%! assert (U, 92.510743229383820, 1e-12)
%! assert (C, 62)

%!test
%! A = [57 59 60 100 59 58 57 58 300 61 62 60 62 58 57];
%! [TF, L, U, C] = isoutlier (A, "gesd", "ThresholdFactor", 5e-10);
%! assert (TF, logical ([0 0 0 0 0 0 0 0 1 0 0 0 0 0 0]))
%! assert (L, 23.976664158788935, 1e-12)
%! assert (U, 100.02333584121110, 1e-12)
%! assert (C, 62)

%!test
%! A = [57 59 60 100 59 58 57 58 300 61 62 60 62 58 57];
%! [TF, L, U, C] = isoutlier (A, "grubbs");
%! assert (TF, logical ([0 0 0 1 0 0 0 0 1 0 0 0 0 0 0]))
%! assert (L, 54.642809574646606, 1e-12)
%! assert (U, 63.511036579199555, 1e-12)
%! assert (C, 59.076923076923080, 1e-12)

%!test
%! A = [57 59 60 100 59 58 57 58 300 61 62 60 62 58 57];
%! [TF, L, U, C] = isoutlier (A, "grubbs", "ThresholdFactor", 0.01);
%! assert (TF, logical ([0 0 0 1 0 0 0 0 1 0 0 0 0 0 0]))
%! assert (L, 54.216083184201850, 1e-12)
%! assert (U, 63.937762969644310, 1e-12)
%! assert (C, 59.076923076923080, 1e-12)

%!test
%! A = [57 59 60 100 59 58 57 58 300 61 62 60 62 58 57];
%! [TF, L, U, C] = isoutlier (A,  "percentiles", [10 90]);
%! assert (TF, logical ([0 0 0 0 0 0 0 0 1 0 0 0 0 0 0]))
%! assert (L, 57)
%! assert (U, 100)
%! assert (C, 78.5)

%!test
%! A = [57 59 60 100 59 58 57 58 300 61 62 60 62 58 57];
%! [TF, L, U, C] = isoutlier (A,  "percentiles", [20 80]);
%! assert (TF, logical ([1 0 0 1 0 0 1 0 1 0 0 0 0 0 1]))
%! assert (L, 57.5)
%! assert (U, 62)
%! assert (C, 59.75)

## Test input validation
%!shared A
%! A = [57 59 60 100 59 58 57 58 300 61 62 60 62 58 57];
%!error<isoutlier: WINDOW must be a positive scalar> ...
%! isoutlier (A, "movmedian", 0);
%!error<isoutlier: WINDOW must be a positive scalar> ...
%! isoutlier (A, "movmedian", []);
%!error<isoutlier: WINDOW must be a positive scalar> ...
%! isoutlier (A, "movmedian", [2 3 4]);
%!error<isoutlier: WINDOW must be a positive integer> ...
%! isoutlier (A, "movmedian", 1.4);
%!error<isoutlier: WINDOW must be a positive scalar> ...
%! isoutlier (A, "movmedian", [0 1]);
%!error<isoutlier: WINDOW must be a positive scalar> ...
%! isoutlier (A, "movmedian", [2 -1]);
%!error<isoutlier: WINDOW must be a positive scalar> ...
%! isoutlier (A, "movmedian", {2 3});
%!error<isoutlier: WINDOW must be a positive scalar> ...
%! isoutlier (A, "movmedian", "char");
%!
%!error<isoutlier: WINDOW must be a positive scalar> ...
%! isoutlier (A, "movmean", 0);
%!error<isoutlier: WINDOW must be a positive scalar> ...
%! isoutlier (A, "movmean", []);
%!error<isoutlier: WINDOW must be a positive scalar> ...
%! isoutlier (A, "movmean", [2 3 4]);
%!error<isoutlier: WINDOW must be a positive integer> ...
%! isoutlier (A, "movmean", 1.4);
%!error<isoutlier: WINDOW must be a positive scalar> ...
%! isoutlier (A, "movmean", [0 1]);
%!error<isoutlier: WINDOW must be a positive scalar> ...
%! isoutlier (A, "movmean", [2 -1]);
%!error<isoutlier: WINDOW must be a positive scalar> ...
%! isoutlier (A, "movmean", {2 3});
%!error<isoutlier: WINDOW must be a positive scalar> ...
%! isoutlier (A, "movmean", "char");
%!
%!error<isoutlier: THRESHOLD must be a two-element> ...
%! isoutlier (A, "percentiles", [-1 90]);
%!error<isoutlier: THRESHOLD must be a two-element> ...
%! isoutlier (A, "percentiles", [10 -90]);
%!error<isoutlier: THRESHOLD must be a two-element> ...
%! isoutlier (A, "percentiles", [90]);
%!error<isoutlier: THRESHOLD must be a two-element> ...
%! isoutlier (A, "percentiles", [90 20]);
%!error<isoutlier: THRESHOLD must be a two-element> ...
%! isoutlier (A, "percentiles", [90 20]);
%!error<isoutlier: THRESHOLD must be a two-element> ...
%! isoutlier (A, "percentiles", [10 20 90]);
%!error<isoutlier: THRESHOLD must be a two-element> ...
%! isoutlier (A, "percentiles", {10 90});
%!error<isoutlier: THRESHOLD must be a two-element> ...
%! isoutlier (A, "percentiles", "char");
%!
%!error<isoutlier: sample points must be a vector.> ...
%! isoutlier (A, "movmean", 5, "SamplePoints", ones(3,15));
%!error<isoutlier: sample points must be a vector.> ...
%! isoutlier (A, "movmean", 5, "SamplePoints", 15);
%!error<isoutlier: sample points must be unique.> ...
%! isoutlier (A, "movmean", 5, "SamplePoints", [1,1:14]);
%!error<isoutlier: sample points must be sorted.> ...
%! isoutlier (A, "movmean", 5, "SamplePoints", [2,1,3:15]);
%!error<isoutlier: sample points must have the same size> ...
%! isoutlier (A, "movmean", 5, "SamplePoints", [1:14]);
%!
%!error<isoutlier: threshold factor must be a nonnegative scalar.> ...
%! isoutlier (A, "movmean", 5, "ThresholdFactor", [1:14]);
%!error<isoutlier: threshold factor must be a nonnegative scalar.> ...
%! isoutlier (A, "movmean", 5, "ThresholdFactor", -1);
%!error<isoutlier: threshold factor must must be in> ...
%! isoutlier (A, "gesd", "ThresholdFactor", 3);
%!error<isoutlier: threshold factor must must be in> ...
%! isoutlier (A, "grubbs", "ThresholdFactor", 3);
%!
%!error<isoutlier: maximum outlier count must be a positive integer scalar.> ...
%! isoutlier (A, "movmean", 5, "MaxNumOutliers", [1:14]);
%!error<isoutlier: maximum outlier count must be a positive integer scalar.> ...
%! isoutlier (A, "movmean", 5, "MaxNumOutliers", -1);
%!error<isoutlier: maximum outlier count must be a positive integer scalar.> ...
%! isoutlier (A, "movmean", 5, "MaxNumOutliers", 0);
%!error<isoutlier: maximum outlier count must be a positive integer scalar.> ...
%! isoutlier (A, "movmean", 5, "MaxNumOutliers", 1.5);
%!
%!error<isoutlier: invalid input argument.> ...
%! isoutlier (A, {"movmean"}, 5, "SamplePoints", [1:15]);
%!error<isoutlier: invalid input argument.> isoutlier (A, {1});
%!error<isoutlier: invalid input argument.> isoutlier (A, true);
%!error<isoutlier: invalid input argument.> isoutlier (A, false);
%!error<isoutlier: DIM must be a positive integer scalar.> isoutlier (A, 0);
%!error<isoutlier: DIM must be a positive integer scalar.> isoutlier (A, [1 2]);
%!error<isoutlier: DIM must be a positive integer scalar.> isoutlier (A, -2);

