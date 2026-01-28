## Copyright (C) 2008 Francesco Potortì <pot@gnu.org>
## Copyright (C) 2024 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {statistics} {@var{D} =} pdist (@var{X})
## @deftypefnx {statistics} {@var{D} =} pdist (@var{X}, @var{Distance})
## @deftypefnx {statistics} {@var{D} =} pdist (@var{X}, @var{Distance}, @var{DistParameter})
##
## Return the distance between any two rows in @var{X}.
##
## @code{@var{D} = pdist (@var{X}} calculates the euclidean distance between
## pairs of observations in @var{X}.  @var{X} must be an @math{MxP} numeric
## matrix representing @math{M} points in @math{P}-dimensional space.  This
## function computes the pairwise distances returned in @var{D} as an
## @math{Mx(M-1)/P} row vector.  Use @code{@var{Z} = squareform (@var{D})} to
## convert the row vector @var{D} into a an @math{MxM} symmetric matrix @var{Z},
## where @qcode{@var{Z}(i,j)} corresponds to the pairwise distance between
## points @qcode{i} and @qcode{j}.
##
## @code{@var{D} = pdist (@var{X}, @var{Y}, @var{Distance})} returns the
## distance between pairs of observations in @var{X} using the metric specified
## by @var{Distance}, which can be any of the following options.
##
## @multitable @columnfractions 0.23 0.02 0.65
## @item @qcode{"euclidean"} @tab @tab Euclidean distance.
## @item @qcode{"squaredeuclidean"} @tab @tab Squared Euclidean distance.
## @item @qcode{"seuclidean"} @tab @tab standardized Euclidean distance.  Each
## coordinate difference between the rows in @var{X} and the query matrix
## @var{Y} is scaled by dividing by the corresponding element of the standard
## deviation computed from @var{X}.  A different scaling vector can be specified
## with the subsequent @var{DistParameter} input argument.
## @item @qcode{"mahalanobis"} @tab @tab Mahalanobis distance, computed using a
## positive definite covariance matrix.  A different covariance matrix can be
## specified with the subsequent @var{DistParameter} input argument.
## @item @qcode{"cityblock"} @tab @tab City block distance.
## @item @qcode{"minkowski"} @tab @tab Minkowski distance.  The default exponent
## is 2.  A different exponent can be specified with the subsequent
## @var{DistParameter} input argument.
## @item @qcode{"chebychev"} @tab @tab Chebychev distance (maximum coordinate
## difference).
## @item @qcode{"cosine"} @tab @tab One minus the cosine of the included angle
## between points (treated as vectors).
## @item @qcode{"correlation"} @tab @tab One minus the sample linear correlation
## between observations (treated as sequences of values).
## @item @qcode{"hamming"} @tab @tab Hamming distance, which is the percentage
## of coordinates that differ.
## @item @qcode{"jaccard"} @tab @tab One minus the Jaccard coefficient, which is
## the percentage of nonzero coordinates that differ.
## @item @qcode{"spearman"} @tab @tab One minus the sample Spearman's rank
## correlation between observations (treated as sequences of values).
## @item @var{@@distfun} @tab @tab Custom distance function handle.  A distance
## function of the form @code{function @var{D2} = distfun (@var{XI}, @var{YI})},
## where @var{XI} is a @math{1xP} vector containing a single observation in
## @math{P}-dimensional space, @var{YI} is an @math{NxP} matrix containing an
## arbitrary number of observations in the same @math{P}-dimensional space, and
## @var{D2} is an @math{NxP} vector of distances, where @qcode{(@var{D2}k)} is
## the distance between observations @var{XI} and @qcode{(@var{YI}k,:)}.
## @end multitable
##
## @code{@var{D} = pdist (@var{X}, @var{Y}, @var{Distance}, @var{DistParameter})}
## returns the distance using the metric specified by @var{Distance} and
## @var{DistParameter}.  The latter one can only be specified when the selected
## @var{Distance} is @qcode{"seuclidean"}, @qcode{"minkowski"}, and
## @qcode{"mahalanobis"}.
##
## @seealso{pdist2, squareform, linkage}
## @end deftypefn

function D = pdist (X, varargin)

  ## Check input data
  if (nargin < 1)
	  error ("pdist: too few input arguments.");
  endif
  if (! isnumeric (X) || isempty (X))
    error ("pdist: X must be a nonempty numeric matrix.");
  endif
  if (ndims (X) != 2)
	  error ("pdist: X must be a two-dimensional matrix.");
  endif
  if (rows (X) < 2)
    D = cast (zeros (1, 0), class (X));
    return;
  endif

  ## Add default values
  Distance = "euclidean";   # Distance metric
  DistParameter = [];       # Distance parameter

  ## Parse additional Distance metric and Distance parameter (if available)
  DMs = {"euclidean", "squaredeuclidean", "seuclidean", ...
         "mahalanobis", "cityblock", "minkowski", "chebychev", ...
         "cosine", "correlation", "hamming", "jaccard", "spearman"};
  if (numel (varargin) > 0)
    if (any (strcmpi (DMs, varargin{1})))
      Distance = tolower (varargin{1});
    elseif (is_function_handle (varargin{1}))
      Distance = varargin{1};
    else
      error ("pdist: invalid value for Distance input argument.");
    endif
  endif
  if (numel (varargin) > 1)
    if (isnumeric (varargin{2}))
      DistParameter = varargin{2};
    else
      error ("pdist: invalid value for DistParameter input argument.");
    endif
  endif

  ## Calculate selected distance
  ## Use loop-based computation to avoid O(n^2) memory for index arrays
  N = rows (X);
  num_pairs = N * (N - 1) / 2;
  D = zeros (1, num_pairs);

  ## Handle build-in distance metric
  if (ischar (Distance))
    switch (Distance)

      case "euclidean"
        ## Precomputation: none needed
        idx = 0;
        for i = 1:(N-1)
          Xi = X(i,:);
          for j = (i+1):N
            idx += 1;
            diffs = Xi - X(j,:);
            D(idx) = sqrt (sum (diffs .^ 2));
          endfor
        endfor

      case "squaredeuclidean"
        ## Precomputation: none needed
        idx = 0;
        for i = 1:(N-1)
          Xi = X(i,:);
          for j = (i+1):N
            idx += 1;
            diffs = Xi - X(j,:);
            D(idx) = sum (diffs .^ 2);
          endfor
        endfor

      case "seuclidean"
        ## Precomputation: compute or validate scale vector
        if (isempty (DistParameter))
          DistParameter = std (X, [], 1);
        else
          if (numel (DistParameter) != columns (X))
            error (strcat ("pdist: DistParameter for standardized", ...
                           " euclidean must be a vector of equal length", ...
                           " to the number of columns in X."));
          endif
          if (any (DistParameter < 0))
            error (strcat ("pdist: DistParameter for standardized", ...
                           " euclidean must be a nonnegative vector."));
          endif
        endif
        DistParameter(DistParameter == 0) = 1;  # fix constant variable
        ## Main loop
        idx = 0;
        for i = 1:(N-1)
          Xi = X(i,:);
          for j = (i+1):N
            idx += 1;
            diffs = (Xi - X(j,:)) ./ DistParameter;
            D(idx) = sqrt (sum (diffs .^ 2));
          endfor
        endfor

      case "mahalanobis"
        ## Precomputation: compute covariance and its inverse
        if (isempty (DistParameter))
          DistParameter = cov (X(! any (isnan (X), 2),:));
        else
          if (columns (DistParameter) != columns (X))
            error (strcat ("pdist: DistParameter for mahalanobis", ...
                           " distance must be a covariance matrix with", ...
                           " the same number of columns as X."));
          endif
          [~, p] = chol (DistParameter);
          if (p != 0)
            error (strcat ("pdist: covariance matrix for mahalanobis", ...
                           " distance must be symmetric and positive", ...
                           " definite."));
          endif
        endif
        ## Catch warning if matrix is close to singular or badly scaled.
        [DP_inv, rc] = inv (DistParameter);
        if (rc < eps)
          msg = sprintf (strcat ("pdist: matrix is close to", ...
                                 " singular or badly scaled.\n RCOND = ", ...
                                 " %e. Results may be inaccurate."), rc);
          warning (msg);
        endif
        ## Main loop
        idx = 0;
        for i = 1:(N-1)
          Xi = X(i,:);
          for j = (i+1):N
            idx += 1;
            dxx = Xi - X(j,:);
            D(idx) = sqrt (sum ((dxx * DP_inv) .* dxx));
          endfor
        endfor

      case "cityblock"
        ## Precomputation: none needed
        idx = 0;
        for i = 1:(N-1)
          Xi = X(i,:);
          for j = (i+1):N
            idx += 1;
            D(idx) = sum (abs (Xi - X(j,:)));
          endfor
        endfor

      case "minkowski"
        ## Precomputation: validate exponent parameter
        if (isempty (DistParameter))
          DistParameter = 2;
        else
          if (! (isnumeric (DistParameter) && isscalar (DistParameter)
                                           && DistParameter > 0))
            error (strcat ("pdist: DistParameter for minkowski distance", ...
                           " must be a positive scalar."));
          endif
        endif
        p_exp = DistParameter;
        p_inv = 1 / DistParameter;
        ## Main loop
        idx = 0;
        for i = 1:(N-1)
          Xi = X(i,:);
          for j = (i+1):N
            idx += 1;
            D(idx) = sum (abs (Xi - X(j,:)) .^ p_exp) .^ p_inv;
          endfor
        endfor

      case "chebychev"
        ## Precomputation: none needed
        idx = 0;
        for i = 1:(N-1)
          Xi = X(i,:);
          for j = (i+1):N
            idx += 1;
            D(idx) = max (abs (Xi - X(j,:)));
          endfor
        endfor

      case "cosine"
        ## Precomputation: compute inverse norms for all points
        sx = sum (X .^ 2, 2) .^ (-1 / 2);
        ## Main loop
        idx = 0;
        for i = 1:(N-1)
          Xi = X(i,:);
          sx_i = sx(i);
          for j = (i+1):N
            idx += 1;
            D(idx) = 1 - sum (Xi .* X(j,:)) * sx_i * sx(j);
          endfor
        endfor

      case "correlation"
        ## Precomputation: none needed (means depend on pairs)
        idx = 0;
        for i = 1:(N-1)
          Xi = X(i,:);
          for j = (i+1):N
            idx += 1;
            Xj = X(j,:);
            mXi = mean (Xi);
            mXj = mean (Xj);
            Xi_c = Xi - mXi;
            Xj_c = Xj - mXj;
            xy = sum (Xi_c .* Xj_c);
            xx = sqrt (sum (Xi_c .^ 2));
            yy = sqrt (sum (Xj_c .^ 2));
            D(idx) = 1 - (xy / (xx * yy));
          endfor
        endfor

      case "hamming"
        ## Precomputation: none needed
        idx = 0;
        for i = 1:(N-1)
          Xi = X(i,:);
          for j = (i+1):N
            idx += 1;
            D(idx) = mean (Xi != X(j,:));
          endfor
        endfor

      case "jaccard"
        ## Precomputation: none needed
        idx = 0;
        for i = 1:(N-1)
          Xi = X(i,:);
          for j = (i+1):N
            idx += 1;
            Xj = X(j,:);
            xy0 = (Xi != 0 | Xj != 0);
            D(idx) = sum ((Xi != Xj) & xy0) / sum (xy0);
          endfor
        endfor

      case "spearman"
        ## Precomputation: compute ranks for all rows
        rX = zeros (size (X));
        for i = 1:N
          rX(i,:) = tiedrank (X(i,:));
        endfor
        rM = (columns (X) + 1) / 2;
        ## Main loop
        idx = 0;
        for i = 1:(N-1)
          rXi = rX(i,:) - rM;
          for j = (i+1):N
            idx += 1;
            rXj = rX(j,:) - rM;
            xy = sum (rXi .* rXj);
            xx = sqrt (sum (rXi .^ 2));
            yy = sqrt (sum (rXj .^ 2));
            D(idx) = 1 - (xy / (xx * yy));
          endfor
        endfor

    endswitch
  endif

  ## Handle a function handle
  if (is_function_handle (Distance))
    ## Check the input output sizes of the user function
    D2 = [];
    try
      D2 = Distance (X(1,:), X([2:end],:));
    catch ME
      error ("pdist: invalid function handle for distance metric.");
    end_try_catch
    Xrows = N - 1;
    if (! isequal (size (D2), [Xrows, 1]))
      error ("pdist: custom distance function produces wrong output size.");
    endif
    ## Evaluate user defined distance metric function
    id_beg = 1;
    for r = 1:Xrows
      id_end = id_beg + (N - r) - 1;
      D(id_beg:id_end) = feval (Distance, X(r,:), X([r+1:end],:));
      id_beg = id_end + 1;
    endfor
  endif

endfunction

## Test output
%!shared xy, t, eucl, x
%! xy = [0 1; 0 2; 7 6; 5 6];
%! t = 1e-3;
%! eucl = @(v,m) sqrt(sumsq(repmat(v,rows(m),1)-m,2));
%! x = [1 2 3; 4 5 6; 7 8 9; 3 2 1];
%!assert (pdist (xy),                 [1.000 8.602 7.071 8.062 6.403 2.000], t);
%!assert (pdist (xy, eucl),           [1.000 8.602 7.071 8.062 6.403 2.000], t);
%!assert (pdist (xy, "euclidean"),    [1.000 8.602 7.071 8.062 6.403 2.000], t);
%!assert (pdist (xy, "seuclidean"),   [0.380 2.735 2.363 2.486 2.070 0.561], t);
%!assert (pdist (xy, "mahalanobis"),  [1.384 1.967 2.446 2.384 1.535 2.045], t);
%!assert (pdist (xy, "cityblock"),    [1.000 12.00 10.00 11.00 9.000 2.000], t);
%!assert (pdist (xy, "minkowski"),    [1.000 8.602 7.071 8.062 6.403 2.000], t);
%!assert (pdist (xy, "minkowski", 3), [1.000 7.763 6.299 7.410 5.738 2.000], t);
%!assert (pdist (xy, "cosine"),       [0.000 0.349 0.231 0.349 0.231 0.013], t);
%!assert (pdist (xy, "correlation"),  [0.000 2.000 0.000 2.000 0.000 2.000], t);
%!assert (pdist (xy, "spearman"),     [0.000 2.000 0.000 2.000 0.000 2.000], t);
%!assert (pdist (xy, "hamming"),      [0.500 1.000 1.000 1.000 1.000 0.500], t);
%!assert (pdist (xy, "jaccard"),      [1.000 1.000 1.000 1.000 1.000 0.500], t);
%!assert (pdist (xy, "chebychev"),    [1.000 7.000 5.000 7.000 5.000 2.000], t);
%!assert (pdist (x), [5.1962, 10.3923, 2.8284, 5.1962, 5.9161, 10.7703], 1e-4);
%!assert (pdist (x, "euclidean"), ...
%!        [5.1962, 10.3923, 2.8284, 5.1962, 5.9161, 10.7703], 1e-4);
%!assert (pdist (x, eucl), ...
%!        [5.1962, 10.3923, 2.8284, 5.1962, 5.9161, 10.7703], 1e-4);
%!assert (pdist (x, "squaredeuclidean"), [27, 108, 8, 27, 35, 116]);
%!assert (pdist (x, "seuclidean"), ...
%!        [1.8071, 3.6142, 0.9831, 1.8071, 1.8143, 3.4854], 1e-4);
%!warning<pdist: matrix is close to singular> ...
%! pdist (x, "mahalanobis");
%!assert (pdist (x, "cityblock"), [9, 18, 4, 9, 9, 18]);
%!assert (pdist (x, "minkowski"), ...
%!        [5.1962, 10.3923, 2.8284, 5.1962, 5.9161, 10.7703], 1e-4);
%!assert (pdist (x, "minkowski", 3), ...
%!        [4.3267, 8.6535, 2.5198, 4.3267, 5.3485, 9.2521], 1e-4);
%!assert (pdist (x, "cosine"), ...
%!        [0.0254, 0.0406, 0.2857, 0.0018, 0.1472, 0.1173], 1e-4);
%!assert (pdist (x, "correlation"), [0, 0, 2, 0, 2, 2], 1e-14);
%!assert (pdist (x, "spearman"), [0, 0, 2, 0, 2, 2], 1e-14);
%!assert (pdist (x, "hamming"), [1, 1, 2/3, 1, 1, 1]);
%!assert (pdist (x, "jaccard"), [1, 1, 2/3, 1, 1, 1]);
%!assert (pdist (x, "chebychev"), [3, 6, 2, 3, 5, 8]);
