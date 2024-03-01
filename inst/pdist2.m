## Copyright (C) 2014-2019 Piotr Dollar <pdollar@gmail.com>
## Copyright (C) 2024 Andreas Bertsatos <abertsatos@biol.uoa.gr>
##
## This file is part of the statistics package for GNU Octave.
##
## This program is free software; you can redistribute it and/or
## modify it under the terms of the GNU General Public License as
## published by the Free Software Foundation; either version 3 of the
## License, or (at your option) any later version.
##
## This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
## General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, see
## <http:##www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn  {statistics} {@var{D} =} pdist2 (@var{X}, @var{Y})
## @deftypefnx {statistics} {@var{D} =} pdist2 (@var{X}, @var{Y}, @var{Distance})
## @deftypefnx {statistics} {@var{D} =} pdist2 (@var{X}, @var{Y}, @var{Distance}, @var{DistParameter})
## @deftypefnx {statistics} {@var{D} =} pdist2 (@dots{}, @var{Name}, @var{Value})
## @deftypefnx {statistics} {[@var{D}, @var{I}] =} pdist2 (@dots{}, @var{Name}, @var{Value})
##
## Compute pairwise distance between two sets of vectors.
##
## @code{@var{D} = pdist2 (@var{X}, @var{Y})} calculates the euclidean distance
## between each pair of observations in @var{X} and @var{Y}.  Let @var{X} be an
## @math{MxP} matrix representing @math{M} points in @math{P}-dimensional space
## and @var{Y} be an @math{NxP} matrix representing another set of points in the
## same space.  This function computes the @math{MxN} distance matrix @var{D},
## where @qcode{@var{D}(i,j)} is the distance between @qcode{@var{X}(i,:)} and
## @qcode{@var{Y}(j,:)}.
##
## @code{@var{D} = pdist2 (@var{X}, @var{Y}, @var{Distance})} returns the
## distance between each pair of observations in @var{X} and @var{Y} using the
## metric specified by @var{Distance}, which can be any of the following options.
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
## @code{@var{D} = pdist2 (@var{X}, @var{Y}, @var{Distance}, @var{DistParameter})}
## returns the distance using the metric specified by @var{Distance} and
## @var{DistParameter}.  The latter one can only be specified when the selected
## @var{Distance} is @qcode{"seuclidean"}, @qcode{"minkowski"}, and
## @qcode{"mahalanobis"}.
##
## @code{@var{D} = pdist2 (@dots{}, @var{Name}, @var{Value})}  for any previous
## arguments, modifies the computation using @var{Name}-@var{Value} parameters.
## @itemize
## @item
## @code{@var{D} = pdist2 (@var{X}, @var{Y}, @var{Distance}, @qcode{"Smallest"},
## @var{K})} computes the distance using the metric specified by
## @var{Distance} and returns the @var{K} smallest pairwise distances to
## observations in @var{X} for each observation in @var{Y} in ascending order.
## @item
## @code{@var{D} = pdist2 (@var{X}, @var{Y}, @var{Distance}, @var{DistParameter},
## @qcode{"Largest"}, @var{K})} computes the distance using the metric specified
## by @var{Distance} and @var{DistParameter} and returns the @var{K} largest
## pairwise distances in descending order.
## @end itemize
##
## @code{[@var{D}, @var{I}] = pdist2 (@dots{}, @var{Name}, @var{Value})} also
## returns the matrix @var{I}, which contains the indices of the observations in
## @var{X} corresponding to the distances in @var{D}.  You must specify either
## @qcode{"Smallest"} or @qcode{"Largest"} as an optional @var{Name}-@var{Value}
## pair pair argument to compute the second output argument.
##
## @seealso{pdist, knnsearch, rangesearch}
## @end deftypefn

function [D, I] = pdist2 (X, Y, varargin)

  ## Check input data
  if (nargin < 2)
	  error ("pdist2: too few input arguments.");
  endif
  if (size (X, 2) != size (Y, 2))
	  error ("pdist2: X and Y must have equal number of columns.");
  endif
  if (ndims (X) != 2 || ndims (Y) != 2)
    error ("pdist2: X and Y must be 2 dimensional matrices.");
  endif

  ## Add default values
  Distance = "euclidean";   # Distance metric
  DistParameter = [];       # Distance parameter
  SortOrder = [];           # Flag for sorting distances to find

  ## Parse additional Distance metric and Distance parameter (if available)
  DMs = {"euclidean", "squaredeuclidean", "seuclidean", ...
         "mahalanobis", "cityblock", "minkowski", "chebychev", ...
         "cosine", "correlation", "hamming", "jaccard", "spearman"};
  if (numel (varargin) > 0)
    if (any (strcmpi (DMs, varargin{1})))
      Distance = tolower (varargin{1});
      varargin(1) = [];
      if (numel (varargin) > 0)
        if (isnumeric (varargin{1}))
          DistParameter = varargin{1};
          varargin(1) = [];
        endif
      endif
    elseif (is_function_handle (varargin{1}))
      Distance = varargin{1};
      varargin(1) = [];
      if (numel (varargin) > 0)
        if (isnumeric (varargin{1}))
          DistParameter = varargin{1};
          varargin(1) = [];
        endif
      endif
    endif
  endif

  ## Parse additional parameters in Name/Value pairs
  parcount = 0;
  while (numel (varargin) > 0)
    if (numel (varargin) < 2)
      error ("pdist2: missing value in optional name/value paired arguments.");
    endif
    switch (tolower (varargin{1}))
      case "smallest"
        SortOrder = "ascend";
        K = varargin{2};
        parcount += 1;
      case "largest"
        SortOrder = "descend";
        K = varargin{2};
        parcount += 1;
      otherwise
        error ("pdist2: invalid NAME in optional pairs of arguments.");
    endswitch
    varargin(1:2) = [];
  endwhile

  ## Check additional arguments
  if (parcount > 1)
    error ("pdist2: you can only use either Smallest or Largest.");
  endif
  if (isempty (SortOrder) && nargout > 1)
    error (strcat (["pdist2: Smallest or Largest must be specified"], ...
                   [" to compute second output."]));
  endif

  ## Calculate selected distance
  [ix, iy] = meshgrid (1:size (X, 1), 1:size (Y, 1));

  ## Handle build-in distance metric
  if (ischar (Distance))
    switch (Distance)
      case "euclidean"
        D = sqrt (sum ((X(ix(:),:) - Y(iy(:),:)) .^ 2, 2));

      case "squaredeuclidean"
        D = sum ((X(ix(:),:) - Y(iy(:),:)) .^ 2, 2);

      case "seuclidean"
        if (isempty (DistParameter))
          DistParameter = std (X, [], 1);
        else
          if (numel (DistParameter) != columns (X))
            error (strcat (["pdist2: DistParameter for standardized"], ...
                           [" euclidean must be a vector of equal length"], ...
                           [" to the number of columns in X."]));
          endif
          if (any (DistParameter < 0))
            error (strcat (["pdist2: DistParameter for standardized"], ...
                           [" euclidean must be a nonnegative vector."]));
          endif
        endif
        DistParameter(DistParameter == 0) = 1;  # fix constant variable
        D = sqrt (sum (((X(ix(:),:) - Y(iy(:),:)) ./ DistParameter) .^ 2, 2));

      case "mahalanobis"
        if (isempty (DistParameter))
          DistParameter = cov (X(! any (isnan (X), 2),:));
        else
          if (columns (DistParameter) != columns (X))
            error (strcat (["pdist2: DistParameter for mahalanobis"], ...
                           [" distance must be a covariance matrix with"], ...
                           [" the same number of columns as X."]));
          endif
          [~, p] = chol (DistParameter);
          if (p != 0)
            error (strcat (["pdist2: covariance matrix for mahalanobis"], ...
                           [" distance must be symmetric and positive"], ...
                           [" definite."]));
          endif
        endif
        ## Catch warning if matrix is close to singular or badly scaled.
        [DP_inv, rc] = inv (DistParameter);
        if (rc < eps)
          msg = sprintf (strcat (["pdist2: matrix is close to"], ...
                                 [" singular or badly scaled.\n RCOND = "], ...
                                 [" %e. Results may be inaccurate."]), rc);
          warning (msg);
        endif
        dxy = X(ix(:),:) - Y(iy(:),:);
        D   = sqrt (sum ((dxy * DP_inv) .* dxy, 2));

      case "cityblock"
        D = sum (abs (X(ix(:),:) - Y(iy(:),:)), 2);

      case "minkowski"
        if (isempty (DistParameter))
          DistParameter = 2;
        else
          if (! (isnumeric (DistParameter) && isscalar (DistParameter)
                                           && DistParameter > 0))
            error (strcat (["pdist2: DistParameter for minkowski distance"],...
                           [" must be a positive scalar."]));
          endif
        endif
        D = sum (abs (X(ix(:),:) - Y(iy(:),:)) .^ DistParameter, 2) .^ ...
                (1 / DistParameter);

      case "chebychev"
        D = max (abs (X(ix(:),:) - Y(iy(:),:)), [], 2);

      case "cosine"
        sx = sum (X .^ 2, 2) .^ (-1 / 2);
        sy = sum (Y .^ 2, 2) .^ (-1 / 2);
        D  = 1 - sum (X(ix(:),:) .* Y(iy(:),:), 2) .* sx(ix(:)) .* sy(iy(:));

      case "correlation"
        mX = mean (X(ix(:),:), 2);
        mY = mean (Y(iy(:),:), 2);
        xy = sum ((X(ix(:),:) - mX) .* (Y(iy(:),:) - mY), 2);
        xx = sqrt (sum ((X(ix(:),:) - mX) .* (X(ix(:),:) - mX), 2));
        yy = sqrt (sum ((Y(iy(:),:) - mY) .* (Y(iy(:),:) - mY), 2));
        D = 1 - (xy ./ (xx .* yy));

      case "hamming"
        D = mean (abs (X(ix(:),:) != Y(iy(:),:)), 2);

      case "jaccard"
        xy0 = (X(ix(:),:) != 0 | Y(iy(:),:) != 0);
        D = sum ((X(ix(:),:) != Y(iy(:),:)) & xy0, 2) ./ sum (xy0, 2);

      case "spearman"
        for i = 1:size (X, 1)
          rX(i,:) = tiedrank (X(i,:));
        endfor
        for i = 1:size (Y, 1)
          rY(i,:) = tiedrank (Y(i,:));
        endfor
        rM = (size (X, 2) + 1) / 2;
        xy = sum ((rX(ix(:),:) - rM) .* (rY(iy(:),:) - rM), 2);
        xx = sqrt (sum ((rX(ix(:),:) - rM) .* (rX(ix(:),:) - rM), 2));
        yy = sqrt (sum ((rY(iy(:),:) - rM) .* (rY(iy(:),:) - rM), 2));
        D = 1 - (xy ./ (xx .* yy));

    endswitch
  endif

  ## Handle a function handle
  if (is_function_handle (Distance))
    ## Check the input output sizes of the user function
    D2 = [];
    try
      D2 = Distance (X(1,:), Y);
    catch ME
      error ("pdist2: invalid function handle for distance metric.");
    end_try_catch
    Yrows = rows (Y);
    if (! isequal (size (D2), [Yrows, 1]))
      error ("pdist2: custom distance function produces wrong output size.");
    endif
    ## Evaluate user defined distance metric function
    Yrows = rows (Y);
    D = zeros (numel (ix), 1);
    id_beg = 1;
    for r = 1:rows (X)
      id_end = id_beg + Yrows - 1;
      D(id_beg:id_end) = feval (Distance, X(r,:), Y);
      id_beg = id_end + 1;
    endfor
  endif

  ## From vector to matrix
  D = reshape (D, size (Y, 1), size (X, 1))';

  if (nargout > 1)
    [D, I] = sort (D', 2, SortOrder);
    K = min (size (D, 2), K);   # fix max K to avoid out of bound error
    D = D(:,1:K)';
    I = I(:,1:K)';
  endif

endfunction


## Test output
%!shared x, y, xx
%! x = [1, 1, 1; 2, 2, 2; 3, 3, 3];
%! y = [0, 0, 0; 1, 2, 3; 0, 2, 4; 4, 7, 1];
%! xx = [1 2 3; 4 5 6; 7 8 9; 3 2 1];
%!test
%! d = sqrt([3, 5, 11, 45; 12, 2, 8, 30; 27, 5, 11, 21]);
%! assert (pdist2 (x, y), d);
%!test
%! d = [5.1962, 2.2361, 3.3166, 6.7082; ...
%!      3.4641, 2.2361, 3.3166, 5.4772];
%! i = [3, 1, 1, 1; 2, 3, 3, 2];
%! [D, I] = pdist2 (x, y, "euclidean", "largest", 2);
%! assert ({D, I}, {d, i}, 1e-4);
%!test
%! d = [1.7321, 1.4142, 2.8284, 4.5826; ...
%!      3.4641, 2.2361, 3.3166, 5.4772];
%! i = [1, 2, 2, 3;2, 1, 1, 2];
%! [D, I] = pdist2 (x, y, "euclidean", "smallest", 2);
%! assert ({D, I}, {d, i}, 1e-4);
%!test
%! yy = [1 2 3;5 6 7;9 5 1];
%! d = [0, 6.1644, 5.3852; 1.4142, 6.9282, 8.7750; ...
%!      3.7417, 7.0711, 9.9499; 6.1644, 10.4881, 10.3441];
%! i = [2, 4, 4; 3, 2, 2; 1, 3, 3; 4, 1, 1];
%! [D, I] = pdist2 (y, yy, "euclidean", "smallest", 4);
%! assert ({D, I}, {d, i}, 1e-4);
%!test
%! yy = [1 2 3;5 6 7;9 5 1];
%! d = [0, 38, 29; 2, 48, 77; 14, 50, 99; 38, 110, 107];
%! i = [2, 4, 4; 3, 2, 2; 1, 3, 3; 4, 1, 1];
%! [D, I] = pdist2 (y, yy, "squaredeuclidean", "smallest", 4);
%! assert ({D, I}, {d, i}, 1e-4);
%!test
%! yy = [1 2 3;5 6 7;9 5 1];
%! d = [0, 3.3256, 2.7249; 0.7610, 3.3453, 4.4799; ...
%!      1.8514, 3.3869, 5.0703; 2.5525, 5.0709, 5.1297];
%! i = [2, 2, 4; 3, 4, 2; 1, 3, 1; 4, 1, 3];
%! [D, I] = pdist2 (y, yy, "seuclidean", "smallest", 4);
%! assert ({D, I}, {d, i}, 1e-4);
%!test
%! d = [2.1213, 4.2426, 6.3640; 1.2247, 2.4495, 4.4159; ...
%!      3.2404, 4.8990, 6.8191; 2.7386, 4.2426, 6.1237];
%! assert (pdist2 (y, x, "mahalanobis"), d, 1e-4);
%!test
%! xx = [1, 3, 4; 3, 5, 4; 8, 7, 6];
%! d = [1.3053, 1.8257, 15.0499; 1.3053, 3.3665, 16.5680];
%! i = [2, 2, 2; 3, 4, 4];
%! [D, I] = pdist2 (y, xx, "mahalanobis", "smallest", 2);
%! assert ({D, I}, {d, i}, 1e-4);
%!test
%! d = [2.5240, 4.1633, 17.3638; 2.0905, 3.9158, 17.0147];
%! i = [1, 1, 3; 4, 3, 1];
%! [D, I] = pdist2 (y, xx, "mahalanobis", "largest", 2);
%! assert ({D, I}, {d, i}, 1e-4);
%!test
%! yy = [1 2 3;5 6 7;9 5 1];
%! d = [0, 8.4853, 18.0416; 2.4495, 10.0995, 19.4808; ...
%!      2.4495, 10.6771, 19.7104; 2.4495, 10.6771, 20.4573];
%! i = [2, 2, 2; 1, 4, 4; 4, 1, 1; 3, 3, 3];
%! [D, I] = pdist2 (y, yy, "mahalanobis", "smallest", 4);
%! assert ({D, I}, {d, i}, 1e-4);
%!test
%! d = [3, 3, 5, 9; 6, 2, 4, 8; 9, 3, 5, 7];
%! assert (pdist2 (x, y, "cityblock"), d);
%!test
%! d = [1, 2, 3, 6; 2, 1, 2, 5; 3, 2, 3, 4];
%! assert (pdist2 (x, y, "chebychev"), d);
%!test
%! d = repmat ([NaN, 0.0742, 0.2254, 0.1472], [3, 1]);
%! assert (pdist2 (x, y, "cosine"), d, 1e-4);
%!test
%! yy = [1 2 3;5 6 7;9 5 1];
%! d = [0, 0, 0.5; 0, 0, 2; 1.5, 1.5, 2; NaN, NaN, NaN];
%! i = [2, 2, 4; 3, 3, 2; 4, 4, 3; 1, 1, 1];
%! [D, I] = pdist2 (y, yy, "correlation", "smallest", 4);
%! assert ({D, I}, {d, i}, eps);
%! [D, I] = pdist2 (y, yy, "spearman", "smallest", 4);
%! assert ({D, I}, {d, i}, eps);
%!test
%! d = [1, 2/3, 1, 1; 1, 2/3, 1, 1; 1, 2/3, 2/3, 2/3];
%! i = [1, 1, 1, 2; 2, 2, 3, 3; 3, 3, 2, 1];
%! [D, I] = pdist2 (x, y, "hamming", "largest", 4);
%! assert ({D, I}, {d, i}, eps);
%! [D, I] = pdist2 (x, y, "jaccard", "largest", 4);
%! assert ({D, I}, {d, i}, eps);
%!test
%! xx = [1, 2, 3, 4; 2, 3, 4, 5; 3, 4, 5, 6];
%! yy = [1, 2, 2, 3; 2, 3, 3, 4];
%! [D, I] = pdist2 (x, y, "euclidean", "Smallest", 4);
%! eucldist = @(v,m) sqrt(sumsq(repmat(v,rows(m),1)-m,2));
%! [d, i] = pdist2 (x, y, eucldist, "Smallest", 4);
%! assert ({D, I}, {d, i});
%!warning<pdist2: matrix is close to singular> ...
%! pdist2 (xx, xx, "mahalanobis");

## Test input validation
%!error<pdist2: too few input arguments.> pdist2 (1)
%!error<pdist2: X and Y must have equal number of columns.> ...
%! pdist2 (ones (4, 5), ones (4))
%!error<pdist2: X and Y must be 2 dimensional matrices.> ...
%! pdist2 (ones (4, 2, 3), ones (3, 2))
%!error<pdist2: missing value in optional name/value paired arguments.> ...
%! pdist2 (ones (3), ones (3), "euclidean", "Largest")
%!error<pdist2: missing value in optional name/value paired arguments.> ...
%! pdist2 (ones (3), ones (3), "minkowski", 3, "Largest")
%!error<pdist2: invalid NAME in optional pairs of arguments.> ...
%! pdist2 (ones (3), ones (3), "minkowski", 3, "large", 4)
%!error<pdist2: you can only use either Smallest or Largest.> ...
%! pdist2 (ones (3), ones (3), "minkowski", 3, "Largest", 4, "smallest", 5)
%!error<pdist2: Smallest or Largest must be specified to compute second output.> ...
%! [d, i] = pdist2(ones (3), ones (3), "minkowski", 3)
%!error<pdist2: DistParameter for standardized euclidean must be a vector of> ...
%! pdist2 (ones (3), ones (3), "seuclidean", 3)
%!error<pdist2: DistParameter for standardized euclidean must be a nonnegative> ...
%! pdist2 (ones (3), ones (3), "seuclidean", [1, -1, 3])
%!error<pdist2: DistParameter for mahalanobis distance must be a covariance> ...
%! pdist2 (ones (3), eye (3), "mahalanobis", eye(2))
%!error<pdist2: covariance matrix for mahalanobis distance must be symmetric> ...
%! pdist2 (ones (3), eye (3), "mahalanobis", ones(3))
%!error<pdist2: DistParameter for minkowski distance must be a positive scalar.> ...
%! pdist2 (ones (3), eye (3), "minkowski", 0)
%!error<pdist2: DistParameter for minkowski distance must be a positive scalar.> ...
%! pdist2 (ones (3), eye (3), "minkowski", -5)
%!error<pdist2: DistParameter for minkowski distance must be a positive scalar.> ...
%! pdist2 (ones (3), eye (3), "minkowski", [1, 2])
%!error<pdist2: invalid function handle for distance metric.> ...
%! pdist2 (ones (3), ones (3), @(v,m) sqrt(repmat(v,rows(m),1)-m,2))
%!error<pdist2: custom distance function produces wrong output size.> ...
%! pdist2 (ones (3), ones (3), @(v,m) sqrt(sum(sumsq(repmat(v,rows(m),1)-m,2))))
