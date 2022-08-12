## Copyright (C) 2022 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn {Function File} @var{d} = manova1 (@var{x}, @var{group})
## @deftypefnx {Function File} @var{d} = manova1 (@var{x}, @var{group}, @var{alpha})
## @deftypefnx {Function File} [@var{d}, @var{p}] = manova1 (@dots{})
## @deftypefnx {Function File} [@var{d}, @var{p}, @var{stats}] = manova1 (@dots{})
##
## One-way multivariate analysis of variance (MANOVA).
##
## @code{@var{d} = manova1 (@var{x}, @var{group}, @var{alpha})} performs a
## one-way MANOVA for comparing the mean vectors of two or more groups of
## multivariate data.
##
## @var{x} is a matrix with each row representing a multivariate observation,
## and each column representing a variable.
##
## @var{group} is a numeric vector, string array, or cell array of strings with
## the same number of rows as @var{x}.  @var{x} values are in the same group if
## they correspond to the same value of GROUP.
##
## @var{alpha} is the scalar significance level and is 0.05 by default.
##
## @var{d} is an estimate of the dimension of the group means.  It is the
## smallest dimension such that a test of the hypothesis that the means lie on
## a space of that dimension is not rejected.  If @var{d} = 0 for example, we
## cannot reject the hypothesis that the means are the same.  If @var{d} = 1, we
## reject the hypothesis that the means are the same but we cannot reject the
## hypothesis that they lie on a line.
##
## @code{[@var{d}, @var{p}] = manova1 (@dots{})} returns P, a vector of p-values
## for testing the null hypothesis that the mean vectors of the groups lie on
## various dimensions.  P(1) is the p-value for a test of dimension 0, P(2) for
## dimension 1, etc.
##
## @code{[@var{d}, @var{p}, @var{stats}] = manova1 (@dots{})} returns a STATS
## structure with the following fields:
##
## @multitable @columnfractions 0.05 0.2 0.75
## @item @tab "W" @tab within-group sum of squares and products matrix
## @item @tab "B" @tab between-group sum of squares and products matrix
## @item @tab "T" @tab total sum of squares and products matrix
## @item @tab "dfW" @tab degrees of freedom for WSSP matrix
## @item @tab "dfB" @tab degrees of freedom for BSSP matrix
## @item @tab "dfT" @tab degrees of freedom for TSSP matrix
## @item @tab "lambda" @tab value of Wilk's lambda (the test statistic)
## @item @tab "chisq" @tab transformation of lambda to a chi-square distribution
## @item @tab "chisqdf" @tab degrees of freedom for chisq
## @item @tab "eigenval" @tab eigenvalues of (WSSP^-1) * BSSP
## @item @tab "eigenvec" @tab eigenvectors of (WSSP^-1) * BSSP; these are the
## coefficients for canonical variables, and they are scaled so the within-group
## variance of C is 1
## @item @tab "canon" @tab canonical variables, equal to XC*eigenvec, where XC
## is X with columns centered by subtracting their means
## @item @tab "mdist" @tab Mahalanobis distance from each point to its group mean
## @item @tab "gmdist" @tab Mahalanobis distances between each pair of group means
## @end multitable
##
## The canonical variables C have the property that C(:,1) is the linear
## combination of the @var{x} columns that has the maximum separation between
## groups, C(:,2) has the maximum separation subject to it being orthogonal to
## C(:,1), and so on.
##
## @end deftypefn
 
function [d, p, stats] = manova1 (x, group, alpha)
  
  ## Check input arguments
  narginchk(2,3)
  nargoutchk(1,3)
  
  ## Validate alpha value if parsed or add default
  if (nargin > 2)
    if (length (alpha) > 1 || ! isreal (alpha))
      error ("manova1: Alpha must be a real scalar.");
    elseif (alpha <= 0 || alpha >= 1)
      error ("manova1: Alpha must be in the range (0,1).");
    endif
  else
    alpha = 0.05;
  endif
  
  ## Convert group to cell array from character array
  if (ischar (group))
    group = cellstr (group);
  endif
  ## Make group a column
  if (size (group, 1) == 1)
    group = group';
  endif

  ## Check for equal size in samples between groups and data
  if (size (group, 1) != size (x, 1))
    error("manova1: Samples in X and groups mismatch.");
  endif

  ## Remove samples (rows) in X and GROUP if there are missing values in X
  no_nan = (sum (isnan (x), 2) == 0);
  x = x(no_nan, :);
  group = group(no_nan, :);
  is_nan = ! no_nan;

  ## Get group names and indices
  [group_idx, group_names] = grp2idx (group);
  ngroups = length (group_names);
  
  ## Remove NaN values from updated GROUP
  no_nan = ! isnan (group_idx);
  if (! all (no_nan))
    group_idx = group_idx(no_nan);
    x = x(no_nan,: );
    is_nan(! is_nan) = ! no_nan;
  endif
  
  ## Get number of samples and variables
  [nsample, nvar] = size(x);
  realgroups = ismember(1:ngroups, group_idx);
  nrgroups = sum(realgroups);

  ## Calculate Total Sum of Squares and Products matrix
  xm = mean (x);
  x = x - xm;
  TSSP = x' * x;
  ## Calculate Within-samples Sum of Squares and Products matrix
  WSSP = zeros (size (TSSP));
  for j = 1:ngroups
     row = find (group_idx == j);
     ## Only meaningful for groups with more than one samples
     if (length (row) > 1)
        group_x = x(row, :);
        group_x = group_x - mean (group_x);
        WSSP = WSSP + group_x' * group_x;
     endif
  endfor
  ## Calculate Between-samples Sum of Squares and Products matrix
  BSSP = TSSP - WSSP;

  ## Instead of simply computing `eig (BSSP / WSSP)` we use Matlab's technique
  ## with chol to insure v' * WSSP * v = I is met
  [R, p] = chol (WSSP);
  if (p > 0)
     error("manova1: Cannot factorize WSSP.");
  endif
  S = R' \ BSSP / R;
  ## Remove asymmetry caused by roundoff
  S = (S + S') / 2;
  [vv, ed] = eig (S);
  v = R \ vv;
  ## Sort in descending order
  [e,ei] = sort (diag (ed));
  ## Check for valid eigevalues
  if (min(e) <= -1)
     error ("manova1: wrong value in eigenvector: singular sum of squares.");
  endif
  
  ## Compute Barlett's statistic for each dimension
  dims = 0:(min (nrgroups - 1, nvar) - 1);
  lambda = flipud (1 ./ cumprod (e + 1));
  lambda = lambda(1 + dims);
  chistat = -(nsample - 1 - (nrgroups + nvar) / 2) .* log (lambda);
  chisqdf = ((nvar - dims) .* (nrgroups - 1 - dims))';
  pp = 1 - chi2cdf (chistat, chisqdf);

  ## Get dimension where we can reject the null hypothesis
  d = dims(pp>alpha);
  if (length(d) > 0)
     d = d(1);
  else
     d = max(dims) + 1;
  end

  ## Create extra outputs as necessary
  if (nargout > 1)
    p = pp;
  endif
  if (nargout > 2)
    stats.W = WSSP;
    stats.B = BSSP;
    stats.T = TSSP;
    stats.dfW = nsample - nrgroups;
    stats.dfB = nrgroups - 1;
    stats.dfT = nsample - 1;
    stats.lambda = lambda;
    stats.chisq = chistat;
    stats.chisqdf = chisqdf;
    ## Reorder to increasing
    stats.eigenval = flipud(e);

    ## Flip so that it is in order of increasing eigenvalues
    v = v(:, flipud (ei));
    ## Re-scale eigenvectors so the within-group variance is 1
    vs = diag((v' * WSSP * v))' ./ (nsample - nrgroups);
    vs(vs<=0) = 1;
    v = v ./ repmat(sqrt(vs), size(v,1), 1);

    ## Flip sign so that the average element is positive
    j = (sum(v) < 0);
    v(:,j) = -v(:,j);
    stats.eigenvec = v;
    canon = x*v;
    if (any(is_nan))
      tmp(~is_nan,:) = canon;
      tmp(is_nan,:) = NaN;
      stats.canon = tmp;
    else
      stats.canon = canon;
    endif

    ## Compute Mahalanobis distances from points to group means
    gmean = nan (ngroups, size (canon, 2));
    gmean(realgroups,:) = grpstats (canon, group_idx);
    mdist = sum ((canon - gmean(group_idx,:)) .^ 2, 2);
    if (any (is_nan))
      stats.mdist(! is_nan) = mdist;
      stats.mdist(is_nan) = NaN;
    else
      stats.mdist = mdist;
    endif
    ## Compute Mahalanobis distances between group means
    stats.gmdist = squareform (pdist (gmean)) .^ 2;
    stats.group_names = group_names;
  endif
endfunction

%!demo
%! load carbig
%! [d,p] = manova1([MPG, Acceleration, Weight, Displacement], Origin)

%!test
%! load carbig
%! [d,p] = manova1([MPG, Acceleration, Weight, Displacement], Origin);
%! assert (d, 3);
%! assert (p, [0, 3.140583347827075e-07, 0.007510999577743149, ...
%!             0.1934100745898493]', [1e-14, 1e-14, 1e-14, 1e-14]');

%!test
%! load carbig
%! [d,p] = manova1([MPG, Acceleration, Weight], Origin);
%! assert (d, 2);
%! assert (p, [0, 0.00516082975137544, 0.1206528056514453]', ...
%!            [1e-14, 1e-14, 1e-14]');