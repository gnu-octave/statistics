## Copyright (C) 2026 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {statistics} {@var{c} =} nancov (@var{x})
## @deftypefnx {statistics} {@var{c} =} nancov (@var{x}, @var{y})
## @deftypefnx {statistics} {@var{c} =} nancov (@dots{}, @var{normalization})
## @deftypefnx {statistics} {@var{c} =} nancov (@dots{}, @var{method})
##
## Compute the covariance matrix while ignoring NaN values.
##
## @code{@var{c} = nancov (@var{x})} returns the covariance matrix of the
## columns of @var{x}, treating each row as an observation, after removing
## @qcode{NaN} values.  If @var{x} is a vector, the scalar variance of its
## non-@qcode{NaN} elements is returned.
##
## @code{@var{c} = nancov (@var{x}, @var{y})}, where @var{x} and @var{y} are of
## equal length, is equivalent to @code{nancov ([@var{x}(:), @var{y}(:)])} and
## returns the 2-by-2 covariance matrix.
##
## @code{@var{c} = nancov (@dots{}, @var{normalization})} specifies the
## normalization.  When @var{normalization} is 0 (default), the covariance is
## normalized by @math{N-1}, where @math{N} is the number of observations used.
## When it is 1, it is normalized by @math{N}.
##
## @code{@var{c} = nancov (@dots{}, @var{method})} selects how @qcode{NaN}
## values are handled.  With @qcode{"complete"} (the default), any row of the
## data that contains a @qcode{NaN} value is removed before the covariance is
## computed.  With @qcode{"pairwise"}, each element @code{@var{c}(i,j)} is
## computed using all rows in which both column @var{i} and column @var{j} are
## non-@qcode{NaN}; the resulting matrix may fail to be positive semidefinite.
##
## @seealso{cov, nanvar, nanstd, nanmean}
## @end deftypefn

function c = nancov (varargin)

  if (nargin < 1)
    print_usage ();
  endif
  x = varargin{1};
  if (! (isnumeric (x) || islogical (x)) || ! isreal (x))
    error ("nancov: X must be a real numeric matrix or vector.");
  endif
  args = varargin(2:end);

  ## Separate a trailing method string from the numeric arguments
  method = 'complete';
  strmask = cellfun (@ischar, args);
  if (any (strmask))
    sopt = args(strmask);
    if (numel (sopt) > 1)
      error ("nancov: only one METHOD option may be specified.");
    endif
    method = lower (sopt{1});
    if (! any (strcmp (method, {'complete', 'pairwise'})))
      error ("nancov: METHOD must be 'complete' or 'pairwise'.");
    endif
    args = args(! strmask);
  endif

  ## Remaining numeric arguments: an optional Y and/or a normalization flag
  y = [];
  nrm = 0;
  if (numel (args) == 1)
    a = args{1};
    if (isscalar (a) && (a == 0 || a == 1))
      nrm = a;
    else
      y = a;
    endif
  elseif (numel (args) == 2)
    y = args{1};
    nrm = args{2};
    if (! (isscalar (nrm) && (nrm == 0 || nrm == 1)))
      error ("nancov: normalization flag must be 0 or 1.");
    endif
  elseif (numel (args) > 2)
    error ("nancov: too many input arguments.");
  endif

  ## Assemble the data matrix (observations in rows, variables in columns)
  if (! isempty (y))
    if (! (isnumeric (y) || islogical (y)) || ! isreal (y))
      error ("nancov: Y must be a real numeric matrix or vector.");
    endif
    if (numel (x) != numel (y))
      error ("nancov: X and Y must have the same number of elements.");
    endif
    X = [x(:), y(:)];
  elseif (isvector (x))
    X = x(:);
  else
    X = x;
  endif

  p = columns (X);

  if (strcmp (method, 'complete'))
    good = all (! isnan (X), 2);
    Xc = X(good, :);
    n = rows (Xc);
    if (n == 0)
      c = NaN (p, p);
      return;
    endif
    Xd = Xc - mean (Xc, 1);
    if (nrm == 1 || n == 1)
      d = n;
    else
      d = n - 1;
    endif
    c = (Xd' * Xd) / d;
  else
    c = zeros (p, p);
    for i = 1:p
      for j = i:p
        rc = ! isnan (X(:,i)) & ! isnan (X(:,j));
        xi = X(rc, i);
        xj = X(rc, j);
        n = numel (xi);
        if (n == 0)
          cij = NaN;
        else
          xi -= mean (xi);
          xj -= mean (xj);
          if (nrm == 1 || n == 1)
            d = n;
          else
            d = n - 1;
          endif
          cij = sum (xi .* xj) / d;
        endif
        c(i,j) = cij;
        c(j,i) = cij;
      endfor
    endfor
  endif

endfunction

%!demo
%! ## Covariance matrix of a data set with missing values (complete-case).
%!
%! x = [1 2 3; 4 5 NaN; 7 NaN 9; 10 11 12; NaN 14 15]
%! c = nancov (x)

%!demo
%! ## The same data set using pairwise deletion of missing values.
%!
%! x = [1 2 3; 4 5 NaN; 7 NaN 9; 10 11 12; NaN 14 15]
%! c = nancov (x, 'pairwise')

## Test output
%!assert (nancov ([1 2 3; 4 5 NaN; 7 NaN 9; 10 11 12; NaN 14 15]), ...
%!        40.5 * ones (3))
%!assert (nancov ([1 2 3; 4 5 NaN; 7 NaN 9; 10 11 12; NaN 14 15], 1), ...
%!        20.25 * ones (3))
%!assert (nancov ([1 2 3; 4 5 NaN; 7 NaN 9; 10 11 12; NaN 14 15], 'pairwise'), ...
%!        [15, 21, 21; 21, 30, 39; 21, 39, 26.25])
%!assert (nancov ([1 2 3; 4 5 NaN; 7 NaN 9; 10 11 12; NaN 14 15], 1, ...
%!        'pairwise'), [11.25, 14, 14; 14, 22.5, 26; 14, 26, 19.6875])
%!assert (nancov ([1 2 3 NaN 5]', [2 NaN 6 8 10]'), [4, 8; 8, 16])
%!assert (nancov ([1 2 3 4 5]'), 2.5)
%!assert (nancov (5), 0)
%!assert (nancov (NaN (3, 2)), NaN (2, 2))

## Test input validation
%!error <Invalid call to nancov> nancov ()
%!error <nancov: X must be a real numeric matrix or vector.> nancov ({1})
%!error <nancov: METHOD must be 'complete' or 'pairwise'.> ...
%! nancov ([1 2; 3 4], 'bogus')
%!error <nancov: X and Y must have the same number of elements.> ...
%! nancov ([1 2 3], [1 2])
%!error <nancov: too many input arguments.> nancov ([1 2], [3 4], 0, 1)
