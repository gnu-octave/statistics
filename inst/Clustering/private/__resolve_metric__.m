## Copyright (C) 2026 Andreas Bertsatos <abertsatos@biol.uoa.gr>
##
## This file is part of the statistics package for GNU Octave.
##
## This program is free software: you can redistribute it and/or
## modify it under the terms of the GNU General Public License as
## published by the Free Software Foundation, either version 3 of the
## License, or (at your option) any later version.
##
## This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; see the file COPYING.  If not, see
## <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn {Private Function} {[@var{dist}, @var{param}] =} __resolve_metric__ (@var{fname}, @var{X}, @var{dist}, @var{P}, @var{C}, @var{S}, @var{metrics})
##
## Validate a distance metric and its parameter for a per-call override.
##
## @var{fname} names the caller for the error messages, @var{X} is the training
## data a default parameter is derived from, and @var{dist} the metric to use.
## @var{P}, @var{C} and @var{S} carry the @qcode{'P'}, @qcode{'Cov'} and
## @qcode{'Scale'} values, each empty when not supplied, and @var{metrics} lists
## the metrics the caller accepts.
##
## The returned @var{param} is the distance parameter belonging to @var{dist}:
## the supplied value where there is one, otherwise the same default the
## constructor would derive from @var{X}.  Nothing is written back to the
## searcher, so an override applies only to the call that asked for it.
##
## @end deftypefn

function [dist, param] = __resolve_metric__ (fname, X, dist, P, C, S, metrics)

  if (! (ischar (dist) && any (strcmpi (dist, metrics))))
    error ("%s: unsupported distance metric '%s'.", fname, dist);
  endif
  dist = lower (dist);

  if (strcmpi (dist, 'minkowski'))
    if (isempty (P))
      param = 2;
    else
      if (! (isscalar (P) && isnumeric (P) && isfinite (P) && P > 0))
        error ("%s: 'P' must be a positive finite scalar.", fname);
      endif
      param = P;
    endif
  elseif (strcmpi (dist, 'seuclidean'))
    if (isempty (S))
      param = std (X, [], 1);
    else
      if (! (isvector (S) && isnumeric (S) && all (S >= 0)
             && numel (S) == columns (X)))
        error (strcat ("%s: 'Scale' must be a nonnegative vector", ...
                       " matching the columns of X."), fname);
      endif
      param = S(:)';
    endif
  elseif (strcmpi (dist, 'mahalanobis'))
    if (isempty (C))
      param = cov (X);
    else
      if (! (isnumeric (C) && ismatrix (C) && rows (C) == columns (C)
             && rows (C) == columns (X)))
        error (strcat ("%s: 'Cov' must be a square matrix matching the", ...
                       " columns of X."), fname);
      endif
      if (! issymmetric (C))
        error ("%s: 'Cov' must be symmetric for mahalanobis.", fname);
      endif
      [~, p] = chol (C);
      if (p != 0)
        error ("%s: 'Cov' must be positive definite for mahalanobis.", fname);
      endif
      param = C;
    endif
  else
    param = [];
  endif

  ## A parameter that belongs to a different metric is a mistake, not a value
  ## to be quietly dropped
  if (! isempty (P) && ! strcmpi (dist, 'minkowski'))
    error ("%s: 'P' applies only to the minkowski metric.", fname);
  endif
  if (! isempty (S) && ! strcmpi (dist, 'seuclidean'))
    error ("%s: 'Scale' applies only to the seuclidean metric.", fname);
  endif
  if (! isempty (C) && ! strcmpi (dist, 'mahalanobis'))
    error ("%s: 'Cov' applies only to the mahalanobis metric.", fname);
  endif

endfunction
