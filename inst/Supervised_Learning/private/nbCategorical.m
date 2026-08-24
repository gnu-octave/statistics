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
##

## -*- texinfo -*-
## @deftypefn {Private Function} {@var{idx} =} nbCategorical (@var{val}, @var{p}, @var{classname})
##
## Resolve the @qcode{CategoricalPredictors} argument into column indices.
##
## @var{val} is what the caller was given: empty for none, the character
## vector @qcode{'all'}, a logical vector with one element per predictor, or a
## numeric vector of column indices.  @var{p} is the number of predictors and
## @var{classname} prefixes any error raised.
##
## @var{idx} is a sorted row vector of indices, without repeats.
##
## @seealso{ClassificationNaiveBayes, fitcnb}
## @end deftypefn

function idx = nbCategorical (val, p, classname)

  if (isempty (val))
    idx = [];
    return;
  endif

  if (ischar (val) && isrow (val))
    if (! strcmpi (val, 'all'))
      error (strcat (classname, ": a character vector", ...
                     " 'CategoricalPredictors' must be 'all'."));
    endif
    idx = 1:p;
    return;
  endif

  if (islogical (val))
    if (numel (val) != p)
      error (strcat (classname, ": a logical 'CategoricalPredictors'", ...
                     " must have one element per predictor."));
    endif
    idx = find (val(:)');
    return;
  endif

  if (! (isnumeric (val) && isvector (val) && isreal (val)))
    error (strcat (classname, ": 'CategoricalPredictors' must be 'all', a", ...
                   " logical vector, or a vector of column indices."));
  endif
  idx = sort (unique (val(:)'));
  if (any (idx != fix (idx)) || any (idx < 1) || any (idx > p))
    error (strcat (classname, ": 'CategoricalPredictors' must be column", ...
                   " indices into X."));
  endif

endfunction
