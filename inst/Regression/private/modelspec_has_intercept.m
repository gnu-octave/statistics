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
## @deftypefn {Private Function} {@var{tf} =} modelspec_has_intercept (@var{modelspec}, @var{intercept_nv})
##
## Report whether a model specification carries an intercept.
##
## @var{modelspec} is a keyword, a numeric terms matrix, or empty, and
## @var{intercept_nv} is the @qcode{'Intercept'} name-value flag.
##
## @var{tf} is true when the fitted model will hold an intercept term.  Every
## keyword specification supplies one, so for those the flag decides alone; a
## numeric terms matrix supplies one only if it holds an all-zero row.
##
## The answer is needed @emph{before} the predictors are encoded, because
## whether a categorical predictor is given all its indicator columns or only
## @math{k - 1} of them depends on it.  It is available that early: an all-zero
## row is all-zero whether or not the matrix carries its optional trailing
## response column, so the encoded width is never needed to decide.
##
## This helper is shared by the @code{LinearModel} and
## @code{GeneralizedLinearModel} classes.
##
## @end deftypefn

function tf = modelspec_has_intercept (modelspec, intercept_nv)

  if (! intercept_nv)
    tf = false;
  elseif (isnumeric (modelspec) && ! isempty (modelspec))
    tf = any (all (double (modelspec) == 0, 2));
  else
    tf = true;
  endif

endfunction
