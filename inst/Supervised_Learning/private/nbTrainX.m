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
## @deftypefn {Private Function} {[@var{X}, @var{Y}] =} nbTrainX (@var{obj})
##
## The data a naive Bayes model was fitted on, as it was fitted.
##
## @qcode{X} and @qcode{Y} are stored as they were supplied, so the rows that
## were dropped for holding a missing value, or for belonging to a class the
## model was not fitted on, are removed here.  Every @code{resub} method goes
## through this, so none of them can disagree with the fit about which
## observations the model actually saw.
##
## @seealso{ClassificationNaiveBayes}
## @end deftypefn

function [X, Y] = nbTrainX (obj)

  X = obj.X;
  Y = obj.Y;
  if (! isempty (obj.RowsUsed))
    X = X(obj.RowsUsed, :);
    Y = Y(obj.RowsUsed, :);
  endif

endfunction
