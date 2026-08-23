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
## @deftypefn {Private Function} {@var{F} =} regFrame (@var{X}, @var{Y}, @var{Weights}, @var{classname})
##
## Validate the data of a regression model and resolve its observation
## weights.
##
## The counterpart of @code{classFrame} for the models with a continuous
## response: check the shapes, drop the rows that are not complete in both
## the predictors and the response, and normalize the weights to sum to one.
##
## Fields: @qcode{X} and @qcode{Y}, the retained data; @qcode{RowsUsed}, a
## logical over the rows as supplied; @qcode{Weights}, the retained weights
## before normalization; @qcode{W}, the same normalized; and @qcode{n} and
## @qcode{p}.
##
## @end deftypefn

function F = regFrame (X, Y, Weights, classname)

  if (! (isnumeric (X) && isreal (X) && ismatrix (X) && ndims (X) == 2))
    error ("%s: invalid values in X.", classname);
  endif
  if (isempty (X))
    error ("%s: X is empty.", classname);
  endif
  if (! (isnumeric (Y) && isreal (Y) && isvector (Y)))
    error ("%s: invalid values in Y.", classname);
  endif
  Y = Y(:);
  if (rows (X) != rows (Y))
    error ("%s: number of rows in X and Y must be equal.", classname);
  endif

  if (isempty (Weights))
    Weights = ones (rows (X), 1);
  else
    Weights = Weights(:);
    if (numel (Weights) != rows (X))
      error ("%s: 'Weights' must have one element per observation.", ...
             classname);
    endif
  endif

  F = struct ();
  F.RowsUsed = ! (isnan (Y) | any (isnan (X), 2));
  F.X = X(F.RowsUsed, :);
  F.Y = Y(F.RowsUsed);
  F.Weights = Weights(F.RowsUsed);
  if (isempty (F.Y))
    error ("%s: no complete observations in the data.", classname);
  endif
  [F.n, F.p] = size (F.X);

  F.W = F.Weights;
  if (sum (F.W) > 0)
    F.W = F.W / sum (F.W);
  endif

endfunction
