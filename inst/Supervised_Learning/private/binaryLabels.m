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
## @deftypefn {Private Function} {@var{labels} =} binaryLabels (@var{ClassNames}, @var{positive})
##
## Class labels of a binary model, in the type of the response.
##
## @var{positive} is true wherever the second class was predicted, and
## @var{labels} comes back with its shape and with the type of
## @var{ClassNames}: a cell array of character vectors, a character matrix,
## a logical array or a numeric array.
##
## @end deftypefn

function labels = binaryLabels (ClassNames, positive)

  positive = logical (positive);

  if (iscellstr (ClassNames))
    labels = cell (size (positive));
    labels(! positive) = ClassNames{1};
    labels(positive) = ClassNames{2};
  elseif (ischar (ClassNames))
    labels = char (zeros (numel (positive), columns (ClassNames)));
    labels(! positive(:),:) = repmat (ClassNames(1,:), ...
                                      sum (! positive(:)), 1);
    labels(positive(:),:) = repmat (ClassNames(2,:), sum (positive(:)), 1);
  elseif (islogical (ClassNames))
    labels = false (size (positive));
    labels(! positive) = ClassNames(1);
    labels(positive) = ClassNames(2);
  else
    labels = zeros (size (positive));
    labels(! positive) = ClassNames(1);
    labels(positive) = ClassNames(2);
  endif

endfunction
