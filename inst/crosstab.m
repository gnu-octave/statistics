## Copyright (C) 1995-2017 Kurt Hornik
## Copyright (C) 2018 John Donoghue
## Copyright (C) 2021 Stefano Guidoni
## Copyright (C) 2022-2025 Andreas Bertsatos <abertsatos@biol.uoa.gr>
## Copyright (C) 2025 Yasin Achengli
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
## @deftypefn  {statistics} {@var{t} =} crosstab (@var{x1}, @var{x2})
## @deftypefnx {statistics} {@var{t} =} crosstab (@var{x1}, @dots{}, @var{xn})
## @deftypefnx {statistics} {[@var{t}, @var{chisq}, @var{p}, @var{labels}] =} crosstab (@dots{})
##
## Create a cross-tabulation (contingency table) @var{t} from data vectors.
##
## The inputs @var{x1}, @var{x2}, ... @var{xn} must be vectors of equal length
## with a data type of numeric, logical, char, or string (cell array).
##
## As additional return values @code{crosstab} returns the chi-square statistics
## @var{chisq}, its p-value @var{p} and a cell array @var{labels}, containing
## the labels of each input argument.
##
## @seealso{grp2idx, tabulate}
## @end deftypefn

function [t, chisq, p, labels] = crosstab (varargin)

  ## check input
  if (nargin < 2)
    print_usage ();
  endif

  ## main - begin
  v_length = [];                    # vector of lengths of input vectors
  reshape_format = [];              # vector of the dimensions of t
  X = [];                           # matrix of the indexed input values
  labels = {};                      # cell array of labels
  coordinates = {};                 # cell array of unique elements

  for i = 1:nargin
    vector = varargin{i};
    ## If char array, convert to numerical vector
    if (ischar (vector) || iscellstr (vector))
      try
        [vector, gnames] = grp2idx (vector);
      catch
        error ("crosstab: x1, x2 ... xn must be vectors.");
      end_try_catch
    else
      if (! isvector (vector))
        error ("crosstab: x1, x2 ... xn must be vectors.");
      endif
      vector = vector(:);
      gnames = cellstr (num2str (vector));
    endif
    v_length(i) = length (vector);
    if (length (unique (v_length)) != 1)
      error ("crosstab: x1, x2 ... xn must be vectors of the same length.");
    endif
    X = [X, vector];
    for h = 1:length (gnames)
      labels{h, i} = gnames{h};
    endfor
    reshape_format(i) = length (unique (vector(!isnan (vector))));
    coordinates(i) = unique (vector(!isnan (vector)));
  endfor

  t = zeros (reshape_format);

  ## Main logic:
  ## For each combination of x1, x2, ... xn, search in unique elements stored 
  ## in coordinates for each dimension and increment the position value in t 
  ## multidimensional matrix (always if there is no NaN element in the combination).
  for idx = 1:size (X, 1)
    if (!any (isnan (X(idx,:))))
      location = zeros (1,size (X, 2));
      for jdx = 1:size (X,2)
        location(jdx) = find (cell2mat (coordinates(jdx)) == X(idx, jdx));
      endfor
      t(num2cell (location){:}) += 1;
    endif
  endfor

  if (nargout > 1)
    if (length (reshape_format) > 1)
      [p, chisq] = chi2test (t);
    endif
  endif
endfunction

## Test input validation
%!error crosstab ()
%!error crosstab (1)
%!error <crosstab: x1, x2 ... xn must be vectors.> crosstab (ones (2), [1 1])
%!error <crosstab: x1, x2 ... xn must be vectors.> crosstab ([1 1], ones (2))
%!error <x1, x2 .* must be vectors of the same length> crosstab ([1], [1 2])
%!error <x1, x2 .* must be vectors of the same length> crosstab ([1 2], [1])
%!test
%! load carbig
%! [table, chisq, p, labels] = crosstab (cyl4, when, org);
%! assert (table(2,3,1), 38);
%! assert (labels{3,3}, "Japan");
%!test
%! load carbig
%! [table, chisq, p, labels] = crosstab (cyl4, when, org);
%! assert (table(2,3,2), 17);
%! assert (labels{1,3}, "USA");
%!test
%! x = [1, 1, 2, 3, 1];
%! y = [1, 2, 5, 3, 1];
%! t = crosstab (x, y);
%! assert (t, [2, 1, 0, 0; 0, 0, 0, 1; 0, 0, 1, 0]);
%!test
%! x = [1, 1, 2, 3, 1];
%! y = [1, 2, 3, 5, 1];
%! t = crosstab (x, y);
%! assert (t, [2, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1]);
%!test
%! x1 = [1, 3, 7, 7, 8];
%! x2 = [4, 2, 1, 1, 1];
%! x3 = [6, 2, 6, 2, NaN];
%! T1 = [0, 0, 0; 0, 1, 0; 1, 0, 0; 0, 0, 0];
%! T2 = [0, 0, 1; 0, 0, 0; 1, 0, 0; 0, 0, 0];
%! T = zeros (4, 3, 2);
%! T(:,:,1) = T1;
%! T(:,:,2) = T2;
%! t = crosstab (x1, x2, x3);
%! assert (t, T);
%!test
%! x = [1, 2, NaN, 1];
%! y = [1, 2, 3, NaN];
%! t = crosstab (x, y);
%! assert (t, [1, 0, 0; 0, 1, 0]);
