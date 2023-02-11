## Copyright (C) 1995-2017 Kurt Hornik
## Copyright (C) 2018 John Donoghue
## Copyright (C) 2021 Stefano Guidoni
## Copyright (C) 2022 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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

  v_length = size (varargin{1}, 1);

  ## main - begin
  v_reshape = [];                    # vector of the dimensions of t
  X = [];                            # matrix of the indexed input values
  labels = {};                       # cell array of labels

  for i = 1 : nargin
    ## If char array, convert to numerical vector
    try
      [vector, gnames] = grp2idx (varargin{i});
    catch
      error ("crosstab: x1, x2 ... xn must be vectors.");
    end_try_catch
    if ((! isvector (vector)) || (v_length != length (vector)))
      error ("crosstab: x1, x2 ... xn must be vectors of the same length.");
    endif
    X = [X, vector];
    for h = 1 : length (gnames)
      labels{h, i} = gnames{h, 1};
    endfor
    v_reshape(i) = length (unique (vector));
  endfor

  v = unique (X(:, nargin));
  t = [];

  ## core logic, this employs a recursive function "crosstab_recursive"
  ## given (x1, x2, x3, ... xn) as inputs
  ## t(i,j,k,...) = sum (x1(:) == v1(i) & x2(:) == v2(j) & ...)
  for i = 1 : length (v)
    t = [t, (crosstab_recursive (nargin - 1,...
      (X(:, nargin) == v(i) | isnan (v(i)) * isnan (X(:, nargin)))))];
  endfor

  t = reshape(t, v_reshape);         # build the nargin-dimensional matrix

  ## additional statistics
  if (length (v_reshape) > 1)
    [p, chisq] = chi2test (t);
  endif
  ## main - end


  ## function: crosstab_recursive
  ## while there are input vectors, let's do iterations over them
  function t_partial = crosstab_recursive (x_idx, t_parent)
    y = X(:, x_idx);
    w = unique (y);

    t_partial = [];
    if (x_idx == 1)
      ## we have reached the last vector,
      ## let the computation begin
      for j = 1 : length (w)
        t_partial = [t_partial, ...
          sum(t_parent & (y == w(j) | isnan (w(j)) * isnan (y)));];
      endfor
    else
      ## if there are more vectors,
      ## just add data and pass it through to the next iteration
      for j = 1 : length (w)
        t_partial = [t_partial, ...
          (crosstab_recursive (x_idx - 1, ...
            (t_parent & (y == w(j) | isnan (w(j)) * isnan (y)))))];
      endfor
    endif
  endfunction
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

