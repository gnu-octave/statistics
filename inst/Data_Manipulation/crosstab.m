## Copyright (C) 1995-2017 Kurt Hornik
## Copyright (C) 2018 John Donoghue
## Copyright (C) 2021 Stefano Guidoni
## Copyright (C) 2022-2025 Andreas Bertsatos <abertsatos@biol.uoa.gr>
## Copyright (C) 2025 Yasin Achengli
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
## FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
## more details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn  {statistics} {@var{t} =} crosstab (@var{x1})
## @deftypefnx {statistics} {@var{t} =} crosstab (@var{x1}, @var{x2})
## @deftypefnx {statistics} {@var{t} =} crosstab (@var{x1}, @dots{}, @var{xn})
## @deftypefnx {statistics} {[@var{t}, @var{chisq}, @var{p}, @var{labels}] =} crosstab (@dots{})
## @deftypefnx {statistics} {[@dots{}] =} crosstab (@dots{}, @var{Name}, @var{Value})
##
## Create a cross-tabulation (contingency table) @var{t} from data vectors.
##
## The inputs @var{x1}, @var{x2}, ... @var{xn} must be vectors of equal length
## with a data type of numeric, logical, char array, categorical, strings,
## datetime, duration, or cell array of character vectors.  Missing values,
## such as @qcode{NaN} or @qcode{NaT}, are not counted by default.  A
## @qcode{NaN} duration is missing too, whereas MATLAB counts it as a group of
## its own.
##
## The following @qcode{Name-Value} pair is supported:
##
## @multitable @columnfractions 0.28 0.72
## @headitem @var{Name} @tab @var{Value}
##
## @item @qcode{'IncludeMissingGroups'} @tab A real double, single or logical
## scalar that is not @qcode{NaN} (default @qcode{false}).  When nonzero or
## @qcode{true}, the missing values are counted as a group of their own, placed
## after the other groups of every input that holds any.  The name is not case
## sensitive and may be shortened to @qcode{'IncludeMissing'}; any other text
## is read as data.  Only its last occurrence counts.
## @end multitable
##
## The label of the missing group is @qcode{NaN} for numeric inputs,
## @qcode{<undefined>} for categorical, @qcode{<missing>} for strings,
## @qcode{@{''@}} for cell arrays of character vectors, @qcode{NaT} for
## datetime, the text of a @qcode{NaN} element for duration, and the blank row
## itself for a char array.
##
## As additional return values @code{crosstab} returns the chi-square statistics
## @var{chisq}, its p-value @var{p} and a cell array @var{labels}, containing
## the labels of each input argument.
##
## With a single input, @var{t} is a column holding the count of each group,
## and @var{chisq} and @var{p} are @code{NaN}, there being no second variable
## to test independence against.  A char row is a single observation.
##
## A categorical input has a row or column for every one of its categories,
## an unused one included, holding zeros.  The chi-square statistic and its
## p-value are computed over the rows and columns that hold observations,
## where MATLAB returns @code{NaN} for both whenever any row or column is
## empty.  With the missing group counted, it follows every category here,
## whereas MATLAB returns a table one row longer than its labels when an unused
## category follows the last used one.
##
## @seealso{grp2idx, tabulate}
## @end deftypefn

function [t, chisq, p, labels] = crosstab (varargin)

  ## check input
  if (nargin < 1)
    print_usage ();
  endif

  ## Trailing Name-Value pairs start at the first text, from the second input
  ## on, that is the option name or its short form, case free; any other text
  ## is data.
  optnames = {"IncludeMissingGroups", "IncludeMissing"};
  isopt = @(s) ischar (s) && isrow (s) && any (strcmpi (s, optnames));
  nvars = nargin;
  incmissing = false;
  k = find (cellfun (isopt, varargin(2:end)), 1);
  if (! isempty (k))
    nvars = k;
    if (mod (nargin - nvars, 2) != 0)
      error ("crosstab: optional arguments must be Name-Value pairs.");
    endif
    for j = nvars+1:2:nargin
      name = varargin{j};
      if (! isopt (name))
        if (ischar (name) && isrow (name))
          error ("crosstab: unknown parameter name '%s'.", name);
        endif
        error ("crosstab: optional arguments must be Name-Value pairs.");
      endif
      value = varargin{j+1};
      if (! (isscalar (value) && (isfloat (value) || islogical (value))
             && isreal (value)))
        error (["crosstab: 'IncludeMissingGroups' must be a real double, ", ...
                "single or logical scalar."]);
      endif
      if (isnan (value))
        error ("crosstab: 'IncludeMissingGroups' cannot be NaN.");
      endif
      incmissing = logical (value);
    endfor
    varargin = varargin(1:nvars);
  endif

  ## main - begin
  v_length = [];                    # vector of lengths of input vectors
  reshape_format = [];              # vector of the dimensions of t
  X = [];                           # matrix of the indexed input values
  labels_data = {};                 # temporary cell array for unique labels
  coordinates = {};                 # cell array of unique elements

  for i = 1:nvars
    vector = varargin{i};
    ## Convert char, cellstr, categorical, string, datetime, or duration to
    ## indexed numeric vector
    if (ischar (vector) || iscellstr (vector) || iscategorical (vector) ||
        isstring (vector) || isdatetime (vector) || isduration (vector))
      try
        [vector, gnames] = grp2idx (vector);
        labels_data{i} = gnames;
      catch
        error ("crosstab: x1, x2 ... xn must be vectors.");
      end_try_catch
    elseif (isnumeric (vector) || islogical (vector))
      if (! isvector (vector))
        error ("crosstab: x1, x2 ... xn must be vectors.");
      endif
      vector = vector(:);
      unique_vals = unique (vector(! isnan (vector)));
      labels_data{i} = arrayfun (@num2str, unique_vals, 'UniformOutput',
                                 false);
      ## X is double: an integer vector would convert the others to its type
      vector = double (vector);
    else
      error ("crosstab: unsupported type for data vector.");
    endif
    v_length(i) = length (vector);
    if (length (unique (v_length)) != 1)
      error ("crosstab: x1, x2 ... xn must be vectors of the same length.");
    endif
    X = [X, vector];
    if (isnumeric (varargin{i}) || islogical (varargin{i}))
      reshape_format(i) = length (unique (vector(! isnan (vector))));
      coordinates(i) = unique (vector(! isnan (vector)));
    else
      ## A group has a row or column for every name it carries, an unused
      ## category included, so the table lines up with LABELS as R2024a's
      ## does.
      reshape_format(i) = numel (labels_data{i});
      coordinates(i) = (1:numel (labels_data{i}))';
    endif
    ## Missing values form a last group of their own when they are counted
    miss = isnan (vector);
    if (incmissing && any (miss))
      labels_data{i}{end+1,1} = crosstab_missing_label_ (varargin{i}, miss);
      reshape_format(i) += 1;
    endif
  endfor

  if (nargout > 3)
    max_rows = 0;
    for i = 1:nvars
      max_rows = max (max_rows, numel (labels_data{i}));
    endfor
    labels = cell (max_rows, nvars);
    for i = 1:nvars
      col_labels = labels_data{i};
      labels(1:numel (col_labels), i) = col_labels;
    endfor
  else
    labels = {};
  endif

  ## A single input gives a column, not a square matrix
  t = zeros ([reshape_format, 1]);

  ## Main logic:
  ## For each combination of x1, x2, ... xn, search in unique elements stored
  ## in coordinates for each dimension and increment the position value in t
  ## multidimensional matrix (always if there is no NaN element in the
  ## combination, or if the missing values are counted, in which case a NaN
  ## goes to the last position of its dimension).
  for idx = 1:size (X, 1)
    if (incmissing || ! any (isnan (X(idx,:))))
      location = zeros (1,size (X, 2));
      for jdx = 1:size (X,2)
        if (isnan (X(idx, jdx)))
          location(jdx) = reshape_format(jdx);
        else
          location(jdx) = find (cell2mat (coordinates(jdx)) == X(idx, jdx));
        endif
      endfor
      t(num2cell (location){:}) += 1;
    endif
  endfor

  ## The test runs over the rows and columns holding observations: an empty
  ## one has no expected count to compare with, and R2024a returns NaN for
  ## the whole table instead.
  if (nargout > 1)
    nd = ndims (t);
    keep = cell (1, nd);
    for d = 1:nd
      tot = t;
      for o = setdiff (1:nd, d)
        tot = sum (tot, o);
      endfor
      keep{d} = find (tot(:) > 0);
    endfor
    tt = t(keep{:});
    if (nvars == 1 || isscalar (tt) || isempty (tt))
      chisq = NaN;
      p = NaN;
    else
      [~, p, st] = chi2test (tt);
      chisq = st.chi2stat;
    endif
  endif
endfunction

## Label of the group that holds the missing values of the input X, whose
## missing elements are flagged by MISS.
function lbl = crosstab_missing_label_ (x, miss)

  if (iscategorical (x))
    lbl = "<undefined>";
  elseif (isstring (x))
    lbl = "<missing>";
  elseif (iscellstr (x))
    lbl = "{''}";
  elseif (ischar (x))
    lbl = x(find (miss, 1),:);
  elseif (isdatetime (x))
    lbl = "NaT";
  elseif (isduration (x))
    lbl = strtrim (char (x(find (miss, 1))));
  else
    lbl = "NaN";
  endif

endfunction

## Test input validation
%!error crosstab ()
%!error <crosstab: x1, x2 ... xn must be vectors.> crosstab (ones (2), [1 1])
%!error <crosstab: x1, x2 ... xn must be vectors.> crosstab ([1 1], ones (2))
%!error <x1, x2 .* must be vectors of the same length> crosstab ([1], [1 2])
%!error <x1, x2 .* must be vectors of the same length> crosstab ([1 2], [1])
%!error <crosstab: unsupported type for data vector.> crosstab ([1 2], {1, 2})
%!error <crosstab: x1, x2 ... xn must be vectors.> ...
%! crosstab (repmat (datetime (2026, 1, 1), 2, 2), [1 2 3 4])
%!error <crosstab: x1, x2 ... xn must be vectors.> ...
%! crosstab (hours ([1 2; 3 4]), [1 2 3 4])
%!error <crosstab: 'IncludeMissingGroups' must be a real double, single or logical scalar.> ...
%! crosstab ([1 2], [1 2], 'IncludeMissingGroups', 'yes')
%!error <crosstab: 'IncludeMissingGroups' must be a real double, single or logical scalar.> ...
%! crosstab ([1 2], [1 2], 'IncludeMissingGroups', int8 (1))
%!error <crosstab: 'IncludeMissingGroups' must be a real double, single or logical scalar.> ...
%! crosstab ([1 2], [1 2], 'IncludeMissingGroups', 1i)
%!error <crosstab: 'IncludeMissingGroups' must be a real double, single or logical scalar.> ...
%! crosstab ([1 2], [1 2], 'IncludeMissingGroups', [])
%!error <crosstab: 'IncludeMissingGroups' cannot be NaN.> ...
%! crosstab ([1 2], [1 2], 'IncludeMissingGroups', NaN)
%!error <crosstab: optional arguments must be Name-Value pairs.> ...
%! crosstab ([1 2], [1 2], 'IncludeMissingGroups')
%!error <crosstab: optional arguments must be Name-Value pairs.> ...
%! crosstab ([1 2], [1 2], 'IncludeMissingGroups', true, 5, true)
%!error <crosstab: unknown parameter name 'Foo'.> ...
%! crosstab ([1 2], [1 2], 'IncludeMissingGroups', true, 'Foo', true)
%!error <crosstab: unknown parameter name 'Include'.> ...
%! crosstab ([1 2], [1 2], 'IncludeMissing', true, 'Include', true)
%!test
%! load carbig
%! [t, chisq, p, labels] = crosstab (cyl4, when, org);
%! assert_equal (t(2,3,1), 38);
%! assert_equal (labels{3,3}, "Japan");
%!test
%! load carbig
%! [t, chisq, p, labels] = crosstab (cyl4, when, org);
%! assert_equal (t(2,3,2), 17);
%! assert_equal (labels{1,3}, "USA");
%!test
%! x = [1, 1, 2, 3, 1];
%! y = [1, 2, 5, 3, 1];
%! t = crosstab (x, y);
%! assert_equal (t, [2, 1, 0, 0; 0, 0, 0, 1; 0, 0, 1, 0]);
%!test
%! x = [1, 1, 2, 3, 1];
%! y = [1, 2, 3, 5, 1];
%! t = crosstab (x, y);
%! assert_equal (t, [2, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1]);
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
%! assert_equal (t, T);
%!test
%! x = [1, 2, NaN, 1];
%! y = [1, 2, 3, NaN];
%! t = crosstab (x, y);
%! assert_equal (t, [1, 0, 0; 0, 1, 0]);

## Test categorical input
%!test
%! x = categorical ({'A', 'B', 'A', 'C', 'B'});
%! y = [1, 2, 1, 3, 2];
%! t = crosstab (x, y);
%! assert_equal (size (t), [3, 3]);
%! assert_equal (t(1, 1), 2);  # A with 1
%! assert_equal (t(2, 2), 2);  # B with 2
%!test
%! x = categorical ({'low', 'med', 'high', 'low', 'med'});
%! y = categorical ({'X', 'Y', 'X', 'Y', 'X'});
%! t = crosstab (x, y);
%! assert_equal (size (t), [3, 2]);
%!test  # MATLAB parity: an unused category keeps its zero row
%! g = categorical ({'hi'; 'hi'; 'lo'; 'lo'; 'hi'; 'lo'}, {'mid', 'lo', 'hi'});
%! w = warning ();
%! warning ('off');
%! [t, chisq, p, labels] = crosstab (g, [1; 2; 5; 6; 3; 7] > 3);
%! warning (w);
%! assert_equal (t, [0, 0; 0, 3; 3, 0]);
%! assert_equal (labels(:,1)', {'mid', 'lo', 'hi'});
%! assert_equal (chisq, 6);
%! assert_equal (p, 0.01430587844, 1e-10);
%!test  # MATLAB parity: rows follow the category order, not the names
%! x = categorical ({'hi'; 'hi'; 'lo'; 'lo'; 'hi'; 'lo'}, {'lo', 'hi'});
%! t = crosstab (x, [1; 2; 5; 6; 3; 7] > 3);
%! assert_equal (t, [0, 3; 3, 0]);
%!test
%! ## Test categorical with numeric
%! x = categorical ([10, 20, 10, 30, 20]);
%! y = [1, 2, 1, 3, 2];
%! t = crosstab (x, y);
%! assert_equal (t, [2, 0, 0; 0, 2, 0; 0, 0, 1]);
%!test
%! smoker = [1 1 0 0 1 0 1 1 0 0 1 0]';
%! gender = [1 0 1 0 1 1 0 0 1 0 0 1]';
%! w = warning ();
%! warning ('off');
%! [t, chisq, p, labels] = crosstab (smoker, gender);
%! warning (w);
%! assert_equal (t, [2 4; 4 2]);
%! assert_equal (chisq, 1.33333333, 1e-8);
%! assert_equal (p, 0.24821308, 1e-8);
%! assert_equal (labels{1,1}, '0');
%! assert_equal (labels{1,2}, '0');
%! assert_equal (labels{2,1}, '1');
%! assert_equal (labels{2,2}, '1');
%!test
%! ## Test for categorical
%! smk_cat = categorical ([0 0 1 1 0 1 0 0 1 1 0 1]');
%! gen_cat = categorical ([0 1 0 1 0 0 1 1 0 1 1 0]');
%! w = warning ();
%! warning ('off');
%! [t, chisq, p] = crosstab (smk_cat, gen_cat);
%! warning (w);
%! assert_equal (t, [2 4; 4 2]);
%! assert_equal (chisq, 1.33333333, 1e-6);
%! assert_equal (p, 0.24821308, 1e-6);
%!test
%! x = [1 1 1 2 2 2 3 3 3]';
%! y = [1 1 1 2 2 2 3 3 3]';
%! w = warning ();
%! warning ('off');
%! [t, chisq, p, labels] = crosstab (x, y);
%! warning (w);
%! assert_equal (t, diag ([3 3 3]));
%! assert_equal (chisq, 18.00000000);
%! assert_equal (p, 0.00123410, 1e-8);
%!test
%! ## The two columns holding only rows with a missing x7 are empty, and the
%! ## test runs over the other four; R2024a returns NaN for chisq and p.
%! x7 = [1 2 3 4 NaN NaN]';
%! y7 = [10 20 30 40 50 60]';
%! w = warning ();
%! warning ('off');
%! [t, chisq, p, labels] = crosstab (x7, y7);
%! warning (w);
%! assert_equal (t, [eye(4), zeros(4, 2)]);
%! assert_equal (chisq, 12);
%! assert_equal (p, 0.2133093051, 1e-10);
%! assert_equal (labels{1,1}, '1');
%! assert_equal (labels{1,2}, '10');
%! assert_equal (labels{2,1}, '2');
%! assert_equal (labels{2,2}, '20');
%! assert_equal (labels{3,1}, '3');
%! assert_equal (labels{3,2}, '30');
%! assert_equal (labels{4,1}, '4');
%! assert_equal (labels{4,2}, '40');
%! assert_equal (isempty (labels{5,1}), true);
%! assert_equal (labels{5,2}, '50');
%! assert_equal (isempty (labels{6,1}), true);
%! assert_equal (labels{6,2}, '60');
%!test
%! a = ones (15,1);
%! b = ones (15,1);
%! w = warning ();
%! warning ('off');
%! [t, chisq, p] = crosstab (a, b);
%! warning (w);
%! assert_equal (t, 15);
%! assert_equal (isnan (chisq), true);
%! assert_equal (isnan (p), true);
%!test
%! ## all NaN → empty table + NaN stats
%! na = NaN (6,1);
%! nb = (1:6)';
%! w = warning ();
%! warning ('off');
%! [t, chisq, p] = crosstab (na, nb);
%! warning (w);
%! assert_equal (all (t(:) == 0), true);
%! assert_equal (isnan (chisq), true);
%! assert_equal (isnan (p), true);
%!test
%! ## single observation → 1×1 table, NaN statistic
%! w = warning ();
%! warning ('off');
%! [t, chisq, p, labels] = crosstab (5, 'Z');
%! warning (w);
%! assert_equal (t, 1);
%! assert_equal (isnan (chisq), true);
%! assert_equal (isnan (p), true);
%! assert_equal (labels{1,1}, '5');
%! assert_equal (labels{1,2}, 'Z');
%!test
%! xx = [1 1 1 1 2 2 2 2]';
%! yy = [10 10 10 10 20 20 20 20]';
%! w = warning ();
%! warning ('off');
%! [t, chisq, p] = crosstab (xx, yy);
%! warning (w);
%! assert_equal (t, [4 0; 0 4]);
%! assert_equal (chisq, 8.00000000);
%! assert_equal (p, 0.00467773, 1e-8);
%!test
%! set1 = repmat ((1:5)', 20, 1);
%! set2 = repmat ([1; 2], 50, 1);
%! w = warning ();
%! warning ('off');
%! [t, chisq, p] = crosstab (set1, set2);
%! warning (w);
%! assert_equal (t, 10 * ones (5,2));
%! assert_equal (chisq, 0);
%! assert_equal (p, 1);
%!test
%! ## 3-way table with NaN
%! a = [1 1 2 2 3 3 1]';
%! b = [1 2 1 2 1 2 1]';
%! c = [1 1 NaN 2 2 1 2]';
%! w = warning ();
%! warning ('off');
%! [t, chisq, p] = crosstab (a, b, c);
%! warning (w);
%! expected(:,:,1) = [1 1; 0 0; 0 1];
%! expected(:,:,2) = [1 0; 0 1; 1 0];
%! assert_equal (t, expected);
%! assert_equal (chisq, 6.00000000);
%! assert_equal (p, 0.53974935, 1e-9);
%!test
%! ## sparse cellstr
%! g = [1 5 7 12 1 5 19]';
%! h = {'A', 'B', 'A', 'C', 'D', 'E', 'A'}';
%! w = warning ();
%! warning ('off');
%! [t, chisq, p] = crosstab (g, h);
%! warning (w);
%! assert_equal (sum (t(:)), 7);
%! assert_equal (chisq, 16.33333333, 1e-8);
%! assert_equal (p, 0.42994852, 1e-8);
%!test
%! ## string array
%! str1 = ['low'; 'high'; 'med'; 'low'; 'high'];
%! str2 = ['X'; 'Y'; 'X'; 'Y'; 'X'];
%! w = warning ();
%! warning ('off');
%! [t, chisq, p] = crosstab (str1, str2);
%! warning (w);
%! assert_equal (t, [1 1; 1 1; 1 0]);
%! assert_equal (chisq, 0.83333333, 1e-8);
%! assert_equal (p, 0.659240631, 1e-9);
%!test
%! ## cellstr
%! c1 = {'A','B','A','C','B','A'}';
%! c2 = {'1','2','1','3','2','1'}';
%! w = warning ();
%! warning ('off');
%! [t, chisq, p] = crosstab (c1, c2);
%! warning (w);
%! assert_equal (t, [3 0 0; 0 2 0; 0 0 1]);
%! assert_equal (chisq, 12.00000000, 1e-14);
%! assert_equal (p, 0.01735127, 1e-8);
%!test
%! ## datetime and duration, against MATLAB R2024a
%! d = repmat ([datetime(2026, 1, 5); datetime(2025, 3, 1); ...
%!              datetime(2026, 10, 13)], 6, 1);
%! d(4) = NaT;
%! h = hours (repmat ([1; 2; 3], 6, 1));
%! w = warning ();
%! warning ('off');
%! [t, chisq, p, labels] = crosstab (d, h);
%! warning (w);
%! assert_equal (t, [0 6 0; 5 0 0; 0 0 6]);
%! assert_equal (chisq, 34, 1e-12);
%! assert_equal (p, 7.4519e-07, 1e-11);
%! assert_equal (labels, {'01-Mar-2025', '1 hr'; '05-Jan-2026', '2 hr'; ...
%!                        '13-Oct-2026', '3 hr'});
%!test
%! ## a NaN duration is missing, where R2024a makes it a group of its own
%! w = warning ();
%! warning ('off');
%! [t, chisq, p, labels] = crosstab (hours ([1; NaN; 2]), [1; 2; 2]);
%! warning (w);
%! assert_equal (t, [1 0; 0 1]);
%! assert_equal (labels, {'1 hr', '1'; '2 hr', '2'});
%!test
%! ## numeric labels without padding, as R2024a's
%! w = warning ();
%! warning ('off');
%! [~, ~, ~, labels] = crosstab ([1.5 2 10.25 1.5], [1 2 3 1]);
%! warning (w);
%! assert_equal (labels, {'1.5', '1'; '2', '2'; '10.25', '3'});
%!test
%! w = warning ();
%! warning ('off');
%! [~, ~, ~, labels] = crosstab ([1e-5 1e6 1e-5], int8 ([1 2 1]));
%! warning (w);
%! assert_equal (labels, {'1e-05', '1'; '1000000', '2'});
%!test
%! [t, chisq, p, labels] = crosstab (zeros (0, 1), zeros (0, 1));
%! assert_equal (t, []);
%! assert_equal (chisq, NaN);
%! assert_equal (p, NaN);
%! assert_equal (labels, cell (0, 2));
%!test
%! w = warning ();
%! warning ('off');
%! [t, chisq, p, labels] = crosstab (int8 ([1 2 10 1]), [-1 0.5 2 -1]);
%! warning (w);
%! assert_equal (t, [2 0 0; 0 1 0; 0 0 1]);
%! assert_equal (chisq, 8, 1e-12);
%! assert_equal (labels, {'1', '-1'; '2', '0.5'; '10', '2'});

## Test a single input
%!test
%! [t, chisq, p, labels] = crosstab ([1 2 NaN 1 2]);
%! assert_equal (t, [2; 2]);
%! assert_equal (chisq, NaN);
%! assert_equal (p, NaN);
%! assert_equal (labels, {'1'; '2'});
%!test
%! [t, chisq, p, labels] = crosstab ([1 1 2 2 2 3]');
%! assert_equal (t, [2; 3; 1]);
%! assert_equal (chisq, NaN);
%! assert_equal (p, NaN);
%! assert_equal (labels, {'1'; '2'; '3'});
%!test
%! [t, chisq, p, labels] = crosstab (1);
%! assert_equal (t, 1);
%! assert_equal (chisq, NaN);
%! assert_equal (labels, {'1'});
%!test
%! [t, ~, ~, labels] = crosstab (NaN);
%! assert_equal (t, zeros (0, 1));
%! assert_equal (labels, cell (0, 1));
%!test
%! [t, ~, ~, labels] = crosstab ('abca');
%! assert_equal (t, 1);
%! assert_equal (labels, {'abca'});
%!test
%! [t, ~, ~, labels] = crosstab (['ab'; 'cd'; 'ab']);
%! assert_equal (t, [2; 1]);
%! assert_equal (labels, {'ab'; 'cd'});
%!test
%! c = categorical ({'a'; 'b'; 'a'}, {'a', 'b', 'c'});
%! [t, ~, ~, labels] = crosstab (c);
%! assert_equal (t, [2; 1; 0]);
%! assert_equal (labels, {'a'; 'b'; 'c'});
%!test
%! t = crosstab ([1 2 NaN 1 2], 'IncludeMissingGroups', true);
%! assert_equal (t, [2; 2; 1]);

## Test IncludeMissingGroups
%!test
%! ## missing values are counted as a group of their own
%! x = [1 2 NaN 1 2];
%! y = [1 1 2 NaN 2];
%! w = warning ();
%! warning ('off');
%! [t, chisq, p, labels] = crosstab (x, y, 'IncludeMissingGroups', true);
%! warning (w);
%! assert_equal (t, [1 0 1; 1 1 0; 0 1 0]);
%! assert_equal (chisq, 3.75, 1e-12);
%! assert_equal (p, 0.440896, 1e-6);
%! assert_equal (labels, {'1', '1'; '2', '2'; 'NaN', 'NaN'});
%!test
%! ## missing values in several inputs share the combined cells
%! x = [1 NaN 2 NaN]';
%! y = [NaN NaN 1 2]';
%! z = [1 1 NaN 2]';
%! w = warning ();
%! warning ('off');
%! [t, chisq, p] = crosstab (x, y, z, 'IncludeMissingGroups', true);
%! warning (w);
%! expected = zeros (3, 3, 3);
%! expected([7 9 15 20]) = 1;
%! assert_equal (t, expected);
%! assert_equal (chisq, 26, 1e-12);
%! assert_equal (p, 0.165812, 1e-6);
%!test
%! ## the missing group of categorical, cellstr and string inputs
%! y = [1; 2; 1; 2];
%! w = warning ();
%! warning ('off');
%! [t1, ~, ~, l1] = crosstab (categorical ({'a'; 'b'; ''; 'a'}), y, ...
%!                            'IncludeMissingGroups', true);
%! [t2, ~, ~, l2] = crosstab ({'a'; 'b'; ''; 'a'}, y, ...
%!                            'IncludeMissingGroups', true);
%! [t3, ~, ~, l3] = crosstab (string ({'a'; 'b'; missing; 'a'}), y, ...
%!                            'IncludeMissingGroups', true);
%! warning (w);
%! assert_equal (t1, [1 1; 0 1; 1 0]);
%! assert_equal (t2, t1);
%! assert_equal (t3, t1);
%! assert_equal ({l1{3,1}, l2{3,1}, l3{3,1}}, ...
%!               {'<undefined>', "{''}", '<missing>'});
%!test
%! ## the missing group of datetime, duration and char inputs
%! y = [1; 2; 1; 2];
%! d = [datetime(2026, 1, 1); NaT; datetime(2026, 1, 1); datetime(2025, 1, 1)];
%! w = warning ();
%! warning ('off');
%! [t1, ~, ~, l1] = crosstab (d, y, 'IncludeMissingGroups', true);
%! [t2, ~, ~, l2] = crosstab (hours ([1; NaN; 2; 1]), y, ...
%!                            'IncludeMissingGroups', true);
%! [~, ~, ~, l3] = crosstab (days ([1; NaN; 2; 1]), y, ...
%!                           'IncludeMissingGroups', true);
%! [t4, ~, ~, l4] = crosstab (['ab'; '  '; 'ab'; 'cd'], y, ...
%!                            'IncludeMissingGroups', true);
%! warning (w);
%! assert_equal (t1, [0 1; 2 0; 0 1]);
%! assert_equal (t2, [1 1; 1 0; 0 1]);
%! assert_equal (t4, [2 0; 0 1; 0 1]);
%! assert_equal ({l1{3,1}, l2{3,1}, l3{3,1}, l4{3,1}}, ...
%!               {'NaT', 'NaN hr', 'NaN days', '  '});
%!test
%! t = crosstab ([1 2 NaN 1 2], [1 1 2 NaN 2], 'includemissing', true);
%! assert_equal (t, [1 0 1; 1 1 0; 0 1 0]);
%!test
%! t = crosstab ([1 2 NaN 1 2], [1 1 2 NaN 2], 'IncludeMissingGroups', ...
%!               single (2));
%! assert_equal (t, [1 0 1; 1 1 0; 0 1 0]);
%!test
%! t = crosstab ([1 2 NaN 1 2], [1 1 2 NaN 2], 'IncludeMissingGroups', ...
%!               true, 'IncludeMissingGroups', false);
%! assert_equal (t, [1 0; 1 1]);
