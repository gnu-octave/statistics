## Copyright (C) 2003 Alberto Terruzzi <t-albert@libero.it>
## Copyright (C) 2023 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {statistics} {} tabulate (@var{x})
## @deftypefnx {statistics} {@var{tbl} =} tabulate (@var{x})
##
## Create a frequency table of unique values in vector @var{x}.
##
## @code{tabulate (x)} displays a frequency table of the data in the vector
## @var{x}.  The input @var{x} can be a numeric vector, a logical vector,
## a character array, a cell array of strings, a categorical array, or a
## string array.
##
## The table displays the value, the number of instances (count), and the
## percentage of that value in @var{x}.  If no output argument is requested,
## the table is displayed in the command window.
##
## @code{@var{tbl} = tabulate (@var{x})} returns the frequency table,
## @var{tbl}, as a numeric matrix when @var{x} is numeric and as a cell array
## otherwise.
##
## If @var{x} is numeric, any missing values (@qcode{NaNs}) are ignored.
## Similarly, undefined elements in categorical arrays and missing elements in
## string arrays are ignored.
##
## If all the elements of @var{x} are positive integers, then the frequency
## table includes 0 counts for the integers between 1 and @qcode{max (@var{x})}
## that do not appear in @var{x}.
##
## For categorical arrays, the frequency table includes 0 counts for any
## categories that are defined but do not appear in @var{x}.
##
## @seealso{bar, pareto}
## @end deftypefn

function tbl = tabulate (x)

  ## Check input for being numeric, string, categorical, cell array or logical
  if (! (isnumeric (x) && (isvector (x) || isempty (x))) && ! (iscellstr (x)
                       && isvector (x)) && ! ischar (x) && ! iscategorical (x)
                       && ! isa (x, "string") && ! islogical (x))
    error (strcat ("tabulate: X must be either a numeric vector, a", ...
                   " vector cell array of strings, a character matrix,", ...
                   " a categorical array, or a string array."));
  endif

  ## Ensure vector input
  if (! ischar (x))
      x = x(:);
  endif

  if (iscategorical (x))
    ## For categorical, we report ALL categories, even if count is 0
    vals = categories (x);
    nc = length (vals);

    ## Count occurrences
    xi = double (x);

    ## Filter out undefined
    valid_mask = ! isnan (xi) & (xi >= 1) & (xi <= nc);
    xi = xi(valid_mask);

    if (isempty (xi))
      counts = zeros (nc, 1);
    else
      counts = accumarray (xi, 1, [nc, 1]);
    endif

    total = sum (counts);
    percents = 100 * counts ./ total;

    ## Output format: Cell array
    out = cell (length (vals), 3);
    out(:,1) = vals;
    out(:,2) = num2cell (counts);
    out(:,3) = num2cell (percents);

  elseif (isa (x, "string"))
    ## Handle string arrays
    x(ismissing (x)) = [];

    ## Convert to cellstr and use grp2idx which is robust
    [idx, vals] = grp2idx (cellstr (x));

    if (isempty (idx))
      counts = [];
      percents = [];
    else
      counts = accumarray (idx, 1);
      total = sum (counts);
      percents = 100 * counts ./ total;
    endif

    ## Output format: Cell array
    vals_cell = vals;
    out = cell (length (vals_cell), 3);
    out(:,1) = vals_cell;
    out(:,2) = num2cell (counts);
    out(:,3) = num2cell (percents);

    if (isempty (idx))
      counts = [];
      percents = [];
    else
      counts = accumarray (idx, 1);
      total = sum (counts);
      percents = 100 * counts ./ total;
    endif

    ## Output format: Cell array
    vals_cell = vals;
    out = cell (length (vals_cell), 3);
    out(:,1) = vals_cell;
    out(:,2) = num2cell (counts);
    out(:,3) = num2cell (percents);

  elseif (islogical (x))
    ## Handle logical arrays
    [vals, ~, idx] = unique (x);
    if (isempty (x))
      counts = [];
      percents = [];
    else
      counts = accumarray (idx, 1);
      total = sum (counts);
      percents = 100 * counts ./ total;
    endif

    vals_cell = cellstr (num2str (double (vals)));

    out = cell (length (vals), 3);
    out(:, 1) = vals_cell;
    out(:, 2) = num2cell (counts);
    out(:, 3) = num2cell (percents);

  elseif (isnumeric (x))
    ## Handle numeric
    ## Remove missing values (NaNs) if numeric
    x(isnan (x)) = [];

    ## Handle positive integers separately
    if (! isempty (x) && all (x == fix (x)) && all (x > 0))
      max_val = max (x);
      vals = (1:max_val)';
      [counts, ~] = hist (x, vals);
      counts = counts(:);
    else
      [vals, ~, idx] = unique (x);
      if (isempty (x))
          counts = [];
      else
          counts = accumarray (idx, 1);
      end
    endif

    if (isempty (counts))
        percents = [];
    else
        percents = 100 * counts ./ sum (counts);
    endif

    ## Output format: Numeric Matrix
    out = [vals, counts, percents];

  else
    ## Handle char and cellstr
    if (ischar (x))
        x = cellstr (x);
    endif

    [idx, vals] = grp2idx (x);

    if (isempty (idx))
        counts = [];
    else
        counts = accumarray (idx, 1);
    endif

    if (isempty (counts))
        percents = [];
    else
        percents = 100 * counts ./ sum (counts);
    endif

    out = cell (length (vals), 3);
    out(:,1) = vals;
    out(:,2) = num2cell (counts);
    out(:,3) = num2cell (percents);
  endif

  if (nargout == 0)
    ## Use table for display if no output requested

    if (isempty (out))
       ## Handle empty case
       disp ("   Value    Count    Percent");
       return;
    endif

    if (isnumeric (out))
       ## Numeric matrix case
       Value = out(:,1);
       Count = out(:,2);
       Percent = out(:,3);
    else
       ## Cell array case
       Value = out(:,1);
       Count = cell2mat (out(:,2));
       Percent = cell2mat (out(:,3));
    endif

    t = table (Value, Count, Percent, "VariableNames", {"Value", "Count", ...
                                                        "Percent"});

    disp (t);
  else
    tbl = out;
  endif

endfunction

%!demo
%! ## Generate a frequency table for a vector of data in a cell array
%! load patients
%!
%! ## Display the first seven entries of the Gender variable
%! gender = Gender(1:7)
%!
%! ## Compute the frequency table that shows the number and
%! ## percentage of Male and Female patients
%! tabulate (Gender)

%!demo
%! ## Create a frequency table for a vector of positive integers
%! load patients
%!
%! ## Display the first seven entries of the Gender variable
%! height = Height(1:7)
%!
%! ## Create a frequency table that shows, in its second and third columns,
%! ## the number and percentage of patients with a particular height.
%! table = tabulate (Height);
%!
%! ## Display the first and last seven entries of the frequency table
%! first = table(1:7,:)
%!
%! last = table(end-6:end,:)

%!demo
%! ## Create a frequency table from a character array
%! load carsmall
%!
%! ## Tabulate the data in the Origin variable, which shows the
%! ## country of origin of each car in the data set
%! tabulate (Origin)

%!demo
%! ## Create a frequency table from a numeric vector with NaN values
%! load carsmall
%!
%! ## The carsmall dataset contains measurements of 100 cars
%! total_cars = length (MPG)
%! ## For six cars, the MPG value is missing
%! missingMPG = length (MPG(isnan (MPG)))
%!
%! ## Create a frequency table using MPG
%! tabulate (MPG)
%! table = tabulate (MPG);
%!
%! ## Only 94 cars were used
%! valid_cars = sum (table(:,2))

%!test
%! load patients
%! table = tabulate (Gender);
%! assert (table{1,1}, "Male");
%! assert (table{2,1}, "Female");
%! assert (table{1,2}, 47);
%! assert (table{2,2}, 53);
%!test
%! load patients
%! table = tabulate (Height);
%! assert (table(end-4,:), [68, 15, 15]);
%! assert (table(end-3,:), [69, 8, 8]);
%! assert (table(end-2,:), [70, 11, 11]);
%! assert (table(end-1,:), [71, 10, 10]);
%! assert (table(end,:), [72, 4, 4]);
%!test
%! ## Test numeric vector including NaNs
%! x = [1; 1; 2; 3; 1; NaN; 2];
%! tbl = tabulate (x);
%! assert (isnumeric (tbl));
%! assert (size (tbl), [3, 3]);
%! assert (tbl(:,1), [1; 2; 3]);
%! assert (tbl(:,2), [3; 2; 1]);
%! assert (tbl(:,3), [50; 33.3333; 16.6667], 3e-4);
%!test
%! ## Test positive integers with gaps
%! x = [1; 3; 3];
%! tbl = tabulate (x);
%! assert (isnumeric (tbl));
%! assert (size (tbl), [3, 3]);
%! assert (tbl(:,1), [1; 2; 3]);
%! assert (tbl(:,2), [1; 0; 2]);
%! assert (tbl(:,3), [33.3333; 0; 66.6667], 3e-4);
%!test
%! ## Test logical inputs (should return cell array with '0'/'1')
%! x = [true; false; true; true];
%! tbl = tabulate (x);
%! assert (iscell (tbl));
%! assert (size (tbl), [2, 3]);
%! assert (tbl(:,1), {'0'; '1'});
%! assert ([tbl{:,2}]', [1; 3]);
%! assert ([tbl{:,3}]', [25; 75]);
%!test
%! ## Test character array
%! x = ['a'; 'b'; 'a'];
%! tbl = tabulate (x);
%! assert (iscell (tbl));
%! assert (size (tbl), [2, 3]);
%! assert (tbl(:,1), {'a'; 'b'});
%! assert ([tbl{:,2}]', [2; 1]);
%!test
%! ## Test cell array of character vectors
%! x = {'a', 'b', 'a'};
%! tbl = tabulate (x);
%! assert (iscell (tbl));
%! assert (size (tbl), [2, 3]);
%! assert (tbl(:,1), {'a'; 'b'});
%! assert ([tbl{:,2}]', [2; 1]);
%!test
%! ## Test string array with missing values
%! x = string ({'a", 'b', 'a'});
%! x(4) = missing;
%! tbl = tabulate (x);
%! assert (iscell (tbl));
%! assert (size (tbl), [2, 3]);
%! assert (tbl(:,1), {'a'; 'b'});
%! assert ([tbl{:,2}]', [2; 1]);
%!test
%! ## Test categorical array with undefined values and vacuous levels
%! x = categorical ({'a', 'a', 'b'}, {'a', 'b', 'c'});
%! tbl = tabulate (x);
%! assert (iscell (tbl));
%! assert (size (tbl), [3, 3]);
%! assert (tbl(:,1), {'a'; 'b'; 'c'});
%! assert ([tbl{:,2}]', [2; 1; 0]);
%! assert ([tbl{:,3}]', [66.6667; 33.3333; 0], 1e-3);
%!test
%! ## Test empty input
%! tbl = tabulate ([]);
%! assert (isempty (tbl));
%!test
%! ## fisheriris (Categorical/CellStr)
%! load fisheriris;
%! tbl = tabulate (species);
%! assert (size (tbl), [3, 3]);
%! assert (tbl(:,1), {'setosa'; 'versicolor'; 'virginica'});
%! assert ([tbl{:,2}]', [50; 50; 50]);
%! assert ([tbl{:,3}]', [33.3333; 33.3333; 33.3333], 1e-4);
%!test
%! ## carsmall (Char/CellStr)
%! load carsmall;
%! tbl = tabulate (Origin);
%! origins = tbl(:,1);
%! counts = [tbl{:,2}];
%! assert (counts(strcmp (origins, 'USA')), 69);
%! assert (counts(strcmp (origins, 'Japan')), 15);
%! assert (counts(strcmp (origins, 'Germany')), 9);
%! assert (counts(strcmp (origins, 'France')), 4);
%! assert (counts(strcmp (origins, 'Sweden')), 2);
%! assert (counts(strcmp (origins, 'Italy')), 1);
%!test
%! ## patients (Logical)
%! load patients;
%! tbl = tabulate (Smoker);
%! assert (size (tbl), [2, 3]);
%! assert (tbl(:,1), {'0'; '1'});
%! assert ([tbl{:,2}]', [66; 34]);
%!test
%! ## patients (String)
%! load patients;
%! tbl = tabulate (Gender);
%! vals = tbl(:,1);
%! counts = [tbl{:,2}];
%! assert (counts(strcmp (vals, 'Male')), 47);
%! assert (counts(strcmp (vals, 'Female')), 53);

%!error<tabulate: X must be either a numeric vector> tabulate (ones (3))
%!error<tabulate: X must be either a numeric vector> tabulate ({1, 2, 3, 4})
%!error<tabulate: X must be either a numeric vector> ...
%! tabulate ({"a", "b"; "a", "c"})
