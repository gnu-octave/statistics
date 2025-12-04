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
## @deftypefnx {statistics} {@var{table} =} tabulate (@var{x})
##
## Calculate a frequency table.
##
## @code{tabulate (x)} displays a frequency table of the data in the vector
## @var{x}.  For each unique value in @var{x}, the tabulate function shows the
## number of instances and percentage of that value in @var{x}.
##
## @code{@var{table} = tabulate (@var{x})} returns the frequency table,
## @var{table}, as a numeric matrix when @var{x} is numeric and as a cell array
## otherwise.  When an output argument is requested, @code{tabulate} does not
## print the frequency table in the command window.
##
## If @var{x} is numeric, any missing values (@qcode{NaNs}) are ignored.
##
## If all the elements of @var{x} are positive integers, then the frequency
## table includes 0 counts for the integers between 1 and @qcode{max (@var{x})}
## that do not appear in @var{x}.
##
## @seealso{bar, pareto}
## @end deftypefn

function table = tabulate (x)

  ## Check input for being either numeric or a cell array
  if (! (isnumeric (x) && isvector (x)) &&
      ! (iscellstr (x) && isvector (x)) &&
      ! ischar (x))
    error (strcat ("tabulate: X must be either a numeric vector, a", ...
                   " vector cell array of strings, or a character matrix."));
  endif

  ## Remove missing values (NaNs) if numeric
  if (isnumeric (x))
    x(isnan (x)) = [];
  endif

  ## Handle positive integers separately
  if (isnumeric (x) && all (x == fix (x)) && all (x > 0))
    [count, value] = hist (x, (1:max (x)));
    posint = true;
  else
    [g, gn, gl] = grp2idx (x);
    [count, value] = hist (g, (1:length (gn)));
    posint = false;
  endif

  ## Calculate percentages
  percent = 100 * count ./ sum (count);

  ## Display results is no output argument
  if (nargout == 0)
    if (posint)
      fprintf ("   Value    Count    Percent\n");
      fprintf ("   %5d    %5d     %6.2f%%\n", value', count', percent');
    else
      valw = max (cellfun ("length", gn));
      valw = max ([5, min([50, valw])]);
      header = sprintf ("  %%%ds    %%5s    %%6s\n", valw);
      result = sprintf ("  %%%ds    %%5d     %%6.2f%%%%\n", valw);
      fprintf (header, "Value", "Count", "Percent");
      for i = 1:length (gn)
        fprintf (result, gn{i}, count(i), percent(i));
      endfor
    endif
  ## Create output table
  else
    if (posint)
      table = [value', count', percent'];
    elseif (isnumeric (x))
      table = [gl, count', percent'];
    else
      table = [gn, num2cell([count', percent'])];
    endif
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

%!error<tabulate: X must be either a numeric vector> tabulate (ones (3))
%!error<tabulate: X must be either a numeric vector> tabulate ({1, 2, 3, 4})
%!error<tabulate: X must be either a numeric vector> ...
%! tabulate ({"a", "b"; "a", "c"})
