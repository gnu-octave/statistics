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
## @deftypefn {Function File} vartestn (@var{x})
## @deftypefnx {Function File} vartestn (@var{x}, @var{group})
## @deftypefnx {Function File} vartestn (@dots{}, @var{name}, @var{value})
## @deftypefnx {Function File} @var{p} = vartestn (@dots{})
## @deftypefnx {Function File} [@var{p}, @var{stats}] = vartestn (@dots{})
## @deftypefnx {Function File} [@var{p}, @var{stats}] = vartestn (@dots{}, @var{name}, @var{value})
##
## Test for equal variances across multiple groups.
##
## @code{@var{h} = vartestn (@var{x})} performs Bartlett's test for equal
## variances for the columns of the matrix @var{x}.  This is a test of the null
## hypothesis that the columns of @var{x} come from normal distributions with
## the same variance, against the alternative that they come from normal
## distributions with different variances.  The result is displayed in a summary
## table of statistics as well as a box plot of the groups.
##
## @code{vartestn (@var{x}, @var{group})} requires a vector @var{x}, and a
## @var{group} argument that is a categorical variable, vector, string array, or
## cell array of strings with one row for each element of @var{x}.  Values of
## @var{x} corresponding to the same value of @var{group} are placed in the same
## group.
##
## @code{vartestn} treats NaNs as missing values, and ignores them.
##
## @code{@var{p} = vartestn (@dots{})} returns the probability of observing the
## given result, or one more extreme, by chance under the null hypothesis that
## all groups have equal variances.  Small values of @var{p} cast doubt on the
## validity of the null hypothesis.
##
## @code{[@var{p}, @var{stats}] = vartestn (@dots{})} returns a structure with
## the following fields:
##
## @multitable @columnfractions 0.05 0.2 0.75
## @item @tab "chistat" @tab -- the value of the test statistic
## @item @tab "df" @tab -- the degrees of freedom of the test
## @end multitable
##
##
## @code{[@var{p}, @var{stats}] = vartestn (@dots{}, @var{name}, @var{value})}
## specifies one or more of the following @var{name}/@var{value} pairs:
##
## @multitable @columnfractions 0.20 0.8
## @item "display" @tab "on" to display a boxplot and table, or "off" to omit
## omit these displays. Default "on".
##
## @item "testtype" @tab One of the following strings to control the type of
## test to perform:
## @end multitable
##
## @multitable @columnfractions 0.03 0.25 0.72
## @item @tab "Bartlett" @tab Bartlett's test (default).
##
## @item @tab "LeveneQuadratic" @tab Levene's test computed by performing anova
## on the squared deviations of the data values from their group means.
##
## @item @tab "LeveneAbsolute" @tab Levene's test computed by performing anova
## on the absolute deviations of the data values from their group means.
##
## @item @tab "BrownForsythe" @tab Brown-Forsythe test computed by performing
## anova on the absolute deviations of the data values from the group medians.
##
## @item @tab "OBrien" @tab O'Brien's modification of Levene's test with W=0.5.
## @end multitable
##
## The classical 'Bartlett' test is sensitive to the assumption that the
## distribution in each group is normal.  The other test types are more robust
## to non-normal distributions, especially ones prone to outliers.  For these
## tests, the STATS output structure has a field named "fstat" containing the
## test statistic, and "df1" and "df2" containing its numerator and denominator
## degrees of freedom.
##
##
## @seealso{vartest, vartest2, anova1, bartlett_test, levene_test}
## @end deftypefn

function [p, stats] = vartestn (x, group, varargin)

  ## Validate input arguments
  if (nargin < 1)
    error ("vartestn: too few input arguments.");
  endif
  if (isscalar (x))
    error ("vartestn: X must be a vector or a matrix.");
  endif
  if (nargin < 2)
    group = [];
  endif
  if (nargin > 1 && any (strcmpi (group, {"display", "testtype"})))
    varargin = [{group} varargin];
    group = [];
  endif
  if (isvector (x) && (nargin < 2 || isempty (group )))
    error ("vartestn: if X is a vector then a group vector is required.");
  endif
  ## Add defaults
  plotdata = true;
  testtype = "Bartlett";
  if (numel (varargin(:)) > 0 && mod (numel (varargin(:)), 2) == 0)
    for idx = 1:2:numel (varargin(:))
      name = varargin{idx};
      value = varargin{idx+1};
      switch (lower (name))
        case "display"
          plotdata = value;
          if (! any (strcmpi (plotdata, {"on", "off"})))
            error ("vartestn: invalid value for display.");
          endif
          if (strcmpi (plotdata, "on"))
            plotdata = true;
          else
            plotdata = false;
          endif
        case "testtype"
          testtype = value;
          if (! any (strcmpi (testtype, {"Bartlett", "LeveneAbsolute", ...
                              "LeveneQuadratic", "BrownForsythe", "OBrien"})))
            error ("vartestn: invalid value for testtype.");
          endif
        otherwise
          error ("vartestn: invalid name for optional arguments.");
      endswitch
    endfor
  elseif (numel (varargin(:)) > 0 && mod (numel (varargin(:)), 2) != 0)
    error ("vartestn: optional arguments must be in name/value pairs.");
  endif
  ## Convert group to cell array from character array, make it a column
  if (! isempty (group) && ischar (group))
    group = cellstr (group);
  endif
  if (size (group, 1) == 1)
    group = group';
  endif
  ## If x is a matrix, convert it to column vector and create a
  ## corresponging column vector for groups
  if (length (x) < prod (size (x)))
    [n, m] = size (x);
    x = x(:);
    gi = reshape (repmat ((1:m), n, 1), n*m, 1);
    if (length (group) == 0)          ## no group names are provided
      group = gi;
    elseif (size (group, 1) == m)     ## group names exist and match columns
      group = group(gi,:);
    else
      error ("vartestn: columns in X and GROUP length do not match.");
    endif
  endif
  ## Check that x and group are the same size
  if (! all (numel (x) == numel (group)))
    error ("vartestn: GROUP must be a vector with the same number of rows as x.");
  endif
  ## Identify NaN values (if any) and remove them from X along with
  ## their corresponding values from group vector
  nonan = ! isnan (x);
  x = x(nonan);
  group = group(nonan, :);
  ## Convert group to indices and separate names
  [group_id, group_names] = grp2idx (group);
  group_id = group_id(:);
  ## Compute group summary statistics
  [group_mean, group_ster, group_size] = grpstats (x, group_id, ...
                                                   {"mean", "sem", "numel"});
  ## Compute group degreed of freedom and variances
  group_DF = group_size - 1;
  groupVAR = group_size .* group_ster .^ 2;
  sum_DF = sum (group_DF);
  ## Caculate pooled variance
  if (sum_DF > 0)
     pooledVAR = sum (group_DF .* groupVAR) / sum_DF;
  else
     pooledVAR = NaN;
  end
  ## Get number of groups
  k = length (group_DF);
  ## Test for equal variance according to specified testtype
  switch (lower (testtype))
    case "bartlett"
      ## Calculate degrees of freedom
      Bdf = max(0, sum (group_DF > 0) - 1);
      ## Get valid groups
      msgroups = group_DF > 0;
      ## For valid groups
      if (Bdf > 0 && sum_DF > 0)
        B = log (pooledVAR) * sum (group_DF) - ...
            sum (group_DF(msgroups) .* log (groupVAR(msgroups)));
        C = 1 + (sum (1 ./ group_DF(msgroups)) - 1 / sum (group_DF)) / (3 * Bdf);
        F = B / C;
      else
        F = NaN;
      endif
      ## Compute p-value
      p = 1 - chi2cdf (F, Bdf);
      testname = "Bartlett's statistic            ";
      if (nargout > 1)
        stats = struct("chisqstat", F, "df", Bdf);
      endif
    case {"leveneabsolute", "levenequadratic"}
      ## Remove single-sample groups
      ssgroups = find (group_size < 2);
      msgroups = ! ismember (group_id, ssgroups);
      ## Center each group with mean
      x_center = x(msgroups) - group_mean(group_id(msgroups));
      ## Get number of valid groups (group size > 1)
      n_groups = length (group_size) - length (ssgroups);
      ## Perform one-way anova and extract results from the anova table
      if (n_groups > 1)
        if (strcmpi (testtype, "LeveneAbsolute"))
         [p, atab] = anova1 (abs (x_center), group_id(msgroups), "off");
         testname = "Levene's statistic (absolute)   ";
        else
         [p, atab] = anova1 (x_center .^ 2, group_id(msgroups), "off");
         testname = "Levene's statistic (quadratic)  ";
        endif
        ## Get F statistic and both degrees of freedom
        F = atab{2,5};
        Bdf = [atab{2,3}, atab{3,3}];
      else
        p = NaN;
        F = NaN;
        Bdf = [0, (length (x_center) - n_groups)];
      endif
      if (nargout > 1)
        stats = struct("fstat", F, "df", Bdf);
      endif
    case "brownforsythe"
      ## Remove single-sample groups
      ssgroups = find (group_size < 2);
      msgroups = ! ismember (group_id, ssgroups);
      ## Calculate group medians
      group_md = grpstats (x, group_id, "median");
      ## Center each group with median
      xcbf = x(msgroups) - group_md(group_id(msgroups));
      ## Get number of valid groups (group size > 1)
      n_groups = length(group_size) - length(ssgroups);
      ## Perform one-way anova and extract results from the anova table
      if (n_groups > 1)
        [p, atab] = anova1 (abs (xcbf), group_id(msgroups), "off");
        ## Get F statistic and both degrees of freedom
        F = atab{2,5};
        Bdf = [atab{2,3}, atab{3,3}];
      else
        p = NaN;
        F = NaN;
        Bdf = [0, (length (xcbf) - n_groups)];
      end
      testname = "Brown-Forsythe statistic        ";
      if (nargout > 1)
        stats = struct("fstat", F, "df", Bdf);
      endif
    case "obrien"
      ## Remove single-sample groups
      ssgroups = find (group_size < 2);
      msgroups = ! ismember (group_id, ssgroups);
      ## Center each group with mean
      x_center = x(msgroups) - group_mean(group_id(msgroups));
      ## Calculate OBrien Z_ij
      xcs = x_center.^2;
      W = 0.5;
      xcw = ((W + group_size(group_id(msgroups)) - 2) .* ...
              group_size(group_id(msgroups)) .* xcs - W .* ...
             (group_size(group_id(msgroups)) - 1) .* ...
              groupVAR(group_id(msgroups))) ./ ...
            ((group_size(group_id(msgroups)) - 1) .* ...
             (group_size(group_id(msgroups)) - 2));
      ## Get number of valid groups (group size > 1)
      n_groups = length(group_size) - length(ssgroups);
      ## Perform one-way anova and extract results from the anova table
      if (n_groups > 1)
        [p, atab] = anova1 (xcw, group_id(msgroups), "off");
        ## Get F statistic and both degrees of freedom
        F = atab{2,5};
        Bdf = [atab{2,3}, atab{3,3}];
      else
        p = NaN;
        F = NaN;
        Bdf = [0, length(xcw)-n_groups];
      end
      testname = "OBrien statistic                ";
      if (nargout > 1)
        stats = struct("fstat", F, "df", Bdf);
      endif
  endswitch
  ## Print Group Summary Table (unless opted out)
  if (nargout == 0 || plotdata)
    groupSTD = sqrt (groupVAR);
    printf ("\n                    Group Summary Table\n\n");
    printf ("Group                        Count        Mean       Std Dev\n");
    printf ("------------------------------------------------------------\n");
    for i = 1:k
      printf ("%-20s  %10i      %9.4f     %1.6f\n", ...
              group_names{i}, group_size(i), group_mean(i), groupSTD(i));
    endfor
    printf ("Pooled Groups         %10i      %9.4f     %1.6f\n", ...
           sum (group_size), mean (group_mean), mean (groupSTD));
    printf ("Pooled valid Groups   %10i      %9.4f     %1.6f\n\n", ...
           sum (group_size(group_id(msgroups))), ...
           mean (group_mean(group_id(msgroups))), ...
           mean (groupSTD(group_id(msgroups))));
    printf ("%s %7.5f\n", testname, F);
    if (numel (Bdf) == 1)
      printf ("Degrees of Freedom      %10i\n", Bdf);
    else
      printf ("Degrees of Freedom      %10i, %3i\n", Bdf(1), Bdf(2));
    endif
    printf ("p-value                          %1.6f\n\n", p);
  endif
  ## Plot data using BOXPLOT (unless opted out)
  if (plotdata)
    boxplot (x, group_id, "Notch", "on", "Labels", group_names);
  endif

endfunction

%!demo
%! ## Test the null hypothesis that the variances are equal across the five
%! ## columns of data in the students’ exam grades matrix, grades.
%!
%! load examgrades
%! vartestn (grades)

%!demo
%! ## Test the null hypothesis that the variances in miles per gallon (MPG) are
%! ## equal across different model years.
%!
%! load carsmall
%! vartestn (MPG, Model_Year)

%!demo
%! ## Use Levene’s test to test the null hypothesis that the variances in miles
%! ## per gallon (MPG) are equal across different model years.
%!
%! load carsmall
%! p = vartestn (MPG, Model_Year, "TestType", "LeveneAbsolute")

%!demo
%! ## Test the null hypothesis that the variances are equal across the five
%! ## columns of data in the students’ exam grades matrix, grades, using the
%! ## Brown-Forsythe test.  Suppress the display of the summary table of
%! ## statistics and the box plot.
%!
%! load examgrades
%! [p, stats] = vartestn (grades, "TestType", "BrownForsythe", "Display", "off")

## Test input validation
%!error<vartestn: too few input arguments.> vartestn ();
%!error<vartestn: X must be a vector or a matrix.> vartestn (1);
%!error<vartestn: if X is a vector then a group vector is required.> ...
%! vartestn ([1, 2, 3, 4, 5, 6, 7]);
%!error<vartestn: if X is a vector then a group vector is required.> ...
%! vartestn ([1, 2, 3, 4, 5, 6, 7], []);
%!error<vartestn: if X is a vector then a group vector is required.> ...
%! vartestn ([1, 2, 3, 4, 5, 6, 7], "TestType", "LeveneAbsolute");
%!error<vartestn: if X is a vector then a group vector is required.> ...
%! vartestn ([1, 2, 3, 4, 5, 6, 7], [], "TestType", "LeveneAbsolute");
%!error<vartestn: invalid value for display.> ...
%! vartestn ([1, 2, 3, 4, 5, 6, 7], [1, 1, 1, 2, 2, 2, 2], "Display", "some");
%!error<vartestn: invalid value for display.> ...
%! vartestn (ones (50,3), "Display", "some");
%!error<vartestn: invalid value for testtype.> ...
%! vartestn (ones (50,3), "Display", "off", "testtype", "some");
%!error<vartestn: optional arguments must be in name/value pairs.> ...
%! vartestn (ones (50,3), [], "som");
%!error<vartestn: invalid name for optional arguments.> ...
%! vartestn (ones (50,3), [], "some", "some");
%!error<vartestn: columns in X and GROUP length do not match.> ...
%! vartestn (ones (50,3), [1, 2], "Display", "off");
## Test results
%!test
%! load examgrades
%! [p, stat] = vartestn (grades, "Display", "off");
%! assert (p, 7.908647337018238e-08, 1e-14);
%! assert (stat.chisqstat, 38.7332, 1e-4);
%! assert (stat.df, 4);
%!test
%! load examgrades
%! [p, stat] = vartestn (grades, "Display", "off", "TestType", "LeveneAbsolute");
%! assert (p, 9.523239714592791e-07, 1e-14);
%! assert (stat.fstat, 8.5953, 1e-4);
%! assert (stat.df, [4, 595]);
%!test
%! load examgrades
%! [p, stat] = vartestn (grades, "Display", "off", "TestType", "LeveneQuadratic");
%! assert (p, 7.219514351897161e-07, 1e-14);
%! assert (stat.fstat, 8.7503, 1e-4);
%! assert (stat.df, [4, 595]);
%!test
%! load examgrades
%! [p, stat] = vartestn (grades, "Display", "off", "TestType", "BrownForsythe");
%! assert (p, 1.312093241723211e-06, 1e-14);
%! assert (stat.fstat, 8.4160, 1e-4);
%! assert (stat.df, [4, 595]);
%!test
%! load examgrades
%! [p, stat] = vartestn (grades, "Display", "off", "TestType", "OBrien");
%! assert (p, 8.235660885480556e-07, 1e-14);
%! assert (stat.fstat, 8.6766, 1e-4);
%! assert (stat.df, [4, 595]);
