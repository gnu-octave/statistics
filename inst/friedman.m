## Copyright (C) 2022 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {statistics} {@var{p} =} friedman (@var{x})
## @deftypefnx {statistics} {@var{p} =} friedman (@var{x}, @var{reps})
## @deftypefnx {statistics} {@var{p} =} friedman (@var{x}, @var{reps}, @var{displayopt})
## @deftypefnx {statistics} {[@var{p}, @var{atab}] =} friedman (@dots{})
## @deftypefnx {statistics} {[@var{p}, @var{atab}, @var{stats}] =} friedman (@dots{})
##
## Performs the nonparametric Friedman's test to compare column effects in a
## two-way layout.  @qcode{friedman} tests the null hypothesis that the column
## effects are all the same against the alternative that they are not all the
## same.
##
## @qcode{friedman} requires one up to three input arguments:
##
## @itemize
## @item
## @var{x} contains the data and it must be a matrix of at least two columns and
## two rows.
## @item
## @var{reps} is the number of replicates for each combination of factor groups.
## If not provided, no replicates are assumed.
## @item
## @var{displayopt} is an optional parameter for displaying the Friedman's ANOVA
## table, when it is 'on' and suppressing the display when it is 'off' (default).
## @end itemize
##
## @qcode{friedman} returns up to three output arguments:
##
## @itemize
## @item
## @var{p} is the p-value of the null hypothesis that all group means are equal.
## @item
## @var{atab} is a table array containing the results of the Friedman's test in
## ANOVA table format.  The table includes columns for Source, SS, df, MS,
## Chi-sq, and Prob>Chi-sq with rows for Columns, [Interaction], Error, and Total.
## @item
## @var{stats} is a structure containing statistics useful for performing
## a multiple comparison of medians with the MULTCOMPARE function.
## @end itemize
##
## If friedman is called without any output arguments, then it prints the results
## in a Friedman's ANOVA table to the standard output.
##
## Examples:
##
## @example
## load popcorn;
## friedman (popcorn, 3);
## @end example
##
##
## @example
## [p, anovatab, stats] = friedman (popcorn, 3);
## disp (p);
## @end example
##
## @seealso{anova2, kruskalwallis, multcompare}
## @end deftypefn

function [p, table_out, stats] = friedman (x, reps, displayopt)

  ## Check for valid number of input arguments
  narginchk (1, 3);
  ## Check for NaN values in X
  if (any (isnan (x(:))))
    error ("friedman: NaN values in input are not allowed.");
  endif
  ## Add defaults
  if (nargin == 1)
    reps = 1;
  endif
  ## Check for correct size of input matrix
  [r, c] = size (x);
  if (r <= 1 || c <= 1)
    error ("friedman: bad size of input matrix.");
  endif
  if (reps > 1)
    r = r / reps;
    if (floor (r) != r)
      error ("friedman: repetitions and observations do not match.");
    endif
  endif
  ## Check for displayopt
  if (nargin < 3)
    displayopt = 'off';
  elseif ! (strcmp (displayopt, 'off') || strcmp (displayopt, 'on'))
    error ("friedman: displayopt must be either 'on' or 'off'.");
  endif
  plotdata = ! (strcmp (displayopt, 'off'));

  ## Prepare a matrix of ranks. Replicates are ranked together.
  m = x;
  sum_R = 0;
  for j = 1:r
    jrows = reps * (j - 1) + (1:reps);
    v = x(jrows,:);
    [R, tieadj] = tiedrank (v(:));
    m(jrows,:) = reshape (R, reps, c);
    sum_R = sum_R + 2 * tieadj;
  endfor

  ## Perform 2-way anova silently
  [p0, anova_table] = anova2 (m, reps, 'off');

  ## Compute Friedman test statistic and p-value
  chi_r = anova_table{2,2};
  sigmasq = c * reps * (reps * c + 1) / 12;
  if (sum_R > 0)
    sigmasq = sigmasq - sum_R / (12 * r * (reps * c - 1));
  endif
  if (chi_r > 0)
    chi_r = chi_r / sigmasq;
  endif
  p = 1 - chi2cdf (chi_r, c - 1);

  ## Create ANOVA table data for output
  if (reps > 1)
    ## When there are replicates, include interaction row
    source_list = {"Columns"; "Interaction"; "Error"; "Total"};
    ss_list = [anova_table{2,2}; anova_table{3,2}; ...
               anova_table{end - 1,2}; anova_table{end,2}];
    df_list = [anova_table{2,3}; anova_table{3,3}; ...
               anova_table{end - 1,3}; anova_table{end,3}];
    ms_list = [anova_table{2,4}; anova_table{3,4}; ...
               anova_table{end - 1,4}; 0];
    chi_sq_list = [chi_r; anova_table{3,5}; 0; 0];
    prob_list = [p; anova_table{3,6}; 0; 0];
  else
    ## When there are no replicates (reps = 1), exclude interaction row
    source_list = {"Columns"; "Error"; "Total"};
    ss_list = [anova_table{2,2}; anova_table{end - 1,2}; anova_table{end,2}];
    df_list = [anova_table{2,3}; anova_table{end - 1,3}; anova_table{end,3}];
    ms_list = [anova_table{2,4}; anova_table{end - 1,4}; 0];
    chi_sq_list = [chi_r; 0; 0];
    prob_list = [p; 0; 0];
  endif

  ## Create output table using datatypes package
  table_out = table (source_list, ss_list, df_list, ms_list, chi_sq_list, ...
                 prob_list, "VariableNames", {"Source", "SS", "df", "MS", ...
                 "Chi_sq", "Prob_Chi_sq"});

  ## Create stats structure (if requested) for MULTCOMPARE
  if (nargout > 2)
    stats.source = 'friedman';
    stats.n = r;
    stats.meanranks = mean (m);
    stats.sigma = sqrt (sigmasq);
  endif

  ## Print results table on screen if no output argument was requested
  if (nargout == 0 || plotdata)
    disp (table_out);
  endif
endfunction


%!demo
%! load popcorn;
%! friedman (popcorn, 3);

%!demo
%! load popcorn;
%! [p, atab] = friedman (popcorn, 3);
%! disp (p);

## testing against popcorn data and results from Matlab
%!test
%! popcorn = [5.5, 4.5, 3.5; 5.5, 4.5, 4.0; 6.0, 4.0, 3.0; ...
%!            6.5, 5.0, 4.0; 7.0, 5.5, 5.0; 7.0, 5.0, 4.5];
%! [p, atab] = friedman (popcorn, 3);
%! assert (p, 0.001028853354594794, 1e-14);
%! assert (atab.SS(1), 99.75, 1e-14);
%! assert (atab.df(1), 2, 0);
%! assert (atab.MS(1), 49.875, 1e-14);
%! assert (atab.Chi_sq(1), 13.75862068965517, 1e-14);
%! assert (atab.Prob_Chi_sq(1), 0.001028853354594794, 1e-14);
%!test
%! popcorn = [5.5, 4.5, 3.5; 5.5, 4.5, 4.0; 6.0, 4.0, 3.0; ...
%!            6.5, 5.0, 4.0; 7.0, 5.5, 5.0; 7.0, 5.0, 4.5];
%! [p, atab, stats] = friedman (popcorn, 3);
%! assert (atab.SS(end), 116, 0);
%! assert (atab.df(end), 17, 0);
%! assert (stats.source, 'friedman');
%! assert (stats.n, 2);
%! assert (stats.meanranks, [8, 4.75, 2.25], 0);
%! assert (stats.sigma, 2.692582403567252, 1e-14);
%!test
%! popcorn = [5.5, 4.5, 3.5; 5.5, 4.5, 4.0; 6.0, 4.0, 3.0; ...
%!            6.5, 5.0, 4.0; 7.0, 5.5, 5.0; 7.0, 5.0, 4.5];
%! s = evalc ('[p, atab] = friedman (popcorn, 3);');
%! assert (isempty (strtrim (s)));
%!test
%! popcorn = [5.5, 4.5, 3.5; 5.5, 4.5, 4.0; 6.0, 4.0, 3.0; ...
%!            6.5, 5.0, 4.0; 7.0, 5.5, 5.0; 7.0, 5.0, 4.5];
%! s = evalc ('[p, atab] = friedman (popcorn, 3, "on");');
%! assert (! isempty (strtrim (s)));
%!test
%! popcorn = [5.5, 4.5, 3.5; 5.5, 4.5, 4.0; 6.0, 4.0, 3.0; ...
%!            6.5, 5.0, 4.0; 7.0, 5.5, 5.0; 7.0, 5.0, 4.5];
%! [p, atab] = friedman (popcorn, 3);
%! assert (size (atab, 1), 4, 0);
%! assert (numel (atab.SS), size (atab, 1), 0);
%!test
%! x = [1, 2, 3; 2, 1, 3; 3, 2, 1];
%! [p, atab] = friedman (x);
%! assert (size (atab, 1), 3, 0);
%! assert (numel (atab.SS), size (atab, 1), 0);

%!error<friedman: displayopt must be either 'on' or 'off'.> ...
%! friedman ([5.5, 4.5, 3.5; 5.5, 4.5, 4.0; 6.0, 4.0, 3.0; 6.5, 5.0, 4.0; ...
%!            7.0, 5.5, 5.0; 7.0, 5.0, 4.5], 3, 'invalid_displayopt');
%!error<friedman: NaN values in input are not allowed.> ...
%! friedman ([1, 2; NaN, 4]);
%!error<friedman: repetitions and observations do not match.> ...
%! friedman ([1,2; 3,4; 5,6], 2);
