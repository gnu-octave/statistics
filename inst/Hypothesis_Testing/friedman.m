## Copyright (C) 2022 Andreas Bertsatos <abertsatos@biol.uoa.gr>
## Copyright (C) 2025 Swayam Shah <swayamshah66@gmail.com>
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
## @deftypefnx {statistics} {[@var{p}, @var{tbl}] =} friedman (@dots{})
## @deftypefnx {statistics} {[@var{p}, @var{tbl}, @var{stats}] =} friedman (@dots{})
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
## table, when it is 'on' (default) and suppressing the display when it is
## 'off'.  MATLAB renders the table in a figure window; this package prints it
## to the standard output, as @code{anova2} does.
## @end itemize
##
## @qcode{friedman} returns up to three output arguments:
##
## @itemize
## @item
## @var{p} is the p-value of the null hypothesis that all group means are equal.
## @item
## @var{tbl} is a cell array containing the results of the Friedman's test in
## ANOVA table format.  Its first row holds the column labels Source, SS, df,
## MS, Chi-sq and Prob>Chi-sq, followed by a row per source: Columns,
## [Interaction], Error and Total.  An entry that does not apply to a row, such
## as the chi-square statistic of the Error row, is empty.
## @item
## @var{stats} is a structure containing statistics useful for performing a
## multiple comparison of medians with the MULTCOMPARE function.
## @end itemize
##
## If friedman is called without any output arguments, then it prints the
## results in a Friedman's ANOVA table to the standard output.
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

function [p, tbl, stats] = friedman (x, reps, displayopt)

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

  ## Check for displayopt.  It is 'on' by default, as it is in MATLAB and in
  ## ANOVA1 and ANOVA2.
  disp_table = true;
  if (nargin == 3)
    if (! any (strcmp (displayopt, {'on', 'off'})))
      error ("friedman: displayopt must be either 'on' or 'off'.");
    endif
    disp_table = strcmp (displayopt, 'on');
  endif

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
    ## ANOVA2 reports Columns, Rows, Interaction, Error and Total; the
    ## Friedman table carries the interaction, row 4, and not the block
    ## effect in row 3.
    source_list = {'Columns'; 'Interaction'; 'Error'; 'Total'};
    ss_list = [anova_table{2,2}; anova_table{4,2}; ...
               anova_table{end - 1,2}; anova_table{end,2}];
    df_list = [anova_table{2,3}; anova_table{4,3}; ...
               anova_table{end - 1,3}; anova_table{end,3}];
    ms_list = {anova_table{2,4}; anova_table{4,4}; ...
               anova_table{end - 1,4}; []};
    chi_sq_list = {chi_r; []; []; []};
    prob_list = {p; []; []; []};
  else
    ## When there are no replicates (reps = 1), exclude interaction row
    source_list = {'Columns'; 'Error'; 'Total'};
    ss_list = [anova_table{2,2}; anova_table{end - 1,2}; anova_table{end,2}];
    df_list = [anova_table{2,3}; anova_table{end - 1,3}; anova_table{end,3}];
    ms_list = {anova_table{2,4}; anova_table{end - 1,4}; []};
    chi_sq_list = {chi_r; []; []};
    prob_list = {p; []; []};
  endif

  ## Create the output table as a cell array with a header row, as MATLAB
  ## does.  Entries that do not apply to a row are empty, not zero, and the
  ## column labels carry the characters a table variable name cannot hold.
  tbl = [{'Source', 'SS', 'df', 'MS', 'Chi-sq', 'Prob>Chi-sq'}; ...
         [source_list, num2cell(ss_list), num2cell(df_list), ...
          ms_list, chi_sq_list, prob_list]];

  ## Create stats structure (if requested) for MULTCOMPARE
  if (nargout > 2)
    stats.source = 'friedman';
    stats.n = r;
    stats.meanranks = mean (m);
    stats.sigma = sqrt (sigmasq);
  endif

  ## Display ANOVA table if opted or no output argument is requested.  MATLAB
  ## renders it in a figure window; this package prints it, as ANOVA2 does.
  if (nargout == 0 || disp_table)
    print_friedman_table (tbl);
  endif

endfunction

## Print the ANOVA table the way ANOVA2 prints its own: one header line, then
## a row per source, with an empty field where the statistic does not apply.
function print_friedman_table (tbl)

  printf ("\n");
  printf ("%-14s %10s %6s %10s %10s %12s\n", tbl{1,:});
  for i = 2:rows (tbl)
    printf ("%-14s", tbl{i,1});
    for j = 2:columns (tbl)
      v = tbl{i,j};
      if (isempty (v))
        printf (" %10s", "");
      elseif (j == 3)
        printf (" %6d", v);
      elseif (j == 6)
        printf (" %12.4f", v);
      else
        printf (" %10.4f", v);
      endif
    endfor
    printf ("\n");
  endfor
  printf ("\n");

endfunction


%!demo
%! load popcorn;
%! friedman (popcorn, 3);

%!demo
%! load popcorn;
%! [p, atab] = friedman (popcorn, 3, 'off');
%! disp (p);

## testing against popcorn data and results from Matlab
%!test
%! popcorn = [5.5, 4.5, 3.5; 5.5, 4.5, 4.0; 6.0, 4.0, 3.0; ...
%!            6.5, 5.0, 4.0; 7.0, 5.5, 5.0; 7.0, 5.0, 4.5];
%! [p, atab] = friedman (popcorn, 3, 'off');
%! assert_equal (p, 0.001028853354594794, 1e-14);
%! assert_equal (atab(1,:), {'Source', 'SS', 'df', 'MS', 'Chi-sq', 'Prob>Chi-sq'});
%! assert_equal (atab{2,1}, 'Columns');
%! assert_equal (atab{2,2}, 99.75, 1e-14);
%! assert_equal (atab{2,3}, 2, 0);
%! assert_equal (atab{2,4}, 49.875, 1e-14);
%! assert_equal (atab{2,5}, 13.75862068965517, 1e-14);
%! assert_equal (atab{2,6}, 0.001028853354594794, 1e-14);
%!test
%! popcorn = [5.5, 4.5, 3.5; 5.5, 4.5, 4.0; 6.0, 4.0, 3.0; ...
%!            6.5, 5.0, 4.0; 7.0, 5.5, 5.0; 7.0, 5.0, 4.5];
%! [p, atab, stats] = friedman (popcorn, 3, 'off');
%! assert_equal (atab{end,1}, 'Total');
%! assert_equal (atab{end,2}, 116, 0);
%! assert_equal (atab{end,3}, 17, 0);
%! assert_equal (stats.source, 'friedman');
%! assert_equal (stats.n, 2);
%! assert_equal (stats.meanranks, [8, 4.75, 2.25], 0);
%! assert_equal (stats.sigma, 2.692582403567252, 1e-14);
%!test
%! ## every row of the table, not only Columns and Total
%! popcorn = [5.5, 4.5, 3.5; 5.5, 4.5, 4.0; 6.0, 4.0, 3.0; ...
%!            6.5, 5.0, 4.0; 7.0, 5.5, 5.0; 7.0, 5.0, 4.5];
%! [p, atab] = friedman (popcorn, 3, 'off');
%! assert_equal (atab{3,1}, 'Interaction');
%! assert_equal (atab{3,2}, 0.083333333333258, 1e-12);
%! assert_equal (atab{3,3}, 2, 0);
%! assert_equal (atab{3,4}, 0.041666666666629, 1e-12);
%! assert_equal (atab{4,1}, 'Error');
%! assert_equal (atab{4,2}, 16.166666666666742, 1e-12);
%! assert_equal (atab{4,3}, 12, 0);
%! assert_equal (atab{4,4}, 1.347222222222229, 1e-12);
%!test
%! ## the interaction carries (c-1)(r-1) degrees of freedom
%! popcorn = [5.5, 4.5, 3.5; 5.5, 4.5, 4.0; 6.0, 4.0, 3.0; ...
%!            6.5, 5.0, 4.0; 7.0, 5.5, 5.0; 7.0, 5.0, 4.5];
%! [p, atab] = friedman (popcorn, 3, 'off');
%! c = 3;  r = 2;
%! assert_equal (atab{3,3}, (c - 1) * (r - 1), 0);
%! assert_equal (atab{2,3} + atab{3,3} + atab{4,3} < atab{5,3}, true);
%!test
%! ## without replicates the table has no interaction row
%! popcorn = [5.5, 4.5, 3.5; 5.5, 4.5, 4.0; 6.0, 4.0, 3.0; ...
%!            6.5, 5.0, 4.0; 7.0, 5.5, 5.0; 7.0, 5.0, 4.5];
%! [p, atab] = friedman (popcorn, 1, 'off');
%! assert_equal (atab(:,1), {'Source'; 'Columns'; 'Error'; 'Total'});
%! assert_equal (atab{2,2}, 12, 1e-12);
%! assert_equal (atab{3,2}, 0, 1e-12);
%! assert_equal (atab{3,3}, 10, 0);
%! assert_equal (atab{4,3}, 17, 0);
%!test
%! popcorn = [5.5, 4.5, 3.5; 5.5, 4.5, 4.0; 6.0, 4.0, 3.0; ...
%!            6.5, 5.0, 4.0; 7.0, 5.5, 5.0; 7.0, 5.0, 4.5];
%! s = evalc ('[p, atab] = friedman (popcorn, 3, "off");');
%! assert_equal (isempty (strtrim (s)), true);
%!test
%! popcorn = [5.5, 4.5, 3.5; 5.5, 4.5, 4.0; 6.0, 4.0, 3.0; ...
%!            6.5, 5.0, 4.0; 7.0, 5.5, 5.0; 7.0, 5.0, 4.5];
%! s = evalc ('[p, atab] = friedman (popcorn, 3, "on");');
%! assert_equal (! isempty (strtrim (s)), true);
%!test
%! ## the table is displayed by default, as in MATLAB and in anova1 and anova2
%! popcorn = [5.5, 4.5, 3.5; 5.5, 4.5, 4.0; 6.0, 4.0, 3.0; ...
%!            6.5, 5.0, 4.0; 7.0, 5.5, 5.0; 7.0, 5.0, 4.5];
%! s = evalc ('[p, atab] = friedman (popcorn, 3);');
%! assert_equal (! isempty (strtrim (s)), true);
%!test
%! popcorn = [5.5, 4.5, 3.5; 5.5, 4.5, 4.0; 6.0, 4.0, 3.0; ...
%!            6.5, 5.0, 4.0; 7.0, 5.5, 5.0; 7.0, 5.0, 4.5];
%! [p, atab] = friedman (popcorn, 3, 'off');
%! assert_equal (size (atab), [5, 6], 0);
%! assert_equal (iscell (atab), true);
%! assert_equal (isempty (atab{end,4}), true);
%!test
%! x = [1, 2, 3; 2, 1, 3; 3, 2, 1];
%! [p, atab] = friedman (x, 1, 'off');
%! assert_equal (size (atab), [4, 6], 0);
%! assert_equal (atab{3,1}, 'Error');
%! assert_equal (isempty (atab{2,5}), false);

%!error<friedman: displayopt must be either 'on' or 'off'.> ...
%! friedman ([5.5, 4.5, 3.5; 5.5, 4.5, 4.0; 6.0, 4.0, 3.0; 6.5, 5.0, 4.0; ...
%!            7.0, 5.5, 5.0; 7.0, 5.0, 4.5], 3, 'invalid_displayopt');
%!error<friedman: NaN values in input are not allowed.> ...
%! friedman ([1, 2; NaN, 4]);
%!error<friedman: repetitions and observations do not match.> ...
%! friedman ([1,2; 3,4; 5,6], 2);
