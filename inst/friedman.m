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
## table, when it is 'on' (default) and suppressing the display when it is 'off'.
## @end itemize
##
## @qcode{friedman} returns up to three output arguments:
##
## @itemize
## @item
## @var{p} is the p-value of the null hypothesis that all group means are equal.
## @item
## @var{atab} is a cell array containing the results in a Friedman's ANOVA table.
## @item
## @var{stats} is a structure containing statistics useful for performing
## a multiple comparison of medians with the MULTCOMPARE function.
## @end itemize
##
## If friedman is called without any output arguments, then it prints the results
## in a one-way ANOVA table to the standard output as if @var{displayopt} is
## 'on'.
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
## [p, anovatab, stats] = friedman (popcorn, 3, "off");
## disp (p);
## @end example
##
## @seealso{anova2, kruskalwallis, multcompare}
## @end deftypefn

function [p, table, stats] = friedman (x, reps, displayopt)

  ## Check for valid number of input arguments
  narginchk (1, 3);
  ## Check for NaN values in X
  if (any (isnan( x(:))))
    error ("NaN values in input are not allowed.");
  endif
  ## Add defaults
  if (nargin == 1)
    reps = 1;
  endif
  ## Check for correct size of input matrix
  [r, c] = size (x);
  if (r <= 1 || c <= 1)
    error ("Bad size of input matrix.");
  endif
  if (reps > 1)
    r = r / reps;
    if (floor (r) != r)
      error ("Repetitions and observations do not match.");
    endif
  endif
  ## Check for displayopt
  if (nargin < 3)
    displayopt = 'on';
  elseif ! (strcmp (displayopt, "off") || strcmp (displayopt, "on"))
    error ("displayopt must be either 'on' (default) or 'off'.");
  endif
  plotdata = ~(strcmp (displayopt, "off"));

  ## Prepare a matrix of ranks. Replicates are ranked together.
  m = x;
  sum_R = 0;
  for j = 1:r
    jrows = reps * (j-1) + (1:reps);
    v = x(jrows,:);
    [R, tieadj] = tiedrank (v(:));
    m(jrows,:) = reshape (R, reps, c);
    sum_R = sum_R + 2 * tieadj;
  endfor

  ## Perform 2-way anova silently
  [p0, table] = anova2 (m, reps, 'off');

  ## Compute Friedman test statistic and p-value
  chi_r = table{2,2};
  sigmasq = c * reps * (reps * c + 1) / 12;
  if (sum_R > 0)
     sigmasq = sigmasq - sum_R / (12 * r * (reps * c - 1));
  endif
  if (chi_r > 0)
     chi_r = chi_r / sigmasq;
  endif
  p = 1 - chi2cdf (chi_r, c - 1);

  ## Remove row info from ANOVA2 table
  table(3,:) = [];
  ## Remove interaction chi-sq and p-value, if there are repetitive measurements
  if (reps > 1)
    table{3,5} = [];
    table{3,6} = [];
  endif
  ## Fix test statistic names
  table{1,5} = "Chi-sq";
  table{1,6} = "Prob>Chi-sq\n";
  ## Fix test statistic values
  table{2,5} = chi_r;
  table{2,6} = p;

  ## Create stats structure (if requested) for MULTCOMPARE
  if (nargout > 2)
    stats.source = 'friedman';
    stats.n = r;
    stats.meanranks = mean (m);
    stats.sigma = sqrt (sigmasq);
  endif

  ## Print results table on screen if no output argument was requested
  if (nargout == 0 || plotdata)
    printf("              Friedman's ANOVA Table\n");
    printf("Source            SS      df        MS    Chi-sq    Prob>Chi-sq\n");
    printf("---------------------------------------------------------------\n");
    printf("Columns      %10.4f %5.0f %10.4f %8.2f %9.4f\n", ...
            table{2,2}, table{2,3}, table{2,4}, table{2,5}, table{2,6});
    if reps > 1
      printf("Interaction  %10.4f %5.0f %10.4f %8.2f %9.4f\n", ...
              table{3,2}, table{3,3}, table{3,4}, table{3,5}, table{3,6});
    endif
    printf("Error        %10.4f %5.0f %10.4f\n", ...
            table{end-1,2}, table{end-1,3}, table{end-1,4});
    printf("Total        %10.4f %5.0f\n", table{end,2}, table{end,3});
  endif
endfunction


%!demo
%! load popcorn;
%! friedman (popcorn, 3);

%!demo
%! load popcorn;
%! [p, atab] = friedman (popcorn, 3, "off");
%! disp (p);

## testing against popcorn data and results from Matlab
%!test
%! popcorn = [5.5, 4.5, 3.5; 5.5, 4.5, 4.0; 6.0, 4.0, 3.0; ...
%!            6.5, 5.0, 4.0; 7.0, 5.5, 5.0; 7.0, 5.0, 4.5];
%! [p, atab] = friedman (popcorn, 3, "off");
%! assert (p, 0.001028853354594794, 1e-14);
%! assert (atab{2,2}, 99.75, 1e-14);
%! assert (atab{2,3}, 2, 0);
%! assert (atab{2,4}, 49.875, 1e-14);
%! assert (atab{2,5}, 13.75862068965517, 1e-14);
%! assert (atab{2,6}, 0.001028853354594794, 1e-14);
%! assert (atab{3,2}, 0.08333333333333215, 1e-14);
%! assert (atab{3,4}, 0.04166666666666607, 1e-14);
%! assert (atab{4,3}, 12, 0);
%!test
%! popcorn = [5.5, 4.5, 3.5; 5.5, 4.5, 4.0; 6.0, 4.0, 3.0; ...
%!            6.5, 5.0, 4.0; 7.0, 5.5, 5.0; 7.0, 5.0, 4.5];
%! [p, atab, stats] = friedman (popcorn, 3, "off");
%! assert (atab{5,2}, 116, 0);
%! assert (atab{5,3}, 17, 0);
%! assert (stats.source, "friedman");
%! assert (stats.n, 2);
%! assert (stats.meanranks, [8, 4.75, 2.25], 0);
%! assert (stats.sigma, 2.692582403567252, 1e-14);
