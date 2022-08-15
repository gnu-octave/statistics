## Copyright (C) 2022 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn {Function File} @var{p} = anova2 (@var{x}, @var{reps})
## @deftypefnx {Function File} @var{p} = anova2 (@var{x}, @var{reps}, @var{displayopt})
## @deftypefnx {Function File} [@var{p}, @var{atab}] = anova2 (@dots{})
## @deftypefnx {Function File} [@var{p}, @var{atab}, @var{stats}] = anova2 (@dots{})
##
## Performs two-way analysis of variance (ANOVA) with balanced designs.  For
## unbalanced designs use @qcode{anovan}.
##
## @qcode{anova2} requires two input arguments with an optional third:
##
## @itemize
## @item
## @var{x} contains the data and it must be a matrix of at least two columns and
## two rows.
## @item
## @var{reps} is the number of replicates for each combination of factor groups.
## @item
## @var{displayopt} is an optional parameter for displaying the ANOVA table,
## when it is 'on' (default) and suppressing the display when it is 'off'.
## @end itemize
##
## @qcode{anova2} returns up to three output arguments:
##
## @itemize
## @item
## @var{p} is the p-value of the null hypothesis that all group means are equal.
## @item
## @var{atab} is a cell array containing the results in a standard ANOVA table.
## @item
## @var{stats} is a structure containing statistics useful for performing
## a multiple comparison of means with the MULTCOMPARE function.
## @end itemize
##
## If anova2 is called without any output arguments, then it prints the results
## in a one-way ANOVA table to the standard output as if @var{displayopt} is
## 'on'.
##
## Examples:
##
## @example
## load popcorn;
## anova2 (popcorn, 3);
## @end example
##
##
## @example
## [p, anovatab, stats] = anova2 (popcorn, 3, "off");
## disp (p);
## @end example
##
## @end deftypefn

function [p, anovatab, stats] = anova2 (x, reps, displayopt)

  ## Check for valid number of input arguments
  narginchk (1, 3);
  ## Check for NaN values in X
  if (any (isnan( x(:))))
    error ("NaN values in input are not allowed.  Use anovan instead.");
  endif
  ## Add defaults
  if (nargin == 1)
    reps = 1;
  endif
  if (nargin < 3)
    displayopt = 'on';
  endif
  plotdata = ~(strcmp (displayopt, 'off'));

  ## Calculate group numbers
  FFGn = size (x, 1) / reps;            ## Number of groups in 1st Factor
  SFGn = size (x, 2);                   ## Number of groups in 2nd Factor

  ## Check for valid repetitions
  if (! (int16 (FFGn) == FFGn))
    error ("The number of rows in X must be a multiple of REPS.");
  else
    idx_s = 1;
    idx_e = reps;
    for i = 1:FFGn
      RIdx(i,:) = [idx_s:idx_e];
      idx_s += reps;
      idx_e += reps;
    endfor
  endif

  ## Calculate group sample sizes
  GTsz = length (x(:));                 ## Number of total samples
  FFGs = prod (size (x(RIdx(1,:),:)));  ## Number of group samples of 1st Factor
  SFGs = size (x, 1);                   ## Number of group samples of 2nd Factor

  ## Calculate group means
  GTmu = sum (x(:)) / GTsz;                 ## Grand mean of groups
  for i = 1:FFGn                            ## Group means of 1st Factor
    FFGm(i) = mean (x(RIdx(i,:),:), "all");
  endfor
  for i = 1:SFGn                            ## Group means of 2nd Factor
    SFGm(i) = mean (x(:,i));
  endfor

  ## Calculate Sum of Squares for 1st and 2nd Factors
  SSR = sum (FFGs * ((FFGm - GTmu) .^ 2));  ## Rows Sum of Squares
  SSC = sum (SFGs * ((SFGm - GTmu) .^ 2));  ## Columns Sum of Squares

  ## Calculate Total Sum of Squares
  SST = (x(:) - GTmu)' * (x(:) - GTmu);

  ## Calculate Sum of Squares Error (Within)
  SSE = 0;
  for i = 1:FFGn
    for j = 1:SFGn
      SSE += sum ((x(RIdx(i,:),j) - mean (x(RIdx(i,:),j))) .^ 2);
    endfor
  endfor

  ## Calculate degrees of freedom and Sum of Squares Interaction (if applicable)
  df_SSR = FFGn - 1;                ## 1st Factor
  df_SSC = SFGn - 1;                ## 2nd Factor
  if (reps == 1)
    df_SSE = df_SSR * df_SSC;       ## No replication, assuming additive model
  else
    df_SSE = GTsz - (FFGn * SFGn);  ## Error with replication
    df_SSI = df_SSR * df_SSC;       ## Interaction: Degrees of Freedom
    SSI = SST - SSR - SSC - SSE;    ## Interaction: Sum of Squares
  endif
  df_tot = GTsz - 1;                ## Total

  ## Calculate Mean Squares, F statistics, and p values
  if (SSE != 0)
    MSE = SSE / df_SSE;           ## Mean Square for Error (Within)
    MSR = SSR / df_SSR;           ## Mean Square for Row Factor
    F_MSR = MSR / MSE;            ## F statistic for Row Factor
    p_MSR = 1 - fcdf (F_MSR, df_SSR, df_SSE);
    MSC = SSC / df_SSC;           ## Mean Square for Column Factor
    F_MSC = MSC / MSE;            ## F statistic for Column Factor
    p_MSC = 1 - fcdf (F_MSC, df_SSC, df_SSE);
    ## With replication
    if (reps > 1)
      MSI = SSI / df_SSI;         ## Mean Square for Interaction
      F_MSI = MSI / MSE;          ## F statistic for Interaction
      p_MSI = 1 - fcdf (F_MSI, df_SSI, df_SSE);
    endif
  else      ## Special cases with no error
    if (df_SSE > 0)
        MSE = 0;                  ## Mean Square for Error (Within)
    else
        MSE = NaN;
    endif
    if (SSR == 0)                 ## No variability in Row Factor
      MSR = 0;
      F_MSR = 1;
      p_MSR = 0;
    else
      MSR = SSR / df_SSR;
      F_MSR = Inf;
      p_MSR = 1;
    endif
    if (SSC == 0)                 ## No variability in Column Factor
      MSC = 0;
      F_MSC = 0;
      p_MSC = 1;
    else
      MSC = SSC / df_SSC;
      F_MSC = Inf;
      p_MSC = 0;
    endif
    if (reps > 1 && SSI == 0)     ## Replication with no Interaction
      MSI = 0;
      F_MSI = 0;
      p_MSI = 1;
    elseif (reps > 1)             ## Replication with Interaction
      MSI = SSI / df_SSI;
      F_MSI = Inf;
      p_MSI = 0;
    endif
  endif

  ## Create p output (if requested)
  if (nargout > 0)
    if (reps > 1)
      p = [p_MSC, p_MSR, p_MSI];
    else
      p = [p_MSC, p_MSR];
    endif
  endif

  ## Create results table (if requested)
  if (nargout > 1 && reps > 1)
    anovatab = {"Source", "SS", "df", "MS", "F", "Prob>F"; ...
                "Columns", SSC, df_SSC, MSC, F_MSC, p_MSC; ...
                "Rows", SSR, df_SSR, MSR, F_MSR, p_MSR; ...
                "Interaction", SSI, df_SSI, MSI, F_MSI, p_MSI; ...
                "Error", SSE, df_SSE, MSE, "", ""; ...
                "Total", SST, df_tot, "", "", ""};
  elseif (nargout > 1 && reps == 1)
    anovatab = {"Source", "SS", "df", "MS", "F", "Prob>F"; ...
                "Columns", SSC, df_SSC, MSC, F_MSC, p_MSC; ...
                "Rows", SSR, df_SSR, MSR, F_MSR, p_MSR; ...
                "Error", SSE, df_SSE, MSE, "", ""; ...
                "Total", SST, df_tot, "", "", ""};
  endif

  ## Create stats structure (if requested) for MULTCOMPARE
  if (nargout > 2)
    stats.source = 'anova2';
    stats.sigmasq = MSE;
    stats.colmeans = SFGm(:)';
    stats.coln = SFGs;
    stats.rowmeans = FFGm(:)';
    stats.rown = FFGs;
    stats.inter = (reps > 1);
    stats.pval = p_MSI;
    stats.df = df_SSE;
  endif

  ## Print results table on screen if no output argument was requested
  if (nargout == 0 || plotdata)
    printf("                      ANOVA Table\n");
    printf("Source             SS      df        MS       F      Prob>F\n");
    printf("-----------------------------------------------------------\n");
    printf("Columns      %10.4f %5.0f %10.4f %8.2f %9.4f\n", ...
            SSC, df_SSC, MSC, F_MSC, p_MSC);
    printf("Rows         %10.4f %5.0f %10.4f %8.2f %9.4f\n", ...
            SSR, df_SSR, MSR, F_MSR, p_MSR);
    if (reps > 1)
      printf("Interaction  %10.4f %5.0f %10.4f %8.2f %9.4f\n", ...
              SSI, df_SSI, MSI, F_MSI, p_MSI);
    endif
    printf("Error        %10.4f %5.0f %10.4f\n", SSE, df_SSE, MSE);
    printf("Total        %10.4f %5.0f\n", SST, df_tot);
  endif
endfunction


%!demo
%! load popcorn;
%! anova2 (popcorn, 3);

%!demo
%! load popcorn;
%! [p, atab] = anova2(popcorn, 3, "off");
%! disp (p);

## testing against popcorn data and results from Matlab
%!test
%! popcorn = [5.5, 4.5, 3.5; 5.5, 4.5, 4.0; 6.0, 4.0, 3.0; ...
%!            6.5, 5.0, 4.0; 7.0, 5.5, 5.0; 7.0, 5.0, 4.5];
%! [p, atab] = anova2 (popcorn, 3, "off");
%! assert (p(1), 7.678957383294716e-07, 1e-14);
%! assert (p(2), 0.0001003738963050171, 1e-14);
%! assert (p(3), 0.7462153966366274, 1e-14);
%! assert (atab{2,5}, 56.700, 1e-14);
%! assert (atab{2,3}, 2, 0);
%! assert (atab{4,2}, 0.08333333333333348, 1e-14);
%! assert (atab{5,4}, 0.1388888888888889, 1e-14);
%!test
%! popcorn = [5.5, 4.5, 3.5; 5.5, 4.5, 4.0; 6.0, 4.0, 3.0; ...
%!            6.5, 5.0, 4.0; 7.0, 5.5, 5.0; 7.0, 5.0, 4.5];
%! [p, atab, stats] = anova2 (popcorn, 3, "off");
%! assert (atab{5,2}, 1.666666666666667, 1e-14);
%! assert (atab{6,2}, 22);
%! assert (stats.source, "anova2");
%! assert (stats.colmeans, [6.25, 4.75, 4]);
%! assert (stats.inter, 1, 0);
%! assert (stats.pval, 0.7462153966366274, 1e-14);
%! assert (stats.df, 12);
