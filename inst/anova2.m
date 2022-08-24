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
## @deftypefnx {Function File} @var{p} = anova2 (@var{x}, @var{reps}, @var{displayopt}, @var{model})
## @deftypefnx {Function File} [@var{p}, @var{atab}] = anova2 (@dots{})
## @deftypefnx {Function File} [@var{p}, @var{atab}, @var{stats}] = anova2 (@dots{})
##
## Performs two-way factorial (crossed) or a nested analysis of variance (ANOVA)  
## for balanced designs. For unbalanced factorial designs use @qcode{anovan}.
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
## @item
## @var{model} is an optional parameter to specify the model type as either:
## 
## @itemize
## @item
## "interaction" or "full" (default): compute both main effects and their
## interaction
##
## @item
## "linear" (default) : compute both main effects without an interaction (e.g.
## balanced randomized block design).
##
## @item
## "nested" : treat the row factor as nested within columns. Note that the row
## factor is considered a random factor in the calculation of the statistics.
##
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

function [p, anovatab, stats] = anova2 (x, reps, displayopt, model)

  ## Check for valid number of input arguments
  narginchk (1, 4);
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
  if (nargin < 4)
    model = "interaction";
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
  if (reps > 1)
    SSE = 0;
    for i = 1:FFGn
      for j = 1:SFGn
        SSE += sum ((x(RIdx(i,:),j) - mean (x(RIdx(i,:),j))) .^ 2);
      endfor
    endfor
  else
    SSE = SST - SSC - SSR;
  endif

  ## Calculate degrees of freedom and Sum of Squares Interaction (if applicable)
  df_SSR = FFGn - 1;                ## 1st Factor
  df_SSC = SFGn - 1;                ## 2nd Factor
  if (reps > 1)
    df_SSE = GTsz - (FFGn * SFGn);  ## Error with replication
    df_SSI = df_SSR * df_SSC;       ## Interaction: Degrees of Freedom
    SSI = SST - SSR - SSC - SSE;    ## Interaction: Sum of Squares
  else
    df_SSE = df_SSR * df_SSC;       ## No replication, assuming additive model
    df_SSI = 0;
    SSI = 0;
  endif
  df_tot = GTsz - 1;                ## Total

  ## Model-specific calculations of sums-of-squares, mean squares and degrees of 
  ## freedom. The calculations are based on equalities for the partitioning of
  ## variance in fully balanced designs.
  switch (lower (model))
    case {"interaction","full"}
      ## TWO-WAY ANOVA WITH INTERACTION (full factorial model)
      ## Sums--of-squares are already partitioned into main effects and
      ## interaction. Just calculate mean-squares and degrees of fredom
      model = "interaction";
      MSE = SSE / df_SSE;           ## Mean Square for Error (Within)
      MSR = SSR / df_SSR;           ## Mean Square for Row Factor
      MS_DENOM = MSE;
      df_DENOM = df_SSE;
    case "linear"
      ## TWO-WAY ANOVA WITHOUT INTERACTION (additive, linear model)
      ## Pool Error and Interaction term
      model = "linear";
      SSE += SSI;
      df_SSE += df_SSI;
      SSI = 0;
      df_SSI = 0;
      reps = 1;                   ## Set reps to 1 to avoid printing interaction
      MSE = SSE / df_SSE;         ## Mean Square for Error (Within)
      MSR = SSR / df_SSR;         ## Mean Square for Row Factor
      MS_DENOM = MSE;
      df_DENOM = df_SSE;
    case "nested"
      ## NESTED ANOVA
      ## Row Factor is nested within Column Factor. Treat Row factor as random.
      ## Pool Row Factor and Interaction term
      model = "nested";
      SSR += SSI;
      df_SSR += df_SSI;
      SSI = 0;
      df_SSI = 0;
      reps = 1;                   ## Set reps to 1 to avoid printing interaction
      MSE = SSE / df_SSE;         ## Mean Square for Error (Within)
      MSR = SSR / df_SSR;         ## Mean Square for Row Factor
      MS_DENOM = MSR;             ## Row factor is random so MSR is denominator
      df_DENOM = df_SSR;          ## Row factor is random so df_SSR is denominator
    otherwise
      error ("model type not recognised");
  endswitch

  ## Calculate F statistics and p values
  MSE = SSE / df_SSE;           ## Mean Square for Error (Within)
  MSR = SSR / df_SSR;           ## Mean Square for Row Factor
  F_MSR = MSR / MSE;            ## F statistic for Row Factor
  p_MSR = 1 - fcdf (F_MSR, df_SSR, df_SSE);
  MSC = SSC / df_SSC;           ## Mean Square for Column Factor
  F_MSC = MSC / MS_DENOM;          ## F statistic for Column Factor
  p_MSC = 1 - fcdf (F_MSC, df_SSC, df_DENOM);

  ## With replication
  if (reps > 1)
    MSI = SSI / df_SSI;         ## Mean Square for Interaction
    F_MSI = MSI / MSE;          ## F statistic for Interaction
    p_MSI = 1 - fcdf (F_MSI, df_SSI, df_SSE);
  else
    MSI = 0;
    F_MSI = 0;
    p_MSI = NaN;                
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
    stats.sigmasq = MS_DENOM; ## MS used to calculate F relating to stats.pval
    stats.colmeans = SFGm(:)';
    stats.coln = SFGs;
    stats.rowmeans = FFGm(:)';
    stats.rown = FFGs;
    stats.inter = (reps > 1);
    if stats.inter
      stats.pval = p_MSI;     ## Interaction p-value if stats.inter is true
    else
      stats.pval = p_MSC;     ## Column Factor p-value if stats.inter is false
    end
    stats.df = df_DENOM;      ## Degrees of freedom used to calculate stats.pval
    stats.model = model;
  endif

  ## Print results table on screen if no output argument was requested
  if (nargout == 0 || plotdata)
    printf("\n                      ANOVA Table\n");
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
    printf("Total        %10.4f %5.0f\n\n", SST, df_tot);
  endif
endfunction


%!demo
%!
%! popcorn = [5.5, 4.5, 3.5; 5.5, 4.5, 4.0; 6.0, 4.0, 3.0; ...
%!            6.5, 5.0, 4.0; 7.0, 5.5, 5.0; 7.0, 5.0, 4.5];
%! anova2 (popcorn, 3);

%!demo
%!
%! popcorn = [5.5, 4.5, 3.5; 5.5, 4.5, 4.0; 6.0, 4.0, 3.0; ...
%!            6.5, 5.0, 4.0; 7.0, 5.5, 5.0; 7.0, 5.0, 4.5];
%! [p, atab] = anova2(popcorn, 3, "off");
%! disp (p);

%!demo
%!
%! data = [4.5924 7.3809 21.322; -0.5488 9.2085 25.0426; ...
%!         6.1605 13.1147 22.66; 2.3374 15.2654 24.1283; ...
%!         5.1873 12.4188 16.5927; 3.3579 14.3951 10.2129; ...
%!         6.3092 8.5986 9.8934; 3.2831 3.4945 10.0203];
%!
%! [p, atab, stats] = anova2 (data,4,"on","nested");


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
%!test
%! data = [4.5924 7.3809 21.322; -0.5488 9.2085 25.0426; ...
%!         6.1605 13.1147 22.66; 2.3374 15.2654 24.1283; ...
%!         5.1873 12.4188 16.5927; 3.3579 14.3951 10.2129; ...
%!         6.3092 8.5986 9.8934; 3.2831 3.4945 10.0203];
%! [p, atab, stats] = anova2 (data,4,"off","nested");
%! assert (atab{2,2}, 745.360306290833, 1e-10);
%! assert (atab{3,2}, 278.01854140125, 1e-10);
%! assert (atab{4,2}, 180.180377467501, 1e-10);
%! assert (atab{5,2}, 1203.55922515958, 1e-10);
%! assert (atab{2,4}, 372.680153145417, 1e-10);
%! assert (atab{3,4}, 92.67284713375, 1e-10);
%! assert (atab{4,4}, 10.0100209704167, 1e-10);
%! assert (atab{2,5}, 4.02146005730833, 1e-10);
%! assert (atab{3,5}, 9.25800729165627, 1e-10);
%! assert (atab{2,6}, 0.141597630656771, 1e-10);
%! assert (atab{3,6}, 0.000636643812875719, 1e-10);
