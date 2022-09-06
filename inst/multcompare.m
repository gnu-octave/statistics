## Copyright (C) 2022 Andrew Penn <A.C.Penn@sussex.ac.uk>
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
## @deftypefn {Function File} @var{C} = multcompare (@var{STATS})
## @deftypefnx {Function File} @var{C} = multcompare (@var{STATS}, "name", @var{value})
## @deftypefnx {Function File} [@var{C}, @var{M}] = multcompare (...)
## @deftypefnx {Function File} [@var{C}, @var{M}, @var{H}] = multcompare (...)
## @deftypefnx {Function File} [@var{C}, @var{M}, @var{H}, @var{GNAMES}] = multcompare (...)
##
## Perform posthoc multiple comparison tests after ANOVA or ANCOVA tests.
##
## @code{@var{C} = multcompare (@var{STATS})} performs a multiple comparison
## using a @var{STATS} structure that is obtained as output from any of
## the following functions:  anovan.
## The return value @var{C} is a matrix with one row per comparison and six
## columns. Columns 1-2 are the indices of the two samples being compared.
## Columns 3-5 are a lower bound, estimate, and upper bound for their
## difference, where the bounds are for 95% confidence intervals. Column 6-8
## are the multiplicity adjusted p-values for each individual comparison, the
## test-statistic and the degrees of freedom. For @qcode{anovan}, the test
## statistic is the t-statistic.
##
## @qcode{multcompare} can take a number of optional parameters as name-value 
## pairs.
##
## @code{[@dots{}] = multcompare (@var{STATS}, "alpha", @var{ALPHA})}
##
## @itemize
## @item
## @var{ALPHA} sets the significance level of null hypothesis significance
## tests to ALPHA, and the central coverage of two-sided confidence intervals to
## 100*(1-@var{ALPHA})%. (Default ALPHA is 0.05).
## @end itemize
##
## @code{[@dots{}] = multcompare (@var{STATS}, "ControlGroup", @var{REF})}
##
## @itemize
## @item
## @var{REF} is the index of the control group to limit comparisons to. The
## index must be a positive integer scalar value. For each dimension (d) listed
## in @var{DIM}, multcompare uses STATS.grpnames@{d@}(idx) as the control group.
## (Default is empty, i.e. [], for full pairwise comparisons)
## @end itemize
##
## @code{[@dots{}] = multcompare (@var{STATS}, "ctype", @var{CTYPE})}
##
## @itemize
## @item
## @var{CTYPE} is the type of comparison test to use. In order of increasing power,
## the choices are: "bonferroni", "scheffe', "holm" (default), "fdr", "lsd". The
## first three control the family-wise error rate. The "fdr" method controls the
## false discovery rate. The final method, "lsd", is Fisher's least significant
## difference, which makes no attempt to control the Type 1 error rate of
## multiple comparisons. The coverage of confidence intervals are only corrected
## for multiple comparisons in the cases where CTYPE is "bonferroni" or "scheffe",
## where control of the Type 1 error rate is for simultaneous inference.
## @end itemize
##
## @code{[@dots{}] = multcompare (@var{STATS}, "dim", @var{DIM})}
##
## @itemize
## @item
## @var{DIM} is a vector specifying the dimension or dimensions over which the
## estimated marginal means are to be calculated. Used only if STATS comes from
## anovan. The value [1 3], for example, computes the estimated marginal mean
## for each combination of the first and third predictor values. The default is
## to compute over the first dimension (i.e. 1). If the specified dimension is,
## or includes, a continuous factor then @qcode{multcompare} will return an error.
## @end itemize
##
## @code{[@dots{}] = multcompare (@var{STATS}, "display", @var{DISPLAY})}
##
## @itemize
## @item
## @var{DISPLAY} is either "on" (the default) to display a graph of the comparisons
## (e.g. difference between means) and their 100*(1-@var{ALPHA})% intervals, or
## "off" to omit the graph. Markers and error bars colored red have multiplicity
## adjusted p-values < ALPHA, otherwise the markers and error bars are blue.
## @end itemize
##
## [@var{C}, @var{M}, @var{H}, @var{GNAMES}] = multcompare (@dots{}) returns
## additional outputs. @var{M} is a matrix where columns 1-2 are the estimated
## marginal means and their standard errors, and columns 3-4 are lower and upper
## bounds of the confidence intervals for the means; the critical value of the
## test statistic is scaled by a factor of 2^(-0.5) before multiplying by the
## standard errors of the group means so that the intervals overlap when the
## difference in means becomes significant at the level @var{ALPHA}. When
## @var{ALPHA} is 0.05, this corresponds to confidence intervals with 83.4%
## central coverage. @var{H} is a handle to the figure containing the graph.
## @var{GNAMES} is a cell array with one row for each group, containing the
## names of the groups.
##
## @seealso{anovan}
## @end deftypefn

function [C, M, H, GNAMES] = multcompare (STATS, varargin)

    if (nargin < 1)
      error (strcat (["multcompare usage: ""multcompare (STATS)""; "], ...
                      [" atleast 1 input arguments required"]));
    endif

    ## Check supplied parameters
    if ((numel (varargin) / 2) != fix (numel (varargin) / 2))
      error ("multcompare: wrong number of arguments.")
    endif
    ALPHA = 0.05;
    REF = [];
    CTYPE = "holm";
    DISPLAY = "on";
    DIM = 1;
    for idx = 3:2:nargin
      name = varargin{idx-2};
      value = varargin{idx-1};
      switch (lower (name))
        case "alpha"
          ALPHA = value;
        case "controlgroup"
          REF = value;
        case {"ctype","CriticalValueType"}
          CTYPE = lower (value);
        case "display"
          DISPLAY = lower (value);
        case {"dim","dimension"}
          DIM = value;
        otherwise
          error (sprintf ("multcompare: parameter %s is not supported", name));
      endswitch
    endfor

    ## Evaluate ALPHA input argument
    if (! isa (ALPHA,"numeric") || numel (ALPHA) != 1)
      error("anovan:alpha must be a numeric scalar value");
    endif
    if ((ALPHA <= 0) || (ALPHA >= 1))
      error("anovan: alpha must be a value between 0 and 1");
    endif

    ## Evaluate CTYPE input argument
    if (! ismember (lower (CTYPE), {"scheffe","bonferroni","holm","fdr","lsd"}))
      error ("multcompare: '%s' is not a supported value for CTYPE", CTYPE)
    endif

    switch (STATS.source)

      case "anovan"

        % Our calculations treat all effects as fixed
        if (! isempty (STATS.random))
          warning (strcat (["multcompare: ignoring random effects"], ... 
                           [" (all effects treated as fixed)"]));
        endif

        if (any (STATS.nlevels(DIM) < 2))
          error (strcat (["multcompare: DIM must specify only categorical"], ...
                         [" factors with 2 or more degrees of freedom."]));
        end

        ## Calculate estimated marginal means
        Nd = numel (DIM);
        n = numel (STATS.resid);
        df = STATS.df;
        dfe = STATS.dfe;
        i = 1 + cumsum(df);
        k = find (sum (STATS.terms(:,DIM), 2) == sum (STATS.terms, 2));
        Nb = 1 + sum(df(k));
        Nt = numel (k);
        L = zeros (n, sum (df) + 1);
        for j = 1:Nt
          L(:, i(k(j)) - df(k(j)) + 1 : i(k(j))) = STATS.X(:,i(k(j)) - ...
                                                   df(k(j)) + 1 : i(k(j)));
        endfor
        L(:,1) = 1;
        U = unique (L, "rows", "stable");
        Ng = size (U, 1);
        idx = zeros (Ng, 1);
        for k = 1:Ng
          idx(k) = find (all (L == U(k, :), 2),1);
        endfor
        gmeans = U * STATS.coeffs(:,1);     # Estimated marginal means
        gcov = U * STATS.vcov * U';
        gvar = diag (gcov);                 # Sampling variance 
        M = cat (2, gmeans, sqrt(gvar));

        ## Create cell array of group names corresponding to each row of m
        GNAMES = cell (Ng, 1);
        for i = 1:Ng
          str = "";
          for j = 1:Nd
            str = sprintf("%s%s=%s,", str, ...
                      num2str(STATS.varnames{DIM(j)}), ...
                      num2str(STATS.grpnames{DIM(j)}{STATS.grps(idx(i),DIM(j))}));
          endfor
          GNAMES{i} = str(1:end-1);
          str = "";
        endfor

        ## Make comparison matrix
        if (isempty (REF))
          ## Pairwise comparisons
          pairs = pairwise (Ng);
        else 
          ## Treatment vs. Control comparisons
          pairs = trt_vs_ctrl (Ng, REF);
        endif
        Np = size (pairs, 1);

        ## Calculate vector t-statistics corresponding to the comparisons. In
        ## balanced ANOVA designs, for the calculation of the t-statistics, the
        ## mean and standard error of the difference can be calculated simply by:
        ##      mean_diff = M(pairs(:,1) - M(pairs(:,2)
        ##      sed = sqrt (M(pairs(:,1),2).^2 + (M(pairs(:,2),2).^2))
        ##      t = mean_diff ./ sed
        ## However, to generalise the calculations for unbalanced N-way ANOVA
        ## we need to take into account correlations, so we use the covariance
        ## matrix of the estimated marginal means instead.
        L = zeros (Np, Ng);
        for j = 1:Np
          L(j, pairs(j,:)) = [1,-1];  # Hypothesis matrix
        endfor
        mean_diff = sum (L * diag (M(:, 1)), 2);
        sed = sqrt (diag (L * gcov * L'));
        t =  mean_diff ./ sed;

      otherwise

        error (strcat (sprintf ("multcompare: the STATS structure from %s", ...
               STATS.source), [" is not currently supported"]))

    endswitch

    ## The test specific code above needs to create the following variables in
    ## order to proceed with the remainder of the function tasks
    ## - Ng: number of groups involved in comparisons
    ## - M: Ng-by-2 matrix of group means (column 1) and standard errors (column 2)
    ## - Np: number of comparisons (pairs of groups being compaired)
    ## - pairs: Np-by-2 matrix of numeric group IDs - each row is a comparison
    ## - sed: vector containing SE of the difference for each comparisons
    ## - t: vector containing t for the difference relating to each comparisons
    ## - dfe: residual/error degrees of freedom
    ## - GNAMES: a cell array containing the names of the groups being compared

    ## If applicable, control simultaneous coverage of confidence intervals
    ## for multiple comparisons by modifying the critical value of the test
    ## statistic
    switch (lower (CTYPE)) 
      case "scheffe"
        critval = sqrt ((Ng - 1) * finv (1 - ALPHA, Ng - 1, dfe));
      case "bonferroni"
        ALPHA = ALPHA / Np;
        critval = tinv (1 - ALPHA / 2, dfe);
      otherwise
        ## No adjustment to confidence interval coverage
        critval = tinv (1 - ALPHA / 2, dfe);
    endswitch

    ## Create matrix of comparisons and calculate confidence intervals and
    ## multiplicity adjusted p-values for the comparisons
    C = zeros (Np, 7);
    C(:,1:2) = pairs;
    C(:,4) = (M(pairs(:, 1),1) - M(pairs(:, 2),1));
    C(:,3) = C(:,4) - sed * critval;
    C(:,5) = C(:,4) + sed * critval;
    p = 2 * (1 - tcdf (abs (t), dfe));
    C(:,6) = feval (CTYPE, p, t, Ng, dfe);
    C(:,7) = t;     # Unlike Matlab, we include the t-statistic
    C(:,8) = dfe;   # Unlike Matlab, we include the degrees of freedom

    ## Calculate confidence intervals of the estimated marginal means with
    ## central coverage such that the intervals start to overlap where the
    ## difference reaches a two-tailed p-value of ALPHA. When ALPHA is 0.05,
    ## central coverage is approximately 83.4%
    M(:,3) = M(:,1) - M(:,2) * critval / sqrt(2);
    M(:,4) = M(:,1) + M(:,2) * critval / sqrt(2);

    ## If requested, plot graph of the difference means for each comparison
    ## with central coverage of confidence intervals at 100*(1-alpha)%
    if (strcmp (DISPLAY, "on"))
      H = figure;
      plot ([0; 0], [0; Np + 1]',"k:"); # Plot vertical dashed line at 0 effect
      set (gca, "Ydir", "reverse")      # Flip y-axis direction
      ylim ([0.5, Np + 0.5]);           # Set y-axis limits
      hold on                           # Plot on the same axis
      for j = 1:Np
        if (C(j,6) < ALPHA)
          ## Plot marker for the difference in means 
          plot (C(j,4), j,"or","MarkerFaceColor", "r");
          ## Plot line for each confidence interval
          plot ([C(j,3), C(j,5)], j * ones(2,1), "r-");   
        else
          ## Plot marker for the difference in means 
          plot (C(j,4), j,"ob","MarkerFaceColor", "b");
          ## Plot line for each confidence interval
          plot ([C(j,3), C(j,5)], j * ones(2,1), "b-");
        endif
      endfor
      hold off
      xlabel (sprintf ("%g%% confidence interval for the difference",...
                       100 * (1 - ALPHA)));
      ylabel ("Row number in matrix of comparisons (C)");  
    endif

endfunction


## Posthoc comparisons

function pairs = pairwise (Ng)

  ## Create pairs matrix for pairwise comparisons
  gid = [1:Ng]';  # Create numeric group ID
  A = ones (Ng, 1) * gid';
  B = tril (gid * ones(1, Ng),-1);
  pairs = [A(:), B(:)];
  ridx = (pairs(:, 2) == 0);
  pairs(ridx, :) = [];

endfunction


function pairs = trt_vs_ctrl (Ng, REF)

  ## Create pairs matrix for comparisons with control (REF)
  gid = [1:Ng]';  # Create numeric group ID
  pairs = zeros (Ng - 1, 2);
  pairs(:, 1) = REF;
  pairs(:, 2) = gid(gid != REF);
        
endfunction


## Methods to control family-wise error rate in multiple comparisons

function padj = scheffe (p, t, Ng, dfe)
  
  padj = 1 - fcdf ((t.^2)/(Ng - 1), Ng - 1, dfe);

endfunction


function padj = bonferroni (p)
  
  ## Bonferroni procedure
  k = numel (p);
  padj = min (p * k, 1.0);

endfunction


function padj = holm (p)

  ## Holm's step-down Bonferroni procedure

  ## Order raw p-values
  [ps, idx] = sort (p, "ascend");
  m = numel(ps);

  ## Implement Holm's step-down Bonferroni procedure
  padj = nan (m,1);
  padj(1) = m * ps(1);
  for i = 2:m
    padj(i) = max (padj(i - 1), (m - i + 1) * ps(i));
  endfor

  ## Reorder the adjusted p-values to match the order of the original p-values
  [jnk, original_order] = sort (idx, "ascend");
  padj = padj(original_order);

  ## Truncate adjusted p-values to 1.0
  padj(padj>1) = 1;

endfunction


function padj = fdr (p)
  
  ## Benjamini-Hochberg procedure to control the false discovery rate (FDR)
  ## This procedure does not control the family-wise error rate

  ## Order raw p-values
  [ps, idx] = sort (p, "ascend");
  m = numel(ps);

  ## Initialize
  m = numel(p);
  padj = nan(m,1);
  alpha = nan(m,1);

  ## Benjamini-Hochberg step-up procedure to control the false discovery rate
  padj = nan (m,1);
  padj(m) = ps(m);
  for j = 1:m-1
    i = m - j;
    padj(i) = min (padj(i + 1), m / i * ps(i));
  endfor

  ## Reorder the adjusted p-values to match the order of the original p-values
  [jnk, original_order] = sort (idx, "ascend");
  padj = padj(original_order);
  
  ## Truncate adjusted p-values to 1.0
  padj(padj>1) = 1;

endfunction


function padj = lsd (p)
  
  ## Fisher's Least Significant Difference
  ## No control of the type I error rate across multiple comparisons
  padj = p;

endfunction


%!demo
%! 
%! # Demonstration using unbalanced one-way ANOVA example from anovan
%!
%! dv =  [ 8.706 10.362 11.552  6.941 10.983 10.092  6.421 14.943 15.931 ...
%!        22.968 18.590 16.567 15.944 21.637 14.492 17.965 18.851 22.891 ...
%!        22.028 16.884 17.252 18.325 25.435 19.141 21.238 22.196 18.038 ...
%!        22.628 31.163 26.053 24.419 32.145 28.966 30.207 29.142 33.212 ...
%!        25.694 ]';
%! g = [1 1 1 1 1 1 1 1 2 2 2 2 2 3 3 3 3 3 3 3 3 ...
%!      4 4 4 4 4 4 4 5 5 5 5 5 5 5 5 5]';
%!
%! [P,ATAB, STATS] = anovan (dv, g, "varnames", "score", "display", "off");
%!
%! [C, M, H, GNAMES] = multcompare (STATS, "dim", 1, "ctype", "holm", ...
%!                                  "ControlGroup", 1, "display", "on")
%!

%!demo
%! 
%! # Demonstration using factorial ANCOVA example from anovan
%!
%! score = [95.6 82.2 97.2 96.4 81.4 83.6 89.4 83.8 83.3 85.7 ...
%! 97.2 78.2 78.9 91.8 86.9 84.1 88.6 89.8 87.3 85.4 ...
%! 81.8 65.8 68.1 70.0 69.9 75.1 72.3 70.9 71.5 72.5 ...
%! 84.9 96.1 94.6 82.5 90.7 87.0 86.8 93.3 87.6 92.4 ...
%! 100. 80.5 92.9 84.0 88.4 91.1 85.7 91.3 92.3 87.9 ...
%! 91.7 88.6 75.8 75.7 75.3 82.4 80.1 86.0 81.8 82.5]';
%! treatment = {"yes" "yes" "yes" "yes" "yes" "yes" "yes" "yes" "yes" "yes" ...
%!              "yes" "yes" "yes" "yes" "yes" "yes" "yes" "yes" "yes" "yes" ...
%!              "yes" "yes" "yes" "yes" "yes" "yes" "yes" "yes" "yes" "yes" ...
%!              "no"  "no"  "no"  "no"  "no"  "no"  "no"  "no"  "no"  "no"  ...
%!              "no"  "no"  "no"  "no"  "no"  "no"  "no"  "no"  "no"  "no"  ...
%!              "no"  "no"  "no"  "no"  "no"  "no"  "no"  "no"  "no"  "no"}';
%! exercise = {"lo"  "lo"  "lo"  "lo"  "lo"  "lo"  "lo"  "lo"  "lo"  "lo"  ...
%!             "mid" "mid" "mid" "mid" "mid" "mid" "mid" "mid" "mid" "mid" ...
%!             "hi"  "hi"  "hi"  "hi"  "hi"  "hi"  "hi"  "hi"  "hi"  "hi"  ...
%!             "lo"  "lo"  "lo"  "lo"  "lo"  "lo"  "lo"  "lo"  "lo"  "lo"  ...
%!             "mid" "mid" "mid" "mid" "mid" "mid" "mid" "mid" "mid" "mid" ...
%!             "hi"  "hi"  "hi"  "hi"  "hi"  "hi"  "hi"  "hi"  "hi"  "hi"}';
%! age = [59 65 70 66 61 65 57 61 58 55 62 61 60 59 55 57 60 63 62 57 ...
%! 58 56 57 59 59 60 55 53 55 58 68 62 61 54 59 63 60 67 60 67 ...
%! 75 54 57 62 65 60 58 61 65 57 56 58 58 58 52 53 60 62 61 61]';
%!
%! [P, ATAB, STATS] = anovan (score, {treatment, exercise, age}, "model", ...
%!                            [1 0 0; 0 1 0; 0 0 1; 1 1 0], "continuous", 3, ...
%!                            "sstype", "h", "display", "off", "contrasts", ...
%!                            {"simple","poly",""});
%!
%! [C, M, H, GNAMES] = multcompare (STATS, "dim", [1 2], "ctype", "holm", ...
%!                                  "display", "on")
%!


%!test
%! 
%! # Test using unbalanced one-way ANOVA example from anovan
%! # Comparison to matlab
%!
%! dv =  [ 8.706 10.362 11.552  6.941 10.983 10.092  6.421 14.943 15.931 ...
%!        22.968 18.590 16.567 15.944 21.637 14.492 17.965 18.851 22.891 ...
%!        22.028 16.884 17.252 18.325 25.435 19.141 21.238 22.196 18.038 ...
%!        22.628 31.163 26.053 24.419 32.145 28.966 30.207 29.142 33.212 ...
%!        25.694 ]';
%! g = [1 1 1 1 1 1 1 1 2 2 2 2 2 3 3 3 3 3 3 3 3 ...
%!      4 4 4 4 4 4 4 5 5 5 5 5 5 5 5 5]';
%!
%! [P,ATAB, STATS] = anovan (dv, g, 'varnames', 'score', 'display', 'off');
%!
%! [C, M, ~, GNAMES] = multcompare (STATS, 'dim', 1, 'ctype', 'lsd', ...
%!                                  'display', 'off');
%! assert (C(1,6), 2.85812420217898e-05, 1e-09);
%! assert (C(2,6), 5.22936741204085e-07, 1e-09);
%! assert (C(3,6), 2.12794763209146e-08, 1e-09);
%! assert (C(4,6), 7.82091664406946e-15, 1e-09);
%! assert (C(5,6), 0.546591417210693, 1e-09);
%! assert (C(6,6), 0.0845897945254446, 1e-09);
%! assert (C(7,6), 9.47436557975328e-08, 1e-09);
%! assert (C(8,6), 0.188873478781067, 1e-09);
%! assert (C(9,6), 4.08974010364197e-08, 1e-09);
%! assert (C(10,6), 4.44427348175241e-06, 1e-09);
%! assert (M(1,1), 10, 1e-09);
%! assert (M(2,1), 18, 1e-09);
%! assert (M(3,1), 19, 1e-09);
%! assert (M(4,1), 21.0001428571429, 1e-09);
%! assert (M(5,1), 29.0001111111111, 1e-09);
%! assert (M(1,2), 1.0177537954095, 1e-09);
%! assert (M(2,2), 1.28736803631001, 1e-09);
%! assert (M(3,2), 1.0177537954095, 1e-09);
%! assert (M(4,2), 1.0880245732889, 1e-09);
%! assert (M(5,2), 0.959547480416536, 1e-09);

