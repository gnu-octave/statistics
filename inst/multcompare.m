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
## @code{@var{C} = multcompare (@var{stats})} performs a multiple comparison
## using a @var{STATS} structure that is obtained as output from any of
## the following functions:  anovan.
## The return value @var{C} is a matrix with one row per comparison and six
## columns. Columns 1-2 are the indices of the two samples being compared.
## Columns 3-5 are a lower bound, estimate, and upper bound for their difference,
## where the bounds are for confidence intervals. Column 6 is the multiplicity
## adjusted p-value for each individual comparison.
##
## @qcode{multcompare} can take a number of optional parameters as name-value 
## pairs.
##
## @code{[@dots{}] = multcompare (@var{STATS}, "alpha", @var{ALPHA})}
##
## @itemize
## @item
## @var{ALPHA} specifies the significance level as ALPHA, and the central 
## coverage of two-sided confidence intervals as 100*(1-@var{ALPHA})% 
## (default ALPHA is 0.05).
## @end itemize
##
## @code{[@dots{}] = multcompare (@var{STATS}, "ctype", @var{DISPLAY})}
##
## @itemize
## @item
## @var{DISPLAY} is either "on" (the default) to display a graph of the comparisons
## (e.g. difference between means) and their 100*(1-@var{ALPHA})% intervals, or
## "off" to omit the graph. Markers and error bars colored red have multiplicity
## adjusted p-values < ALPHA, otherise the markers and error bars are blue.
## @end itemize
##
## @code{[@dots{}] = multcompare (@var{STATS}, "ctype", @var{CTYPE})}
##
## @itemize
## @item
## @var{CTYPE} is the type of comparison test to use. In order of increasing power,
## the choices are "scheffe', "bonferroni", "holm" (default), "fdr", "lsd". The
## first three control the family-wise error rate. The "fdr" method contraols
## false discovery rate. The final method, "lsd" is Fisher's least significant
## difference, which makes no attempt to control the Type 1 error rate of
## multiple comparisons. The coverage of confidence intervals are only corrected
## for multiple comparisons in the case where CTYPE is "scheffe" or "bonferroni".
## @end itemize
##
## @code{[@dots{}] = multcompare (@var{STATS}, "dim", @var{DIM})}
##
## @itemize
## @item
## @var{DIM} is a vector specifying the dimension or dimensions over which the
## estimated marginal means are to be calculated. Used only if STATS comes from
## anovan. The default is to compute over the first dimension associated with a
## categorical (non-continuous) factor. The value [1 3], for example, computes
## the estimated marginal mean for each combination of the first and third
## predictor values.
## @end itemize
##
## @code{[@dots{}] = multcompare (@var{STATS}, "ControlGroup", @var{REF})}
##
## @itemize
## @item
## @var{REF} is the index of the control group for Dunnett's test, specified as 
## a positive integer value. For each dimension (d) listed in @var{DIM},
## multcompare uses STATS.grpnames@{d@}(idx) as the control group.
## @end itemize
##
## [@var{C}, @var{M}, @var{H}, @var{GNAMES}] = multcompare (@dots{}) returns
## additional outputs. @var{M} is a matrix with columns equal to the estimated
## marginal means and their standard errors. @var{H} is a handle to the figure
## containing the graph. @var{GNAMES} is a cell array with one row for each
## group, containing the names of the groups.
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
          CTYPE = value;
        case "display"
          DISPLAY = value;
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
    if (~ismember (lower (CTYPE),{"scheffe","bonferroni","holm","fdr","lsd"}))
      error ("multcompare: '%s' is not a supported value for CTYPE", CTYPE)
    endif

    switch (STATS.source)

      case "anovan"

        % Our calculations treat all effects as fixed
        if (~isempty (STATS.ems))
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
          L(:, i(k(j))-df(k(j)) + 1 : i(k(j))) = STATS.X(:,i(k(j)) - ...
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
    
        ## Calculate vector t-statistics corresponding to the comparisons
        L = zeros (Np, Ng);
        t = nan (Np, 1);
        sed = nan (Np, 1);
        for j = 1:Np
          L(j, pairs(j,:)) = [1,-1];
          sed(j) = sqrt (L(j,:) * gcov * L(j,:)');
          t(j) = (M(pairs(j,1),1) - M(pairs(j,2),1)) / sed(j);
        endfor

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
        critval = sqrt ((Ng - 1) * finv (1 - ALPHA, Ng-1, dfe));
      case "bonferroni"
        ALPHA = ALPHA / Np;
        critval = tinv (1 - ALPHA / 2, STATS.dfe);
      otherwise
        ## No adjustment
        critval = tinv (1 - ALPHA / 2, STATS.dfe);
    endswitch

    ## reate comparisons matrix and calculate confidence intervals and
    ## multiplicity adjusted p-values for the comparisons
    C = zeros (Np, 6);
    C(:,1:2) = pairs;
    C(:,4) = (M(pairs(:, 1),1) - M(pairs(:, 2),1));
    C(:,3) = C(:,4) - sed * critval;
    C(:,5) = C(:,4) + sed * critval;
    p = 2 * (1 - tcdf (abs (t), dfe));
    C(:,6) = feval (CTYPE, p, t, Ng, dfe);

    ## If requested, plot graph of the difference means for each comparison
    ## with central coverage of confidence intervals at 100*(1-alpha)%
    if (strcmpi (DISPLAY, "on"))
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
      ylabel ("Comparison matrix (C) row number");  
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
  pairs(:, 2) = gid(gid ~= REF);
        
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
    i = m-j;
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
