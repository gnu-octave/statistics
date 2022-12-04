## Copyright (C) 2022 Andrew Penn <A.C.Penn@sussex.ac.uk>
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
## @deftypefn {Function File} @var{C} = multcompare (@var{STATS})
## @deftypefnx {Function File} @var{C} = multcompare (@var{STATS}, "name", @var{value})
## @deftypefnx {Function File} [@var{C}, @var{M}] = multcompare (...)
## @deftypefnx {Function File} [@var{C}, @var{M}, @var{H}] = multcompare (...)
## @deftypefnx {Function File} [@var{C}, @var{M}, @var{H}, @var{GNAMES}] = multcompare (...)
## @deftypefnx {Function File} @var{padj} = multcompare (@var{p})
## @deftypefnx {Function File} @var{padj} = multcompare (@var{p}, "ctype", @var{CTYPE})
##
## Perform posthoc multiple comparison tests or p-value adjustments to control
## the family-wise error rate (FWER) or false discovery rate (FDR).
##
## @code{@var{C} = multcompare (@var{STATS})} performs a multiple comparison
## using a @var{STATS} structure that is obtained as output from any of
## the following functions: anova1, anova2, anovan, kruskalwallis, and friedman.
## The return value @var{C} is a matrix with one row per comparison and six
## columns. Columns 1-2 are the indices of the two samples being compared.
## Columns 3-5 are a lower bound, estimate, and upper bound for their
## difference, where the bounds are for 95% confidence intervals. Column 6-8 are
## the multiplicity adjusted p-values for each individual comparison, the test
## statistic and the degrees of freedom.
## All tests by multcompare are two-tailed.
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
## @var{CTYPE} is the type of comparison test to use. In order of increasing
## power, the choices are: "bonferroni", "scheffe", "mvt", "holm" (default),
## "hochberg", "fdr", "lsd". The first five methods control the family-wise
## error rate.  The "fdr" method controls false discovery rate (by the original
## Benjamini-Hochberg step-up procedure). The final method, "lsd" (or "none"),
## makes no attempt to control the Type 1 error rate of multiple comparisons. 
## The coverage of confidence intervals are only corrected for multiple
## comparisons in the cases where @var{CTYPE} is "bonferroni", "scheffe" or
## "mvt", which control the Type 1 error rate for simultaneous inference.
##
## The "mvt" method uses the multivariate t distribution to assess the
## probability or critical value of the maximum statistic across the tests,
## thereby accounting for correlations among comparisons in the control of the
## family-wise error rate with simultaneous inference.  In the case of pairwise
## comparisons, it simulates Tukey's (or the Games-Howell) test, in the case of
## comparisons with a single control group, it simulates Dunnett's test.
## @var{CTYPE} values "tukey-kramer" and "hsd" are recognised but set the value
## of @var{CTYPE} and @var{REF} to "mvt" and empty respectively.  A @var{CTYPE}
## value "dunnett" is recognised but sets the value of @var{CTYPE} to "mvt", and
## if @var{REF} is empty, sets @var{REF} to 1.  Since the algorithm uses a Monte
## Carlo method (of 1e+06 random samples), you can expect the results to
## fluctuate slightly with each call to multcompare and the calculations may be
## slow to complete for a large number of comparisons.  If the parallel package
## is installed and loaded, @qcode{multcompare} will automatically accelerate
## computations by parallel processing.  Note that p-values calculated by the
## "mvt" are truncated at 1e-06.
## @end itemize
##
## @code{[@dots{}] = multcompare (@var{STATS}, "df", @var{DF})}
##
## @itemize
## @item
## @var{DF} is an optional scalar value to set the number of degrees of freedom
## in the calculation of p-values for the multiple comparison tests. By default,
## this value is extracted from the @var{STATS} structure of the ANOVA test, but
## setting @var{DF} maybe necessary to approximate Satterthwaite correction if
## @qcode{anovan} was performed using weights.
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
## or includes, a continuous factor then @qcode{multcompare} will return an
## error.
## @end itemize
##
## @code{[@dots{}] = multcompare (@var{STATS}, "estimate", @var{ESTIMATE})}
##
## @itemize
## @item
## @var{ESTIMATE} is a string specifying the estimates to be compared when
## computing multiple comparisons after anova2; this argument is ignored by
## anovan and anova1. Accepted values for @var{ESTIMATE} are either "column"
## (default) to compare column means, or "row" to compare row means. If the
## model type in anova2 was "linear" or "nested" then only "column" is accepted
## for @var{ESTIMATE} since the row factor is assumed to be a random effect.
## @end itemize
##
## @code{[@dots{}] = multcompare (@var{STATS}, "display", @var{DISPLAY})}
##
## @itemize
## @item
## @var{DISPLAY} is either "on" (the default) to display a graph of the
## comparisons (e.g. difference between means) and their 100*(1-@var{ALPHA})%
## intervals, or "off" to omit the graph.  Markers and error bars colored red
## have multiplicity adjusted p-values < ALPHA, otherwise the markers and error
## bars are blue.
## @end itemize
##
## @code{[@var{C}, @var{M}, @var{H}, @var{GNAMES}] = multcompare (@dots{})}
## returns additional outputs. @var{M} is a matrix where columns 1-2 are the
## estimated marginal means and their standard errors, and columns 3-4 are lower
## and upper bounds of the confidence intervals for the means; the critical
## value of the test statistic is scaled by a factor of 2^(-0.5) before
## multiplying by the standard errors of the group means so that the intervals
## overlap when the difference in means becomes significant at approximately
## the level @var{ALPHA}. When @var{ALPHA} is 0.05, this corresponds to
## confidence intervals with 83.4% central coverage.  @var{H} is a handle to the
## figure containing the graph.  @var{GNAMES} is a cell array with one row for
## each group, containing the names of the groups.
##
## @code{@var{padj} = multcompare (@var{p})} calculates and returns adjusted
## p-values (@var{padj}) using the Holm-step down Bonferroni procedure to
## control the family-wise error rate.
##
## @code{@var{padj} = multcompare (@var{p}, "ctype", @var{CTYPE})} calculates
## and returns adjusted p-values (@var{padj}) computed using the method
## @var{CTYPE}. In order of increasing power, @var{CTYPE} for p-value adjustment
## can be either "bonferroni", "holm" (default), "hochberg", and "fdr". See
## above for further information about the @var{CTYPE} methods.
##
## @seealso{anova1, anova2, anovan, kruskalwallis, friedman, fitlm}
## @end deftypefn

function [C, M, H, GNAMES] = multcompare (STATS, varargin)

    if (nargin < 1)
      error (strcat (["multcompare usage: ""multcompare (ARG)""; "], ...
                      [" atleast 1 input argument required"]));
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
    ESTIMATE = "column";
    DFE = [];
    for idx = 3:2:nargin
      name = varargin{idx-2};
      value = varargin{idx-1};
      switch (lower (name))
        case "alpha"
          ALPHA = value;
        case {"controlgroup","ref"}
          REF = value;
        case {"ctype","CriticalValueType"}
          CTYPE = lower (value);
        case "display"
          DISPLAY = lower (value);
        case {"dim","dimension"}
          DIM = value;
        case "estimate"
          ESTIMATE = lower (value);
        case {"df","dfe"}
          DFE = value;
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
    if (ismember (CTYPE, {"tukey-kramer", "hsd"}))
      CTYPE = "mvt";
      REF = [];
    elseif (strcmp (CTYPE, "dunnett"))
      CTYPE = "mvt";
      if (isempty (REF))
        REF = 1;
      endif
    elseif (strcmp (CTYPE, "none"))
      CTYPE = "lsd";
    endif
    if (! ismember (CTYPE, ...
                    {"bonferroni","scheffe","mvt","holm","hochberg","fdr","lsd"}))
      error ("multcompare: '%s' is not a supported value for CTYPE", CTYPE)
    endif

    ## Evaluate DFE input argument
    if (! isempty (DFE))
      if (! isscalar (DFE))
        error ("multcompare: df must be a scalar value.");
      endif
      if (!(DFE > 0) || isinf (DFE))
        error ("multcompare: df must be a positive finite value.");
      endif
    endif

    ## If STATS is numeric, assume it is a vector of p-values
    if (isnumeric (STATS))
      if (! ismember (CTYPE, {"bonferroni","holm","hochberg","fdr"}))
        error ("multcompare: '%s' is not a supported p-adjustment method", CTYPE)
      endif
      p = STATS;
      if (all (size (p) > 1))
        error ("multcompare: p-values must be a vector")
      endif
      padj = feval (CTYPE, p);
      if (size (p, 1) > 1)
        C = padj;
      else 
        C = padj';
      endif
      return
    endif

    ## Perform test specific calculations
    switch (STATS.source)

      case "anova1"

        ## Make matrix of requested comparisons (pairs)
        ## Also return the corresponding hypothesis matrix (L)
        n = STATS.n(:);
        Ng = numel (n);
         if (isempty (REF))
          ## Pairwise comparisons
          [pairs, L] = pairwise (Ng);
        else
          ## Treatment vs. Control comparisons
          [pairs, L] = trt_vs_ctrl (Ng, REF);
        endif
        Np = size (pairs, 1);

        switch (STATS.vartype)

          case "equal"

            ## Calculate estimated marginal means and their standard errors
            gmeans = STATS.means(:);
            gvar = (STATS.s^2) ./ n;       # Sampling variance
            gcov = diag (gvar);
            Ng = numel (gmeans);
            M = zeros (Ng, 4);
            M(:,1:2)  = cat (2, gmeans, sqrt(gvar));

            ## Get the error degrees of freedom from anova1 output
            if (isempty (DFE))
              DFE = STATS.df;
            endif

          case "unequal"

            ## Error checking
            if (strcmp (CTYPE, "scheffe"))
              error (strcat (["multcompare: the CTYPE value 'scheffe'"], ...
                 [" does not support tests with varying degrees of freedom "]));
            endif

            ## Calculate estimated marginal means and their standard errors
            gmeans = STATS.means(:);
            gvar = STATS.vars(:) ./ n;     # Sampling variance
            gcov = diag (gvar);
            Ng = numel (gmeans);
            M = zeros (Ng, 4);
            M(:,1:2) = cat (2, gmeans, sqrt (gvar));

            ## Calculate Welch's corrected degrees of freedom
            if (isempty (DFE))
              DFE = sum (gvar(pairs), 2).^2 ./ ...
                    sum ((gvar(pairs).^2 ./ (n(pairs) - 1)), 2);
            endif

        endswitch

        ## Calculate t statistics corresponding to the comparisons defined in L
        [mean_diff, sed, t] = tValue (gmeans, gcov, L);

        ## Calculate correlation matrix
        vcov = L * gcov * L';
        R = cov2corr (vcov);

        ## Create cell array of group names corresponding to each row of m
        GNAMES = STATS.gnames;

      case "anova2"

        ## Fetch estimate specific information from the STATS structure
        switch (ESTIMATE)
          case {"column","columns","col","cols"}
            gmeans = STATS.colmeans(:);
            Ng = numel (gmeans);
            n = STATS.coln;
          case {"row","rows"}
            if (ismember (STATS.model, {"linear","nested"}))
              error (strcat (["multcompare: no support for the row factor"],...
                             [" (random effect) in a 'nested' or 'linear' anova2 model"]));
            endif
            gmeans = STATS.rowmeans(:);
            Ng = numel (gmeans);
            n = STATS.rown;
        endswitch

        ## Make matrix of requested comparisons (pairs)
        ## Also return the corresponding hypothesis matrix (L)
         if (isempty (REF))
          ## Pairwise comparisons
          [pairs, L, R] = pairwise (Ng);
        else
          ## Treatment vs. Control comparisons
          [pairs, L, R] = trt_vs_ctrl (Ng, REF);
        endif
        Np = size (pairs, 1);

        ## Calculate estimated marginal means and their standard errors
        gvar = ((STATS.sigmasq) / n) * ones (Ng, 1);  # Sampling variance
        gcov = diag (gvar);
        M = zeros (Ng, 4);
        M(:,1:2) = cat (2, gmeans, sqrt (gvar));

        ## Get the error degrees of freedom from anova2 output
        if (isempty (DFE))
          DFE = STATS.df;
        endif

        ## Calculate t statistics corresponding to the comparisons defined in L
        [mean_diff, sed, t] = tValue (gmeans, gcov, L);

        ## Create character array of group names corresponding to each row of m
        GNAMES = cellstr (num2str ([1:Ng]'));

      case {"anovan","fitlm"}

        ## Our calculations treat all effects as fixed
        if (ismember (STATS.random, DIM))
          warning (strcat (["multcompare: ignoring random effects"], ...
                           [" (all effects treated as fixed)"]));
        endif

        ## Check what type of factor is requested in DIM
        if (any (STATS.nlevels(DIM) < 2))
          error (strcat (["multcompare: DIM must specify only categorical"], ...
                         [" factors with 2 or more degrees of freedom."]));
        endif

        ## Check that all continuous variables were centered
        if (any (STATS.continuous - STATS.center_continuous))
          error (strcat (["use a STATS structure from a model refit with"], ...
                         [" sum-to-zero contrast coding, e.g. ""simple"""]))
        endif

        ## Calculate estimated marginal means and their standard errors
        Nd = numel (DIM);
        n = numel (STATS.resid);
        df = STATS.df;
        if (isempty (DFE))
          DFE = STATS.dfe;
        endif
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
        M = zeros (Ng, 4);
        M(:,1:2) = cat (2, gmeans, sqrt(gvar));

        ## Create cell array of group names corresponding to each row of m
        GNAMES = cell (Ng, 1);
        for i = 1:Ng
          str = "";
          for j = 1:Nd
            str = sprintf("%s%s=%s, ", str, ...
                      num2str(STATS.varnames{DIM(j)}), ...
                      num2str(STATS.grpnames{DIM(j)}{STATS.grps(idx(i),DIM(j))}));
          endfor
          GNAMES{i} = str(1:end-2);
          str = "";
        endfor

        ## Make matrix of requested comparisons (pairs)
        ## Also return the corresponding hypothesis matrix (L)
        if (isempty (REF))
          ## Pairwise comparisons
          [pairs, L] = pairwise (Ng);
        else
          ## Treatment vs. Control comparisons
          [pairs, L] = trt_vs_ctrl (Ng, REF);
        endif
        Np = size (pairs, 1);

        ## Calculate t statistics corresponding to the comparisons defined in L
        [mean_diff, sed, t] = tValue (gmeans, gcov, L);

        ## Calculate correlation matrix.
        vcov = L * gcov * L';
        R = cov2corr (vcov);

      case "friedman"

        ## Get stats from structure
        gmeans = STATS.meanranks(:);
        Ng = length (gmeans);
        sigma = STATS.sigma;

        ## Make group names
        GNAMES = strjust (num2str ((1:Ng)'), "left");

        ## Make matrix of requested comparisons (pairs)
        ## Also return the corresponding hypothesis matrix (L)
        if (isempty (REF))
          ## Pairwise comparisons
          [pairs, L, R] = pairwise (Ng);
        else
          ## Treatment vs. Control comparisons
          [pairs, L, R] = trt_vs_ctrl (Ng, REF);
        endif
        Np = size (pairs, 1);

        ## Calculate covariance matrix
        gcov = ((sigma ^ 2) / STATS.n) * eye (Ng);

        ## Create matrix with group means and standard errors
        M = cat (2, gmeans, sqrt (diag (gcov)));

        ## Calculate t statistics corresponding to the comparisons defined in L
        [mean_diff, sed, t] = tValue (gmeans, gcov, L); # z-statistic (not t)

        ## Calculate degrees of freedom from number of groups
        if (isempty (DFE))
          DFE = inf;  # this is a z-statistic so infinite degrees of freedom
        endif

      case "kruskalwallis"

        ## Get stats from structure
        gmeans = STATS.meanranks(:);
        sumt = STATS.sumt;
        Ng = length (gmeans);
        n = STATS.n(:);
        N = sum (n);

        ## Make group names
        GNAMES = STATS.gnames;
        ## Make matrix of requested comparisons (pairs)
        ## Also return the corresponding hypothesis matrix (L)
        if (isempty (REF))
          ## Pairwise comparisons
          [pairs, L] = pairwise (Ng);
        else
          ## Treatment vs. Control comparisons
          [pairs, L] = trt_vs_ctrl (Ng, REF);
        endif
        Np = size (pairs, 1);

        ## Calculate covariance matrix
        gcov = diag (((N * (N + 1) / 12) - (sumt / (12 * (N - 1)))) ./ n);

        ## Create matrix with group means and standard errors
        M = cat (2, gmeans, sqrt (diag (gcov)));

        ## Calculate t statistics corresponding to the comparisons defined in L
        [mean_diff, sed, t] = tValue (gmeans, gcov, L); # z-statistic (not t)

        ## Calculate correlation matrix
        vcov = L * gcov * L';
        R = cov2corr (vcov);

        ## Calculate degrees of freedom from number of groups
        if (isempty (DFE))
          DFE = inf;  # this is a z-statistic so infinite degrees of freedom
        endif

      otherwise

        error (strcat (sprintf ("multcompare: the STATS structure from %s", ...
               STATS.source), [" is not currently supported"]))

    endswitch

    ## The test specific code above needs to create the following variables in
    ## order to proceed with the remainder of the function tasks
    ## - Ng: number of groups involved in comparisons
    ## - M: Ng-by-2 matrix of group means (col 1) and standard errors (col 2)
    ## - Np: number of comparisons (pairs of groups being compaired)
    ## - pairs: Np-by-2 matrix of numeric group IDs - each row is a comparison
    ## - R: correlation matrix for the requested comparisons
    ## - sed: vector containing SE of the difference for each comparisons
    ## - t: vector containing t for the difference relating to each comparisons
    ## - DFE: residual/error degrees of freedom
    ## - GNAMES: a cell array containing the names of the groups being compared

    ## Create matrix of comparisons and calculate confidence intervals and
    ## multiplicity adjusted p-values for the comparisons.
    C = zeros (Np, 7);
    C(:,1:2) = pairs;
    C(:,4) = (M(pairs(:, 1),1) - M(pairs(:, 2),1));
    C(:,7) = t;     # Unlike Matlab, we include the t statistic
    C(:,8) = DFE;   # Unlike Matlab, we include the degrees of freedom
    if (any (isinf (DFE)))
      p = 2 * (1 - normcdf (abs (t)));
    else
      p = 2 * (1 - tcdf (abs (t), DFE));
    endif
    [C(:,6), critval, C(:,8)] = feval (CTYPE, p, t, Ng, DFE, R, ALPHA);
    C(:,3) = C(:,4) - sed .* critval;
    C(:,5) = C(:,4) + sed .* critval;

    ## Calculate confidence intervals of the estimated marginal means with
    ## central coverage such that the intervals start to overlap where the
    ## difference reaches a two-tailed p-value of ALPHA. When ALPHA is 0.05,
    ## central coverage is approximately 83.4%
    if (! isscalar(DFE))
      # Upper bound critval (corresponding to lower bound DFE)
      critval = max (critval);
    endif
    M(:,3) = M(:,1) - M(:,2) .* critval / sqrt(2);
    M(:,4) = M(:,1) + M(:,2) .* critval / sqrt(2);

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
    else
      H = [];
    endif

endfunction


## Posthoc comparisons

function [pairs, L, R] = pairwise (Ng)

  ## Create pairs matrix for pairwise comparisons
  gid = [1:Ng]';  # Create numeric group ID
  A = ones (Ng, 1) * gid';
  B = tril (gid * ones(1, Ng),-1);
  pairs = [A(:), B(:)];
  ridx = (pairs(:, 2) == 0);
  pairs(ridx, :) = [];

  ## Calculate correlation matrix (required for CTYPE "mvt")
  Np = size (pairs, 1);
  L = zeros (Np, Ng);
  for j = 1:Np
    L(j, pairs(j,:)) = [1,-1];  # Hypothesis matrix
  endfor
  R = corr (L'); # Correlation matrix

endfunction


function [pairs, L, R] = trt_vs_ctrl (Ng, REF)

  ## Create pairs matrix for comparisons with control (REF)
  gid = [1:Ng]';  # Create numeric group ID
  pairs = zeros (Ng - 1, 2);
  pairs(:, 1) = REF;
  pairs(:, 2) = gid(gid != REF);

  ## Calculate correlation matrix (required for CTYPE "mvt")
  Np = size (pairs, 1);
  L = zeros (Np, Ng);
  for j = 1:Np
    L(j, pairs(j,:)) = [1,-1];  # Hypothesis matrix
  endfor
  R = corr (L'); # Correlation matrix

endfunction

function [mn, se, t] = tValue (gmeans, gcov, L)

  ## Calculate means, standard errors and t (or z) statistics
  ## corresponding to the comparisons defined in L.
  mn = sum (L * diag (gmeans), 2);
  se = sqrt (diag (L * gcov * L'));
  t =  mn ./ se;

endfunction

function R = cov2corr (vcov)

   ## Convert covariance matrix to correlation matrix
   sed = sqrt (diag (vcov));
   R = vcov ./ (sed * sed');
   R = (R + R') / 2; # This step ensures that the matrix is positive definite

endfunction


## Methods to control family-wise error rate in multiple comparisons

function [padj, critval, dfe] = scheffe (p, t, Ng, dfe, R, ALPHA)

  ## Calculate the p-value
  if (isinf (dfe))
    padj = 1 - chi2cdf (t.^2, Ng - 1);
  else
    padj = 1 - fcdf ((t.^2) / (Ng - 1), Ng - 1, dfe);
  endif

  ## Calculate critical value at Scheffe-adjusted ALPHA level
  if (isinf (dfe))
    tmp = chi2inv (1 - ALPHA, Ng - 1) / (Ng - 1);
  else
    tmp = finv (1 - ALPHA, Ng - 1, dfe);
  end
  critval = sqrt ((Ng - 1) * tmp);

endfunction


function [padj, critval, dfe] = bonferroni (p, t, Ng, dfe, R, ALPHA)

  ## Bonferroni procedure
  Np = numel (p);
  padj = min (p * Np, 1.0);

  ## If requested, calculate critical value at Bonferroni-adjusted ALPHA level
  if (nargout > 1)
    critval = tinv (1 - ALPHA / Np * 0.5, dfe);
  endif

endfunction

function [padj, critval, dfe] = mvt (p, t, Ng, dfe, R, ALPHA)

  ## Monte Carlo simulation of the maximum test statistic in random samples
  ## generated from a multivariate t distribution. This method accounts for
  ## correlations among comparisons. This method simulates Tukey's test in the
  ## case of pairwise comparisons or Dunnett's tests in the case of trt_vs_ctrl.
  ## The "mvt" method is equivalent to methods used in the following R packages:
  ##   - emmeans: the "mvt" adjust method in functions within emmeans
  ##   - glht: the "single-step" adjustment in the multcomp.function

  ## Lower bound for error degrees of freedom to ensure type 1 error rate isn't
  ## exceeded for any test
  if (! isscalar(dfe))
    dfe = max (1, round (min (dfe)));
    fprintf ("Note: df set to %u (lower bound)\n", dfe);
  endif

  ## Check if we can use parallel processing to accelerate computations
  pat = '^parallel';
  software = pkg('list');
  names = cellfun (@(S) S.name, software, 'UniformOutput', false);
  status = cellfun (@(S) S.loaded, software, 'UniformOutput', false);
  index = find (! cellfun (@isempty, regexpi (names, pat)));
  if (! isempty (index))
    if (logical (status{index}))
      PARALLEL = true;
    else
      PARALLEL = false;
    endif
  else
    PARALLEL = false;
  endif

  ## Generate the distribution of (correlated) t statistics under the null, and
  ## calculate the maximum test statistic for each random sample. Computations
  ## are performed in chunks to prevent memory issues when the number of
  ## comparisons is large.
  chunkSize = 1000;
  numChunks = 1000;
  nsim = chunkSize * numChunks;
  if (isinf (dfe))
    # Multivariate z-statistics
    func = @(jnk) max (abs (mvnrnd (0, R, chunkSize)'), [], 1);
  else
    # Multivariate t-statistics
    func = @(jnk) max (abs (mvtrnd (R, dfe, chunkSize)'), [], 1);
  endif
  if (PARALLEL)
    maxT = cell2mat (parcellfun (nproc, func, ...
                                 cell (1, numChunks), 'UniformOutput', false));
  else
    maxT = cell2mat (cellfun (func, cell (1, numChunks), 'UniformOutput', false));
  endif

  ## Calculate multiplicity adjusted p-values (two-tailed)
  padj = max (sum (bsxfun (@ge, maxT, abs (t)), 2) / nsim, nsim^-1);

  ## Calculate critical value adjusted by the maxT procedure
  critval = quantile (maxT, 1 - ALPHA);

endfunction


function [padj, critval, dfe] = holm (p, t, Ng, dfe, R, ALPHA)

  ## Holm's step-down Bonferroni procedure

  ## Order raw p-values
  [ps, idx] = sort (p, "ascend");
  Np = numel (ps);

  ## Implement Holm's step-down Bonferroni procedure
  padj = nan (Np,1);
  padj(1) = Np * ps(1);
  for i = 2:Np
    padj(i) = max (padj(i - 1), (Np - i + 1) * ps(i));
  endfor

  ## Reorder the adjusted p-values to match the order of the original p-values
  [jnk, original_order] = sort (idx, "ascend");
  padj = padj(original_order);

  ## Truncate adjusted p-values to 1.0
  padj(padj>1) = 1;

  ## If requested, calculate critical value at ALPHA
  ## No adjustment to confidence interval coverage
  if (nargout > 1)
    critval = tinv (1 - ALPHA / 2, dfe);
  endif

endfunction


function [padj, critval, dfe] = hochberg (p, t, Ng, dfe, R, ALPHA)

  ## Hochberg's step-up Bonferroni procedure

  ## Order raw p-values
  [ps, idx] = sort (p, "ascend");
  Np = numel (ps);

  ## Implement Hochberg's step-down Bonferroni procedure
  padj = nan (Np,1);
  padj(Np) = ps(Np);
  for j = 1:Np-1
    i = Np - j;
    padj(i) = min (padj(i + 1), (Np -i + 1) * ps(i));
  endfor

  ## Reorder the adjusted p-values to match the order of the original p-values
  [jnk, original_order] = sort (idx, "ascend");
  padj = padj(original_order);

  ## Truncate adjusted p-values to 1.0
  padj(padj>1) = 1;

  ## If requested, calculate critical value at ALPHA
  ## No adjustment to confidence interval coverage
  if (nargout > 1)
    critval = tinv (1 - ALPHA / 2, dfe);
  endif

endfunction


function [padj, critval, dfe] = fdr (p, t, Ng, dfe, R, ALPHA)

  ## Benjamini-Hochberg procedure to control the false discovery rate (FDR)
  ## This procedure does not control the family-wise error rate

  ## Order raw p-values
  [ps, idx] = sort (p, "ascend");
  Np = numel (ps);

  ## Initialize
  padj = nan (Np,1);
  alpha = nan (Np,1);

  ## Benjamini-Hochberg step-up procedure to control the false discovery rate
  padj = nan (Np,1);
  padj(Np) = ps(Np);
  for j = 1:Np-1
    i = Np - j;
    padj(i) = min (padj(i + 1), Np / i * ps(i));
  endfor

  ## Reorder the adjusted p-values to match the order of the original p-values
  [jnk, original_order] = sort (idx, "ascend");
  padj = padj(original_order);

  ## Truncate adjusted p-values to 1.0
  padj(padj>1) = 1;

  ## If requested, calculate critical value at ALPHA
  ## No adjustment to confidence interval coverage
  if (nargout > 1)
    critval = tinv (1 - ALPHA / 2, dfe);
  endif

endfunction


function [padj, critval, dfe] = lsd (p, t, Ng, dfe, R, ALPHA)

  ## Fisher's Least Significant Difference
  ## No control of the type I error rate across multiple comparisons
  padj = p;

  ## Calculate critical value at ALPHA
  ## No adjustment to confidence interval coverage
  critval = tinv (1 - ALPHA / 2, dfe);

endfunction


%!demo
%!
%! ## Demonstration using unbalanced one-way ANOVA example from anovan
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
%! ## Demonstration using factorial ANCOVA example from anovan
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

%!demo
%!
%! ## Demonstration using one-way ANOVA from anovan, with fit by weighted least
%! ## squares to account for heteroskedasticity.
%!
%! g = [1, 1, 1, 1, 1, 1, 1, 1, ...
%!      2, 2, 2, 2, 2, 2, 2, 2, ...
%!      3, 3, 3, 3, 3, 3, 3, 3]';
%!
%! y = [13, 16, 16,  7, 11,  5,  1,  9, ...
%!      10, 25, 66, 43, 47, 56,  6, 39, ...
%!      11, 39, 26, 35, 25, 14, 24, 17]';
%!
%! [P,ATAB,STATS] = anovan(y, g, "display", "off");
%! fitted = STATS.X * STATS.coeffs(:,1); # fitted values
%! b = polyfit (fitted, abs (STATS.resid), 1);
%! v = polyval (b, fitted);  # Variance as a function of the fitted values
%! [P,ATAB,STATS] = anovan (y, g, "weights", v.^-1, "display", "off");
%! [C, M] =  multcompare (STATS, "display", "on", "ctype", "mvt")

%!test
%!
%! ## Tests using unbalanced one-way ANOVA example from anovan and anova1
%!
%! ## Test for anovan - compare pairwise comparisons with matlab for CTYPE "lsd"
%!
%! dv =  [ 8.706 10.362 11.552  6.941 10.983 10.092  6.421 14.943 15.931 ...
%!        22.968 18.590 16.567 15.944 21.637 14.492 17.965 18.851 22.891 ...
%!        22.028 16.884 17.252 18.325 25.435 19.141 21.238 22.196 18.038 ...
%!        22.628 31.163 26.053 24.419 32.145 28.966 30.207 29.142 33.212 ...
%!        25.694 ]';
%! g = [1 1 1 1 1 1 1 1 2 2 2 2 2 3 3 3 3 3 3 3 3 ...
%!      4 4 4 4 4 4 4 5 5 5 5 5 5 5 5 5]';
%!
%! [P, ATAB, STATS] = anovan (dv, g, "varnames", "score", "display", "off");
%! [C, M, H, GNAMES] = multcompare (STATS, "dim", 1, "ctype", "lsd", ...
%!                                  "display", "off");
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
%!
%! ## Compare "fdr" adjusted p-values to those obtained using p.adjust in R
%!
%! [C, M, H, GNAMES] = multcompare (STATS, "dim", 1, "ctype", "fdr", ...
%!                                  "display", "off");
%! assert (C(1,6), 4.08303457454140e-05, 1e-09);
%! assert (C(2,6), 1.04587348240817e-06, 1e-09);
%! assert (C(3,6), 1.06397381604573e-07, 1e-09);
%! assert (C(4,6), 7.82091664406946e-14, 1e-09);
%! assert (C(5,6), 5.46591417210693e-01, 1e-09);
%! assert (C(6,6), 1.05737243156806e-01, 1e-09);
%! assert (C(7,6), 2.36859139493832e-07, 1e-09);
%! assert (C(8,6), 2.09859420867852e-01, 1e-09);
%! assert (C(9,6), 1.36324670121399e-07, 1e-09);
%! assert (C(10,6), 7.40712246958735e-06, 1e-09);
%!
%! ## Compare "hochberg" adjusted p-values to those obtained using p.adjust in R
%!
%! [C, M, H, GNAMES] = multcompare (STATS, "dim", 1, "ctype", "hochberg", ...
%!                                  "display", "off");
%! assert (C(1,6), 1.14324968087159e-04, 1e-09);
%! assert (C(2,6), 3.13762044722451e-06, 1e-09);
%! assert (C(3,6), 1.91515286888231e-07, 1e-09);
%! assert (C(4,6), 7.82091664406946e-14, 1e-09);
%! assert (C(5,6), 5.46591417210693e-01, 1e-09);
%! assert (C(6,6), 2.53769383576334e-01, 1e-09);
%! assert (C(7,6), 6.63205590582730e-07, 1e-09);
%! assert (C(8,6), 3.77746957562134e-01, 1e-09);
%! assert (C(9,6), 3.27179208291358e-07, 1e-09);
%! assert (C(10,6), 2.22213674087620e-05, 1e-09);
%!
%! ## Compare "holm" adjusted p-values to those obtained using p.adjust in R
%!
%! [C, M, H, GNAMES] = multcompare (STATS, "dim", 1, "ctype", "holm", ...
%!                                  "display", "off");
%! assert (C(1,6), 1.14324968087159e-04, 1e-09);
%! assert (C(2,6), 3.13762044722451e-06, 1e-09);
%! assert (C(3,6), 1.91515286888231e-07, 1e-09);
%! assert (C(4,6), 7.82091664406946e-14, 1e-09);
%! assert (C(5,6), 5.46591417210693e-01, 1e-09);
%! assert (C(6,6), 2.53769383576334e-01, 1e-09);
%! assert (C(7,6), 6.63205590582730e-07, 1e-09);
%! assert (C(8,6), 3.77746957562134e-01, 1e-09);
%! assert (C(9,6), 3.27179208291358e-07, 1e-09);
%! assert (C(10,6), 2.22213674087620e-05, 1e-09);
%!
%! ## Compare "scheffe" adjusted p-values to those obtained using 'scheffe' in Matlab
%!
%! [C, M, H, GNAMES] = multcompare (STATS, "dim", 1, "ctype", "scheffe", ...
%!                                  "display", "off");
%! assert (C(1,6), 0.00108105386141085, 1e-09);
%! assert (C(2,6), 2.7779386789517e-05, 1e-09);
%! assert (C(3,6), 1.3599854038198e-06, 1e-09);
%! assert (C(4,6), 7.58830197867751e-13, 1e-09);
%! assert (C(5,6), 0.984039948220281, 1e-09);
%! assert (C(6,6), 0.539077018557706, 1e-09);
%! assert (C(7,6), 5.59475764460574e-06, 1e-09);
%! assert (C(8,6), 0.771173490574105, 1e-09);
%! assert (C(9,6), 2.52838425729905e-06, 1e-09);
%! assert (C(10,6), 0.000200719143889168, 1e-09);
%!
%! ## Compare "bonferroni" adjusted p-values to those obtained using p.adjust in R
%!
%! [C, M, H, GNAMES] = multcompare (STATS, "dim", 1, "ctype", "bonferroni", ...
%!                                  "display", "off");
%! assert (C(1,6), 2.85812420217898e-04, 1e-09);
%! assert (C(2,6), 5.22936741204085e-06, 1e-09);
%! assert (C(3,6), 2.12794763209146e-07, 1e-09);
%! assert (C(4,6), 7.82091664406946e-14, 1e-09);
%! assert (C(5,6), 1.00000000000000e+00, 1e-09);
%! assert (C(6,6), 8.45897945254446e-01, 1e-09);
%! assert (C(7,6), 9.47436557975328e-07, 1e-09);
%! assert (C(8,6), 1.00000000000000e+00, 1e-09);
%! assert (C(9,6), 4.08974010364197e-07, 1e-09);
%! assert (C(10,6), 4.44427348175241e-05, 1e-09);
%!
%! ## Test for anova1 ("equal")- comparison of results from Matlab
%!
%! [P, ATAB, STATS] = anova1 (dv, g, "off", "equal");
%! [C, M, H, GNAMES] = multcompare (STATS, "ctype", "lsd", "display", "off");
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
%!
%! ## Test for anova1 ("unequal") - comparison with results from GraphPad Prism 8
%! [P, ATAB, STATS] = anova1 (dv, g, "off", "unequal");
%! [C, M, H, GNAMES] = multcompare (STATS, "ctype", "lsd", "display", "off");
%! assert (C(1,6), 0.001247025266382, 1e-09);
%! assert (C(2,6), 0.000018037115146, 1e-09);
%! assert (C(3,6), 0.000002974595187, 1e-09);
%! assert (C(4,6), 0.000000000786046, 1e-09);
%! assert (C(5,6), 0.5693192886650109, 1e-09);
%! assert (C(6,6), 0.110501699029776, 1e-09);
%! assert (C(7,6), 0.000131226488700, 1e-09);
%! assert (C(8,6), 0.1912101409715992, 1e-09);
%! assert (C(9,6), 0.000005385256394, 1e-09);
%! assert (C(10,6), 0.000074089106171, 1e-09);

%!test
%!
%! ## Test for anova2 ("interaction") - comparison with results from Matlab for column effect
%! popcorn = [5.5, 4.5, 3.5; 5.5, 4.5, 4.0; 6.0, 4.0, 3.0; ...
%!            6.5, 5.0, 4.0; 7.0, 5.5, 5.0; 7.0, 5.0, 4.5];
%! [P, ATAB, STATS] = anova2 (popcorn, 3, "off");
%! [C, M, H, GNAMES] = multcompare (STATS, "estimate", "column",...
%!                                  "ctype", "lsd", "display", "off");
%! assert (C(1,6), 1.49311100811177e-05, 1e-09);
%! assert (C(2,6), 2.20506904243535e-07, 1e-09);
%! assert (C(3,6), 0.00449897860490058, 1e-09);
%! assert (M(1,1), 6.25, 1e-09);
%! assert (M(2,1), 4.75, 1e-09);
%! assert (M(3,1), 4, 1e-09);
%! assert (M(1,2), 0.152145154862547, 1e-09);
%! assert (M(2,2), 0.152145154862547, 1e-09);
%! assert (M(3,2), 0.152145154862547, 1e-09);

%!test
%!
%! ## Test for anova2 ("linear") - comparison with results from GraphPad Prism 8
%! words = [10 13 13; 6 8 8; 11 14 14; 22 23 25; 16 18 20; ...
%!          15 17 17; 1 1 4; 12 15 17;  9 12 12;  8 9 12];
%! [P, ATAB, STATS] = anova2 (words, 1, "off", "linear");
%! [C, M, H, GNAMES] = multcompare (STATS, "estimate", "column",...
%!                                  "ctype", "lsd", "display", "off");
%! assert (C(1,6), 0.000020799832702, 1e-09);
%! assert (C(2,6), 0.000000035812410, 1e-09);
%! assert (C(3,6), 0.003038942449215, 1e-09);

%!test
%!
%! ## Test for anova2 ("nested") - comparison with results from GraphPad Prism 8
%! data = [4.5924 7.3809 21.322; -0.5488 9.2085 25.0426; ...
%!         6.1605 13.1147 22.66; 2.3374 15.2654 24.1283; ...
%!         5.1873 12.4188 16.5927; 3.3579 14.3951 10.2129; ...
%!         6.3092 8.5986 9.8934; 3.2831 3.4945 10.0203];
%! [P, ATAB, STATS] = anova2 (data, 4, "off", "nested");
%! [C, M, H, GNAMES] = multcompare (STATS, "estimate", "column",...
%!                                  "ctype", "lsd", "display", "off");
%! assert (C(1,6), 0.261031111511073, 1e-09);
%! assert (C(2,6), 0.065879755907745, 1e-09);
%! assert (C(3,6), 0.241874613529270, 1e-09);

%!shared visibility_setting
%! visibility_setting = get (0, "DefaultFigureVisible");

%!test
%! set (0, "DefaultFigureVisible", "off");
%!
%! ## Test for kruskalwallis - comparison with results from MATLAB
%! data = [3,2,4; 5,4,4; 4,2,4; 4,2,4; 4,1,5; ...
%!         4,2,3; 4,3,5; 4,2,4; 5,2,4; 5,3,3];
%! group = [1:3] .* ones (10,3);
%! [P, ATAB, STATS] = kruskalwallis (data(:), group(:), "off");
%! C = multcompare (STATS, "ctype", "lsd");
%! assert (C(1,6), 0.000163089828959986, 1e-09);
%! assert (C(2,6), 0.630298044801257, 1e-09);
%! assert (C(3,6), 0.00100567660695682, 1e-09);
%! C = multcompare (STATS, "ctype", "bonferroni", "display", "off");
%! assert (C(1,6), 0.000489269486879958, 1e-09);
%! assert (C(2,6), 1, 1e-09);
%! assert (C(3,6), 0.00301702982087047, 1e-09);
%! C = multcompare(STATS, "ctype", "scheffe");
%! assert (C(1,6), 0.000819054880289573, 1e-09);
%! assert (C(2,6), 0.890628039849261, 1e-09);
%! assert (C(3,6), 0.00447816059021654, 1e-09);
%! set (0, "DefaultFigureVisible", visibility_setting);

%!test
%! set (0, "DefaultFigureVisible", "off");
%! ## Test for friedman - comparison with results from MATLAB
%! popcorn = [5.5, 4.5, 3.5; 5.5, 4.5, 4.0; 6.0, 4.0, 3.0; ...
%!            6.5, 5.0, 4.0; 7.0, 5.5, 5.0; 7.0, 5.0, 4.5];
%! [P, ATAB, STATS] = friedman (popcorn, 3, "off");
%! C = multcompare(STATS, "ctype", "lsd");
%! assert (C(1,6), 0.227424558028569, 1e-09);
%! assert (C(2,6), 0.0327204848315735, 1e-09);
%! assert (C(3,6), 0.353160353315988, 1e-09);
%! C = multcompare(STATS, "ctype", "bonferroni");
%! assert (C(1,6), 0.682273674085708, 1e-09);
%! assert (C(2,6), 0.0981614544947206, 1e-09);
%! assert (C(3,6), 1, 1e-09);
%! C = multcompare(STATS, "ctype", "scheffe", "display", "off");
%! assert (C(1,6), 0.482657360384373, 1e-09);
%! assert (C(2,6), 0.102266573027672, 1e-09);
%! assert (C(3,6), 0.649836502233148, 1e-09);
%! set (0, "DefaultFigureVisible", visibility_setting);

%!test
%! set (0, "DefaultFigureVisible", "off");
%! ## Test for fitlm - same comparisons as for first anovan example
%! y =  [ 8.706 10.362 11.552  6.941 10.983 10.092  6.421 14.943 15.931 ...
%!        22.968 18.590 16.567 15.944 21.637 14.492 17.965 18.851 22.891 ...
%!        22.028 16.884 17.252 18.325 25.435 19.141 21.238 22.196 18.038 ...
%!        22.628 31.163 26.053 24.419 32.145 28.966 30.207 29.142 33.212 ...
%!        25.694 ]';
%! X = [1 1 1 1 1 1 1 1 2 2 2 2 2 3 3 3 3 3 3 3 3 4 4 4 4 4 4 4 5 5 5 5 5 5 5 5 5]';
%! [TAB,STATS] = fitlm (X,y,"linear","categorical",1,"display","off");
%! [C, M] = multcompare(STATS, "ctype", "lsd");
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
%! set (0, "DefaultFigureVisible", visibility_setting);
