## Copyright (C) 2022 Andrew Penn <A.C.Penn@sussex.ac.uk>
## Copyright (C) 2022 Andreas Bertsatos <abertsatos@biol.uoa.gr>
## Copyright (C) 2021 Christian Scholz
## Copyright (C) 2003-2005 Andy Adler <adler@ncf.ca>
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
## @deftypefn {Function File} [@var{P}, @var{T}, @var{STATS}, @var{TERMS}] = anovan (@var{Y}, @var{GROUP})
## @deftypefnx {Function File} [@var{P}, @var{T}, @var{STATS}, @var{TERMS}] = anovan (@var{Y}, @var{GROUP}, "name", @var{value})
##
## Perform a multi (N)-way analysis of (co)variance (ANOVA or ANCOVA) to
## evaluate the effect of one or more categorical or continuous predictors
## on a continuous outcome. The algorithms used make @qcode{anovan} suitable
## for balanced or unbalanced factorial (crossed) designs. By default, @qcode{anovan}
## treats all factors as fixed. Examples of function usage can be found by
## entering the command @code{demo anovan}.
##
## Data is a single vector @var{Y} with groups specified by a corresponding
## matrix or cell array of group labels @var{GROUP}, where each column of
## @var{GROUP} has the same number of rows as @var{Y}. For example, if @var{Y}
## = [1.1;1.2]; @var{GROUP} = [1,2,1; 1,5,2]; then observation 1.1 was measured
## under conditions 1,2,1 and observation 1.2 was measured under conditions
## 1,5,2. Note that @var{GROUP} do not need to be sequentially numbered.
##
## @qcode{anovan} can take a number of optional parameters as name-value pairs.
##
## @code{[@dots{}] = anovan (@var{Y}, @var{GROUP}, "continuous", @var{continuous})}
##
## @itemize
## @item
## @var{continuous} is a vector of indices indicating which of the columns (i.e.
## factors) in @var{GROUP} should be treated as continuous predictors rather
## than as categorical predictors. The relationship between continuous predictors 
## and the outcome should be linear.
## @end itemize
##
## @code{[@dots{}] = anovan (@var{Y}, @var{GROUP}, "random", @var{random})}
##
## @itemize
## @item
## @var{random} is a vector of indices indicating which of the columns (i.e.
## factors) in @var{GROUP} should be treated as random effects rather than fixed
## effects. Octave @qcode{anovan} provides only basic support for random effects.
## Specifically, since all F-statistics in @qcode{anovan} are calculated using
## the mean-squared error (MSE), any interaction terms containing a random effect
## are dropped from the model term definitions and their associated variance 
## is pooled with the residual, unexplained variance making up the MSE. Variable 
## names for random factors are appended with a ' symbol.
## @end itemize
##
## @code{[@dots{}] = anovan (@var{Y}, @var{GROUP}, "model", @var{modeltype})}
##
## @var{modeltype} can specified as one of the following:
##
## @itemize
## @item
## "linear" (default) : compute N main effects with no interactions.
##
## @item
## "interaction" : compute N effects and N*(N-1) two-factor interactions
##
## @item
## "full" : compute the N main effects and interactions at all levels
##
## @item
## a scalar integer : representing the maximum interaction order
##
## @item
## a matrix of term definitions : each row is a term and each column is a factor
## @end itemize
##
## @example
## -- Example:
## A two-way ANOVA with interaction would be: [1 0; 0 1; 1 1]
## @end example
##
## @code{[@dots{}] = anovan (@var{Y}, @var{GROUP}, "sstype", @var{sstype})}
##
## @var{sstype} can specified as one of the following:
##
## @itemize
## @item
## 1 : Type I sequential sums-of-squares.
##
## @item
## 2 or 'h': Type II partially sequential (or hierarchical) sums-of-squares
##
## @item
## 3 (default): Type III partial, constrained or marginal sums-of-squares
##
## @end itemize
##
## @code{[@dots{}] = anovan (@var{Y}, @var{GROUP}, "varnames", @var{varnames})}
##
## @itemize
## @item
## @var{varnames} must be a cell array of strings with each element containing a
## factor name for each column of GROUP.  By default (if not parsed as optional
## argument), @var{varnames} are "X1","X2","X3", etc.
## @end itemize
##
## @code{[@dots{}] = anovan (@var{Y}, @var{GROUP}, "display", @var{dispopt})}
##
## @itemize
## @item
## @var{dispopt} can be either "on" (default) or "off" and switches the display
## of the ANOVA table.
## @end itemize
##
## @code{[@dots{}] = anovan (@var{Y}, @var{GROUP}, "contrasts", @var{contrasts})}
##
## @itemize
## @item
## @var{contrasts} must be a cell array, where each cell contains a matrix
## of the contrast coding scheme (i.e. the generalized inverse of contrast
## weights) for the corresponding column (i.e. factor) in @var{GROUP} (unless
## @var{GROUP} contains only one factor in which case the contrast matrix need
## not be provided within in a cell array). Rows in the contrast matrices
## correspond to factor levels in the order that they first appear in the
## corresponding @var{GROUP} column (and as they appear in @var{STATS}.grpnames).
## If cells are left empty, then the default contrasts are applied, specifically
## deviation effect coding. Cells corresponding to continuous factors are ignored
## and can be left empty. Note that SSTYPE 3 cannot be calculated if columns in
## the contrast matrices do not sum to zero (at single precision), in which case
## @qcode{anovan} will fall back to SSTYPE 2.
## @end itemize
##
## @qcode{anovan} can return up to four output arguments:
##
## @var{P} = anovan (@dots{}) returns a vector of p-values, one for each term.
##
## [@var{P}, @var{T}] = anovan (@dots{}) returns a cell array containing the
## ANOVA table.
##
## [@var{P}, @var{T}, @var{STATS}] = anovan (@dots{}) returns a structure
## containing additional statistics, including degrees of freedom and effect
## sizes for each term in the linear model, the design matrix, the variance-covariance
## matrix, model residuals, and the mean squared error. The columns of @var{STATS}.coeffs
## (from left-to-right) report the model coefficients, standard errors, t-statistics
## and p-values, which relate to the contrasts. The number appended to each term
## name in @var{STATS}.coeffnames relates to the column number in the corresponding
## contrast matrix for that factor.
##
## [@var{P}, @var{T}, @var{STATS}, @var{TERMS}] = anovan (@dots{}) returns the
## model term definitions.
##
##
## @end deftypefn

##  Author: Andrew Penn <a.c.penn@sussex.ac.uk>
##  Includes some code by: Andy Adler <adler@ncf.ca>, Christian Scholz and Andreas Bertsatos
##

function [P, T, STATS, TERMS] = anovan (Y, GROUP, varargin)

    if (nargin <= 1)
      error (strcat (["anovan usage: ""anovan (Y, GROUP)""; "], ...
                      [" atleast 2 input arguments required"]));
    endif

    ## Check supplied parameters
    if ((numel (varargin) / 2) != fix (numel (varargin) / 2))
      error ("anovan: wrong number of arguments.")
    endif
    MODELTYPE = "linear";
    DISPLAY = "on";
    SSTYPE = 3;
    VARNAMES = [];
    CONTINUOUS = [];
    RANDOM = [];
    CONTRASTS = {};
    for idx = 3:2:nargin
      name = varargin{idx-2};
      value = varargin{idx-1};
      switch (lower (name))
        case "model"
          MODELTYPE = value;
        case "continuous"
          CONTINUOUS = value;
        case "random"
          RANDOM = value;
        case "nested"
          error (strcat (["anovan: nested ANOVA is not supported. Please use"], ... 
                         [" anova2 for fully balanced nested ANOVA designs."]));
        case "sstype"
          SSTYPE = value;
        case "varnames"
          VARNAMES = value;
        case "display"
          DISPLAY = value;
        case "contrasts"
          CONTRASTS = value;
        otherwise
          error (sprintf ("anovan: parameter %s is not supported", name));
      endswitch
    endfor
    if (isnumeric (CONTINUOUS))
      if (any (CONTINUOUS != abs (fix (CONTINUOUS))))
        error (strcat (["anovan: the value provided for the CONTINUOUS"], ...
                       [" parameter must be a positive integer"]));
      endif
    else
      error (strcat (["anovan: the value provided for the CONTINUOUS"], ...
                     [" parameter must be numeric"]));
    endif

    ## Accomodate for different formats for GROUP
    ## GROUP can be a matrix of numeric identifiers of a cell arrays
    ## of strings or numeric idenitiers
    N = size (GROUP,2); # number of anova "ways"
    n = numel (Y);      # total number of observations
    if (prod (size (Y)) != n)
      error ("anovan: for ""anovan (Y, GROUP)"", Y must be a vector");
    endif
    if (numel (unique (CONTINUOUS)) > N)
      error (strcat (["anovan: the number of factors assigned as"], ...
                 [" continuous cannot exceed the number of factors in GROUP"]));
    endif
    if (any ((CONTINUOUS > N) || any (CONTINUOUS <= 0)))
      error (strcat (["anovan: one or more indices provided in the"], ...
                 [" value for the continuous parameter are out of range"]));
    endif
    cont_vec = false (1, N);
    cont_vec(CONTINUOUS) = true;
    if (iscell (GROUP))
      if (size(GROUP, 1) == 1)
        for j = 1:N
          if (isnumeric (GROUP{j}))
            if (any (CONTINUOUS == j))
              tmp(:,j) = num2cell (GROUP{j});
            else
              tmp(:,j) = cellstr (num2str (GROUP{j}));
            endif
          else
            if (any (CONTINUOUS == j))
              error ("anovan: continuous factors must be a numeric datatype");
            endif
            tmp(:,j) = GROUP{j};
          endif
        endfor
        GROUP = tmp;
      else
        for j = 1:N
          tmp(:,j) = cellstr (char (GROUP{j}));
        endfor
      endif
    endif
    if (size (GROUP,1) != n)
      error ("anovan: GROUP must be a matrix of the same number of rows as Y");
    endif
    if (! isempty (VARNAMES))
      if (iscell (VARNAMES))
        if (all (cellfun (@ischar, VARNAMES)))
          nvarnames = numel(VARNAMES);
        else
          error (strcat (["anovan: all variable names must be character"], ...
                         [" or character arrays"]));
        endif
      elseif (ischar (VARNAMES))
        nvarnames = 1;
        VARNAMES = {VARNAMES};
      elseif (isstring (VARNAMES))
        nvarnames = 1;
        VARNAMES = {char(VARNAMES)};
      else
        error (strcat (["anovan: varnames is not of a valid type. Must be"], ...
               [" cell array of character arrays, character array or string"]));
      endif
    else
      nvarnames = N;
      VARNAMES = arrayfun(@(x) ["X",num2str(x)], 1:N, "UniformOutput", 0);
    endif
    if (nvarnames != N)
      error (strcat (["anovan: number of variable names is not equal"], ...
                     ["  to number of grouping variables"]));
    endif

    ## Evaluate random argument (if applicable)
    if (! isempty(RANDOM))
      if (isnumeric (RANDOM))
        if (any (RANDOM != abs (fix (RANDOM))))
          error (strcat (["anovan: the value provided for the RANDOM"], ...
                         [" parameter must be a positive integer"]));
        endif
      else
        error (strcat (["anovan: the value provided for the RANDOM"], ...
                     ["  parameter must be numeric"]));
      endif
      if (numel (RANDOM) > N)
        error (strcat (["anovan: the number of elements in RANDOM cannot"], ...
               [" exceed the number of columns in GROUP."]));
      endif
      if (max (RANDOM) > N)
        error (strcat (["anovan: the indices listed in RANDOM cannot"], ...
               [" exceed the number of columns in GROUP."]));
      endif
      for v = 1:N
        if (ismember (v, RANDOM))
          VARNAMES{v} = strcat (VARNAMES{v},"'");
        endif
      endfor
    endif

    ## Evaluate contrasts (if applicable)
    if isempty (CONTRASTS)
      CONTRASTS = cell (1, N);
    else
      if (! iscell (CONTRASTS))
        CONTRASTS = {CONTRASTS};
      endif
      if (numel (CONTRASTS) != N)
        error (strcat (["anovan: the number of contrast matrices in"], ...
                 ["  CONTRASTS does not match the number of main effects"]));
      endif
      for i = 1:N
        if (! isempty (CONTRASTS{i}))
          ## Check that the columns sum to 0
          if (any (abs (sum (CONTRASTS{i})) > eps("single")) ...
                        && strcmpi (num2str (SSTYPE), "3"))
            ## Contrasts must sum to 0 for SSTYPE 3.
            ## If they don't, use SSTYPE 2 instead
            SSTYPE = 2;
          endif
        endif
      endfor
    endif

    ## Remove NaN or non-finite observations
    excl = logical (isnan(Y) + isinf(Y));
    Y(excl) = [];
    GROUP(excl,:) = [];
    if (size (Y, 1) == 1)
      Y = Y.';         # if Y is a row vector, make it a column vector
    endif
    n = numel (Y);     # recalculate total number of observations

    ## Evaluate model type input argument and create terms matrix if not provided
    msg = strcat (["anovan: the number of columns in the term definitions"], ...
                  [" cannot exceed the number of columns of GROUP"]);
    if (ischar (MODELTYPE))
      switch (lower (MODELTYPE))
        case "linear"
          MODELTYPE = 1;
        case "interaction"
          MODELTYPE = 2;
        case "full"
          MODELTYPE = N;
        otherwise
          error ("anovan: model type not recognised");
      endswitch
    endif
    if (isscalar (MODELTYPE))
      TERMS = cell (MODELTYPE,1);
      v = false (1, N);
      switch (lower (MODELTYPE))
        case 1
          ## Create term definitions for an additive linear model
          TERMS = eye (N);
        case 2
          ## Create term definitions for a model with two factor interactions
          Nx = nchoosek (N, 2);
          TERMS = zeros (N + Nx, N);
          TERMS(1:N,:) = eye (N);
          for j = 1:N
            for i = j:N-1
              TERMS(N+j+i-1,j) = 1;
              TERMS(N+j+i-1,i+1) = 1;
            endfor
          endfor
        otherwise
          if (MODELTYPE > N)
            error (msg);
          endif
          ## Create term definitions for a full model
          Nx = zeros (1, N-1);
          Nx = 0;
          for k = 1:N
            Nx = Nx + nchoosek(N,k);
          endfor
          for j = 1:MODELTYPE
            v(1:j) = 1;
            TERMS(j) = flipud (unique (perms (v), "rows"));
          endfor
          TERMS = cell2mat (TERMS);
      endswitch
      TERMS = logical (TERMS);
    else
      ## Assume that the user provided a suitable matrix of term definitions
      if (size (MODELTYPE, 2) > N)
        error (msg);
      endif
      if (! all (ismember (MODELTYPE(:), [0,1])))
        error (strcat (["anovan: elements of the model terms matrix"], ...
                       [" must be either 0 or 1"]));
      endif
      TERMS = logical (MODELTYPE);
    endif
    ## Evaluate terms matrix
    Ng = sum (TERMS, 2);
    if (any (diff (Ng) < 0))
      error (strcat (["anovan: the model terms matrix must list main"], ...
                     [" effects above/before interactions"]));
    endif
    Nm = sum (Ng == 1);
    ## Drop terms that include interactions with random effects.
    drop = any (bsxfun (@and, TERMS(:,RANDOM), (Ng > 1)), 2);
    TERMS (drop, :) = [];
    Ng(drop) = [];
    Nx = sum (Ng > 1);
    Nt = Nm + Nx;
    if (any (any (TERMS(1:Nm,:), 1) != any (TERMS, 1)))
      error (strcat (["anovan: all factors involved in interactions"], ...
                     [" must have a main effect"]));
    endif

    ## Calculate total sum-of-squares
    ct  = sum (Y)^2 / n;   % correction term
    sst = sum (Y.^2) - ct;
    dft = n - 1;

    ## Fit linear models, and calculate sums-of-squares for ANOVA
    switch (lower (SSTYPE))
      case 1
        ## Type I sequential sums-of-squares (SSTYPE = 1)
        R = sst;
        ss = zeros (Nt,1);
        [X, grpnames, nlevels, df, termcols, coeffnames, vmeans, gid, ...
         CONTRASTS] = make_design_matrix (GROUP, TERMS, CONTINUOUS, ...
                      CONTRASTS, VARNAMES, n, Nm, Nx, Ng);
        for j = 1:Nt
          XS = cell2mat (X(1:j+1));
          [b, sse] = lmfit (XS, Y);
          ss(j) = R - sse;
          R = sse;
        endfor
        [b, sse, resid, ucov] = lmfit (XS, Y);
        sstype_char = "I";
      case {2,'h'}
        ## Type II (partially sequential, or hierarchical) sums-of-squares
        ss = zeros (Nt,1);
        [X, grpnames, nlevels, df, termcols, coeffnames, vmeans, gid, ...
         CONTRASTS] = make_design_matrix (GROUP, TERMS, CONTINUOUS, ...
                      CONTRASTS, VARNAMES, n, Nm, Nx, Ng);
        for j = 1:Nt
          i = find (TERMS(j,:));
          k = cat (1, 1, 1 + find (any (!TERMS(:,i),2)));
          XS = cell2mat (X(k));
          [jnk, R1] = lmfit (XS, Y);
          k = cat (1, j+1, k);
          XS = cell2mat (X(k));
          [jnk, R2] = lmfit (XS, Y);
          ss(j) = R1 - R2;
        endfor
        [b, sse, resid, ucov] = lmfit (cell2mat (X), Y);
        sstype_char = "II";
      case 3
        ## Type III (partial, constrained or marginal) sums-of-squares
        ss = zeros (Nt, 1);
        [X, grpnames, nlevels, df, termcols, coeffnames, vmeans, gid, ...
         CONTRASTS] = make_design_matrix (GROUP, TERMS, CONTINUOUS, ...
                      CONTRASTS, VARNAMES, n, Nm, Nx, Ng);
        [b, sse, resid, ucov] = lmfit (cell2mat (X), Y);
        for j = 1:Nt
          XS = cell2mat (X(1:Nt+1 != j+1));
          [jnk, R, resid] = lmfit (XS, Y);
          ss(j) = R - sse;
        endfor
        sstype_char = "III";
      otherwise
        error ("anovan: sstype value not supported");
    endswitch
    dfe = dft - sum (df);
    ms = ss ./ df;
    mse = sse / dfe;
    eta_sq = ss ./ sst;
    partial_eta_sq = ss ./ (ss + sse);
    F = ms / mse;
    P = 1 - fcdf (F, df, dfe);

    ## Prepare cell array containing the ANOVA table (T)
    T = cell (Nt + 3, 7);
    T(1,:) = {"Source", "Sum Sq.", "d.f.", "Mean Sq.", "Eta Sq.", "F", "Prob>F"};
    T(2:Nt+1,2:7) = num2cell ([ss df ms partial_eta_sq F P]);
    T(end-1,1:4) = {"Error", sse, dfe, mse};
    T(end,1:3) = {"Total", sst, dft};
    for i = 1:Nt
      str = sprintf ("%s*", VARNAMES{find (TERMS(i,:))});
      T(i+1,1) = str(1:end-1);
    endfor

    ## If requested, prepare stats output structure
    ## Note that the information provided by STATS is not sufficient
    ## for multcompare function as yet
    if (nargout > 2)

      ## Calculate a standard error, t-statistic and p-value for each
      ## of the regression coefficients
      se = sqrt (diag (ucov) * mse);
      t =  b ./ se;
      coeff_stats = zeros (1 + sum (df), 4);
      coeff_stats(:,1) = b;                                # coefficients
      coeff_stats(:,2) = se;                               # standard errors
      coeff_stats(:,3) = t;                                # t-statistics
      coeff_stats(:,4) = 2 * (1 - (tcdf (abs (t), dfe)));  # p-values

      STATS = struct ("source","anovan", ...
                      "resid", resid, ...
                      "coeffs", coeff_stats, ...
                      "Rtr", [], ...           # Not used by Octave
                      "rowbasis", [], ...      # Not used by Octave
                      "dfe", dfe, ...
                      "mse", mse, ...
                      "nullproject", [], ...   # Not used by Octave
                      "terms", TERMS, ...
                      "nlevels", nlevels, ...
                      "continuous", cont_vec, ...
                      "vmeans", vmeans, ...
                      "termcols", termcols, ...
                      "coeffnames", {vertcat(coeffnames{:})}, ...
                      "vars", [], ...          # Not used by Octave
                      "varnames", {VARNAMES}, ...
                      "grpnames", {grpnames}, ...
                      "vnested", [], ...       # Not used since "nested" argument name is not supported
                      "ems", [], ...           # Not used since "nested" argument name is not supported
                      "denom", [], ...         # Not used since interactions with random effects is not supported
                      "dfdenom", [], ...       # Not used since interactions with random effects is not supported
                      "msdenom", [], ...       # Not used since interactions with random effects is not supported
                      "varest", [], ...        # Not used since interactions with random effects is not supported
                      "varci", [], ...         # Not used since interactions with random effects is not supported
                      "txtdenom", [], ...      # Not used since interactions with random effects is not supported
                      "txtems", [], ...        # Not used since interactions with random effects is not supported
                      "rtnames", [], ...       # Not used since interactions with random effects is not supported
                      ## Additional STATS fields used exclusively by Octave
                      "df", df, ...
                      "contrasts", {CONTRASTS}, ...
                      "X", sparse (cell2mat (X)), ...
                      "vcov", sparse (ucov * mse), ...
                      "grps", gid, ...
                      "eta_squared", eta_sq, ...
                      "partial_eta_squared", partial_eta_sq);
    endif

    ## Print ANOVA table
    switch (lower (DISPLAY))
      case "on"
        ## Get dimensions of the ANOVA table
        [nrows, ncols] = size (T);
        ## Print table
        fprintf("\nANOVA table (Type %s sums-of-squares):\n\n", sstype_char);
        fprintf("Source                   Sum Sq.    d.f.    Mean Sq.  R Sq.            F  Prob>F\n");
        fprintf("********************************************************************************\n");
        for i = 1:Nt
          str = T{i+1,1};
          l = numel(str);  # Needed to truncate source term name at 18 characters
          ## Format and print the statistics for each model term
          ## Format F statistics and p-values in APA style
          if (P(i) < 0.001)
            fprintf ("%-20s  %10.5g  %6d  %10.5g  %4.3f  %11.2f   <.001 \n", ...
                      str(1:min(18,l)), T{i+1,2:end-1});
          elseif (P(i) < 1.0)
            fprintf ("%-20s  %10.5g  %6d  %10.5g  %4.3f  %11.2f   .%03u \n", ...
                      str(1:min(18,l)), T{i+1,2:end-1}, round (P(i) * 1e+03));
          else
            fprintf ("%-20s  %10.5g  %6d  %10.5g  %4.3f  %11.2f   1.000 \n", ...
                      str(1:min(18,l)), T{i+1,2:end-1});
          endif
        endfor
        fprintf("Error                 %10.5g  %6d  %10.5g\n", T{end-1,2:4});
        fprintf("Total                 %10.5g  %6d \n", T{end,2:3});
        fprintf("\n");
      case "off"
        ## do nothing
      otherwise
        error ("anovan: wrong value for ""display"" parameter.");
    endswitch

endfunction


function [X, levels, nlevels, df, termcols, coeffnames, vmeans, gid, ...
          CONTRASTS] = make_design_matrix (GROUP, TERMS, CONTINUOUS, ...
                       CONTRASTS, VARNAMES, n, Nm, Nx, Ng)

  ## Returns a cell array of the design matrix for each term in the model

  ## Fetch factor levels from each column (i.e. factor) in GROUP
  levels = cell (Nm, 1);
  gid = zeros (n, Nm);
  nlevels = zeros (Nm, 1);
  df = zeros (Nm + Nx, 1);
  termcols = ones (1 + Nm + Nx, 1);
  for j = 1:Nm
    if (any (j == CONTINUOUS))
      ## Continuous predictor
      nlevels(j) = 1;
      termcols(j+1) = 1;
      df(j) = 1;
      if iscell (GROUP(:,j))
        gid(:,j) = cell2mat ([GROUP(:,j)]);
      else
        gid(:,j) = GROUP(:,j);
      end
    else
      ## Categorical predictor
      levels{j} = unique (GROUP(:,j), "stable");
      if isnumeric (levels{j})
        levels{j} = num2cell (levels{j});
      endif
      nlevels(j) = numel (levels{j});
      for k = 1:nlevels(j)
        gid(ismember (GROUP(:,j),levels{j}{k}),j) = k;
      endfor
      termcols(j+1) = nlevels(j);
      df(j) = nlevels(j) - 1;
    endif
  endfor

  ## Create design matrix
  ## Prepare design matrix columns for the main effects
  X = cell (1, 1 + Nm + Nx);
  X(1) = ones (n, 1);
  coeffnames = cell (1, 1 + Nm + Nx);
  coeffnames(1) = "(Intercept)";
  vmeans = zeros (Nm, 1);
  for j = 1:Nm
    if (any (j == CONTINUOUS))
      ## Continuous predictor
      if iscell (GROUP(:,j))
        X(1+j) = cell2mat (GROUP(:,j));
      else
        X(1+j) = GROUP(:,j);
      end
      vmeans(j) = mean ([X{1+j}]);
      #X(1+j) = [X{1+j}] - vmeans(j); # Like Matlab, we won't center continuous factors
      if (isempty (CONTRASTS{j}))
        CONTRASTS{j} = [];
      end
      ## Create names of the coefficients relating to continuous main effects
      coeffnames{1+j} = VARNAMES{j};
    else
      ## Categorical predictor
      if (isempty (CONTRASTS{j}))
        CONTRASTS{j} = contr_sum (nlevels(j));
      else
        ## Check that the contrast matrix provided is the correct size
        if (! all (size (CONTRASTS{j},1) == nlevels(j)))
          error (strcat (["anovan: the number of rows in the contrast"], ...
                       [" matrices should equal the number of factor levels"]));
        endif
        if (! all (size (CONTRASTS{j},2) == df(j)))
          error (strcat (["anovan: the number of columns in each contrast"], ...
                  [" matrix should equal the degrees of freedom (i.e."], ...
                  [" number of levels minus 1) for that factor"]));
        endif
        if (! all (any (CONTRASTS{j})))
          error (strcat (["anovan: a contrast must be coded in each"], ...
                         [" column of the contrast matrices"]));
        endif
      endif
      C = CONTRASTS{j};
      func = @(x) x(gid(:,j));
      X(1+j) = cell2mat (cellfun (func, num2cell (C, 1), "UniformOutput", false));
      ## Create names of the coefficients relating to continuous main effects
      coeffnames{1+j} = cell (df(j),1);
      for v = 1:df(j)
        coeffnames{1+j}{v} = sprintf ("%s_%u", VARNAMES{j}, v);
      endfor
    endif
  endfor
  ## If applicable, prepare design matrix columns for all the interaction terms
  if (Nx > 0)
    row = TERMS((Ng > 1),:);
    for i = 1:Nx
      I = 1 + find (row(i,:));
      df(Nm+i) = prod (df(I-1));
      termcols(1+Nm+i) = prod (df(I-1) + 1);
      tmp = ones (n,1);
      for j = 1:numel(I);
        tmp = num2cell (tmp, 1);
        for k = 1:numel(tmp)
          tmp(k) = bsxfun (@times, tmp{k}, X{I(j)});
        endfor
        tmp = cell2mat(tmp);
      endfor
      X{1+Nm+i} = tmp;
      coeffnames{1+Nm+i} = cell (df(Nm+i),1);
      for v = 1:df(Nm+i)
        coeffnames{1+Nm+i}{v} = strcat (sprintf ("%s_", VARNAMES{I-1}), num2str (v));
      endfor
    endfor
  endif

  function C = contr_sum (N)

    ## Create contrast matrix (of doubles) using deviation effect coding
    ## These contrasts are centered (i.e. sum to 0)
    C =  cat (1, eye (N-1), - (ones (1,N-1)));

  endfunction

endfunction


function [b, sse, resid, ucov] = lmfit (X, Y)

  ## Get model coefficients by solving the linear equation by QR decomposition
  ## (this achieves the same thing as b = X \ Y)
  ## The number of free parameters (i.e. intercept + coefficients) is equal
  ## to n - dfe
  [Q, R] = qr (X, 0);
  b = R \ Q' * Y;

  ## Get fitted values
  fit = X * b;
  ## Get residuals from the fit
  resid = Y - fit;
  ## Calculate residual sums-of-squares
  sse = sum ((resid).^2);
  ## Calculate unscaled covariance matrix
  if (nargout > 3)
    ucov = R \ Q' / X';
  end

endfunction

%!demo
%!
%! # Two-sample unpaired test on independent samples (equivalent to Student's
%! # t-test). Note that the absolute value of t-statistic can be obtained by
%! # taking the square root of the reported F statistic. In this example,
%! # t = sqrt (1.44) = 1.20.
%!
%! score = [54 23 45 54 45 43 34 65 77 46 65]';
%! gender = {"male" "male" "male" "male" "male" "female" "female" "female" ...
%!           "female" "female" "female"}';
%!
%! [P, T] = anovan (score, gender, "display", "on", "varnames", "gender");

%!demo
%!
%! # Two-sample paired test on dependent or matched samples equivalent to a
%! # paired t-test. As for the first example, the t-statistic can be obtained by
%! # taking the square root of the reported F statistic. Note that the interaction
%! # between treatment x subject was dropped from the full model by assigning 
%! # subject as a random factor (').
%!
%! score = [4.5 5.6; 3.7 6.4; 5.3 6.4; 5.4 6.0; 3.9 5.7]';
%! treatment = {"before" "after"; "before" "after"; "before" "after";
%!              "before" "after"; "before" "after"}';
%! subject = {"GS" "GS"; "JM" "JM"; "HM" "HM"; "JW" "JW"; "PS" "PS"}';
%!
%! [P, T] = anovan (score(:), {treatment(:), subject(:)}, "model", "full", ...
%!                  "random", 2, "sstype", 2, "display", "on", ...
%!                  "varnames", {"treatment", "subject"});

%!demo
%!
%! # One-way ANOVA on the data from a study on the strength of structural beams,
%! # in Hogg and Ledolter (1987) Engineering Statistics. New York: MacMillan
%!
%! strength = [82 86 79 83 84 85 86 87 74 82 ...
%!            78 75 76 77 79 79 77 78 82 79]';
%! alloy = {"st","st","st","st","st","st","st","st",...
%!          "al1","al1","al1","al1","al1","al1",...
%!          "al2","al2","al2","al2","al2","al2"}';
%!
%! [P, T] = anovan (strength, alloy, "display", "on", "varnames", "alloy");

%!demo
%!
%! # One-way repeated measures ANOVA on the data from a study on the number of
%! # words recalled by 10 subjects for three time condtions, in Loftus & Masson
%! # (1994) Psychon Bull Rev. 1(4):476-490, Table 2. Note that the interaction
%! # between seconds x subject was dropped from the full model by assigning 
%! # subject as a random factor (').
%!
%! words = [10 13 13; 6 8 8; 11 14 14; 22 23 25; 16 18 20; ...
%!          15 17 17; 1 1 4; 12 15 17;  9 12 12;  8 9 12];
%! seconds = [1 2 5; 1 2 5; 1 2 5; 1 2 5; 1 2 5; ...
%!            1 2 5; 1 2 5; 1 2 5; 1 2 5; 1 2 5;];
%! subject = [ 1  1  1;  2  2  2;  3  3  3;  4  4  4;  5  5  5; ...
%!             6  6  6;  7  7  7;  8  8  8;  9  9  9; 10 10 10];
%!
%! [P, T] = anovan (words(:), {seconds(:), subject(:)}, "model", "full", ...
%!                  "random", 2, "sstype", 2, "display", "on", ...
%!                  "varnames", {"seconds", "subject"});

%!demo
%!
%! # Balanced two-way ANOVA with interaction on the data from a study of popcorn
%! # brands and popper types, in Hogg and Ledolter (1987) Engineering Statistics.
%! # New York: MacMillan
%!
%! popcorn = [5.5, 4.5, 3.5; 5.5, 4.5, 4.0; 6.0, 4.0, 3.0; ...
%!            6.5, 5.0, 4.0; 7.0, 5.5, 5.0; 7.0, 5.0, 4.5];
%! brands = {"Gourmet", "National", "Generic"; ...
%!           "Gourmet", "National", "Generic"; ...
%!           "Gourmet", "National", "Generic"; ...
%!           "Gourmet", "National", "Generic"; ...
%!           "Gourmet", "National", "Generic"; ...
%!           "Gourmet", "National", "Generic"};
%! popper = {"oil", "oil", "oil"; "oil", "oil", "oil"; "oil", "oil", "oil"; ...
%!           "air", "air", "air"; "air", "air", "air"; "air", "air", "air"};
%!
%! [P, T] = anovan (popcorn(:), {brands(:), popper(:)}, "display", "on",...
%!                  "model", "full", "varnames", {"brands", "popper"});

%!demo
%!
%! # Unbalanced two-way ANOVA (2x2) on the data from a study on the effects of
%! # gender and having a college degree on salaries of company employees,
%! # in Maxwell, Delaney and Kelly (2018): Chapter 7, Table 15
%!
%! salary = [24 26 25 24 27 24 27 23 15 17 20 16, ...
%!           25 29 27 19 18 21 20 21 22 19]';
%! gender = {"f" "f" "f" "f" "f" "f" "f" "f" "f" "f" "f" "f"...
%!           "m" "m" "m" "m" "m" "m" "m" "m" "m" "m"}';
%! degree = [1 1 1 1 1 1 1 1 0 0 0 0 1 1 1 0 0 0 0 0 0 0]';
%!
%! [P, T] = anovan (salary, {gender, degree}, "model", "full", ...
%!                  "sstype", 3, "display", "on", "varnames", ...
%!                  {"gender", "degree"});

%!demo
%!
%! # Unbalanced two-way ANOVA (3x2) on the data from a study of the effect of
%! # adding sugar and/or milk on the tendency of coffee to make people babble,
%! # in from Navarro (2019): 16.10
%!
%! sugar = {"real" "fake" "fake" "real" "real" "real" "none" "none" "none" ...
%!          "fake" "fake" "fake" "real" "real" "real" "none" "none" "fake"}';
%! milk = {"yes" "no" "no" "yes" "yes" "no" "yes" "yes" "yes" ...
%!         "no" "no" "yes" "no" "no" "no" "no" "no" "yes"}';
%! babble = [4.6 4.4 3.9 5.6 5.1 5.5 3.9 3.5 3.7...
%!           5.6 4.7 5.9 6.0 5.4 6.6 5.8 5.3 5.7]';
%!
%! [P, T] = anovan (babble, {sugar, milk}, "model", "full", "sstype", 3,...
%!                 "display", "on", "varnames", {"sugar", "milk"});

%!demo
%!
%! # Unbalanced three-way ANOVA (3x2x2) on the data from a study of the effects
%! # of three different drugs, biofeedback and diet on patient blood pressure,
%! # adapted* from Maxwell, Delaney and Kelly (2018): Chapter 8, Table 12
%! # * Missing values introduced to make the sample sizes unequal to test the
%! #   calculation of different types of sums-of-squares
%!
%! drug = {"X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" ...
%!         "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X";
%!         "Y" "Y" "Y" "Y" "Y" "Y" "Y" "Y" "Y" "Y" "Y" "Y" ...
%!         "Y" "Y" "Y" "Y" "Y" "Y" "Y" "Y" "Y" "Y" "Y" "Y";
%!         "Z" "Z" "Z" "Z" "Z" "Z" "Z" "Z" "Z" "Z" "Z" "Z" ...
%!         "Z" "Z" "Z" "Z" "Z" "Z" "Z" "Z" "Z" "Z" "Z" "Z"};
%! feedback = [1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0;
%!             1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0;
%!             1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0];
%! diet = [0 0 0 0 0 0 1 1 1 1 1 1 0 0 0 0 0 0 1 1 1 1 1 1;
%!         0 0 0 0 0 0 1 1 1 1 1 1 0 0 0 0 0 0 1 1 1 1 1 1;
%!         0 0 0 0 0 0 1 1 1 1 1 1 0 0 0 0 0 0 1 1 1 1 1 1];
%! BP = [170 175 165 180 160 158 161 173 157 152 181 190 ...
%!       173 194 197 190 176 198 164 190 169 164 176 175;
%!       186 194 201 215 219 209 164 166 159 182 187 174 ...
%!       189 194 217 206 199 195 171 173 196 199 180 NaN;
%!       180 187 199 170 204 194 162 184 183 156 180 173 ...
%!       202 228 190 206 224 204 205 199 170 160 NaN NaN];
%!
%! [P, T] = anovan (BP(:), {drug(:), feedback(:), diet(:)}, ...
%!                 "model", "full", "sstype", 3, "display", "on", ...
%!                 "varnames", {"drug", "feedback", "diet"});

%!demo
%!
%! # Balanced three-way ANOVA (2x2x2) with one of the factors being a blocking
%! # factor. The data is from a randomized block design study on the effects
%! # of antioxidant treatment on glutathione-S-transferase (GST) levels in
%! # different mouse strains, from Festing (2014), ILAR Journal, 55(3):427-476.
%! # Note that all interactions involving block were dropped from the full model 
%! # by assigning block as a random factor (').
%!
%! measurement = [444 614 423 625 408  856 447 719 ...
%!                764 831 586 782 609 1002 606 766]';
%! block = [1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2]';
%! strain= {"NIH","NIH","BALB/C","BALB/C","A/J","A/J","129/Ola","129/Ola",...
%!          "NIH","NIH","BALB/C","BALB/C","A/J","A/J","129/Ola","129/Ola"}';
%! treatment={"C" "T" "C" "T" "C" "T" "C" "T" "C" "T" "C" "T" "C" "T" "C" "T"}';
%!
%! [P, T] = anovan (measurement/10, {block, strain, treatment}, "sstype", 2, ...
%!                  "model", "full", "random", 1, "display", "on", ...
%!                  "varnames", {"block", "strain", "treatment"});

%!demo
%!
%! # One-way ANCOVA on data from a study of the additive effects of species
%! # and temperature on chirpy pulses of crickets, from Stitch, The Worst Stats
%! # Text eveR
%!
%! pulse = [67.9 65.1 77.3 78.7 79.4 80.4 85.8 86.6 87.5 89.1 ...
%!          98.6 100.8 99.3 101.7 44.3 47.2 47.6 49.6 50.3 51.8 ...
%!          60 58.5 58.9 60.7 69.8 70.9 76.2 76.1 77 77.7 84.7]';
%! temp = [20.8 20.8 24 24 24 24 26.2 26.2 26.2 26.2 28.4 ...
%!         29 30.4 30.4 17.2 18.3 18.3 18.3 18.9 18.9 20.4 ...
%!         21 21 22.1 23.5 24.2 25.9 26.5 26.5 26.5 28.6]';
%! species = {'ex' 'ex' 'ex' 'ex' 'ex' 'ex' 'ex' 'ex' 'ex' 'ex' 'ex' ...
%!            'ex' 'ex' 'ex' 'niv' 'niv' 'niv' 'niv' 'niv' 'niv' 'niv' ...
%!            'niv' 'niv' 'niv' 'niv' 'niv' 'niv' 'niv' 'niv' 'niv' 'niv'};
%!
%! [P, T, STATS] = anovan (pulse, {species, temp}, 'model', 'linear', ...
%!                         'continuous', 2, 'sstype', 'h', 'display', 'on', ...
%!                         'varnames', {"species", "temp"});

%!demo
%!
%! # Factorial ANCOVA on data from a study of the effects of treatment and
%! # exercise on stress reduction score after adjusting for age. Data from R
%! # datarium package).
%!
%! score = [95.6 82.2 97.2 96.4 81.4 83.6 89.4 83.8 83.3 85.7 ...
%!          97.2 78.2 78.9 91.8 86.9 84.1 88.6 89.8 87.3 85.4 ...
%!          81.8 65.8 68.1 70.0 69.9 75.1 72.3 70.9 71.5 72.5 ...
%!          84.9 96.1 94.6 82.5 90.7 87.0 86.8 93.3 87.6 92.4 ...
%!          100. 80.5 92.9 84.0 88.4 91.1 85.7 91.3 92.3 87.9 ...
%!          91.7 88.6 75.8 75.7 75.3 82.4 80.1 86.0 81.8 82.5]';
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
%!        58 56 57 59 59 60 55 53 55 58 68 62 61 54 59 63 60 67 60 67 ...
%!        75 54 57 62 65 60 58 61 65 57 56 58 58 58 52 53 60 62 61 61]';
%!
%! [P, T, STATS] = anovan (score, {treatment, exercise, age}, ...
%!                 'model', [1 0 0; 0 1 0; 0 0 1; 1 1 0], ...
%!                 'continuous', 3, 'sstype', 'h', 'display', 'on', ...
%!                 'varnames', {'treatment','exercise','age'});

%!demo
%!
%! # Unbalanced one-way ANOVA with custom, orthogonal contrasts. The ANOVA
%! # table is displayed followed by a matrix with columns from left-to-right
%! # corresponding to regression coefficients and their standard errors,
%! # t-statistics and p-values
%!
%! dv =  [ 8.706 10.362 11.552  6.941 10.983 10.092  6.421 14.943 15.931 ...
%!        22.968 18.590 16.567 15.944 21.637 14.492 17.965 18.851 22.891 ...
%!        22.028 16.884 17.252 18.325 25.435 19.141 21.238 22.196 18.038 ...
%!        22.628 31.163 26.053 24.419 32.145 28.966 30.207 29.142 33.212 ...
%!        25.694 ]';
%! g = [1 1 1 1 1 1 1 1 2 2 2 2 2 3 3 3 3 3 3 3 3 ...
%!      4 4 4 4 4 4 4 5 5 5 5 5 5 5 5 5]';
%! C = [ 0.4001601  0.3333333  0.5  0.0
%!       0.4001601  0.3333333 -0.5  0.0
%!       0.4001601 -0.6666667  0.0  0.0
%!      -0.6002401  0.0000000  0.0  0.5
%!      -0.6002401  0.0000000  0.0 -0.5];
%!
%! [P,T,STATS] = anovan (dv, g, "contrasts", {C}, "varnames", "score", ...
%!                       "display", "on");
%! disp (STATS.coeffnames)
%! disp (STATS.coeffs)

## Test 1 for anovan example 1
## Test compares anovan to results from MATLAB's anovan and ttest2 functions
%!test
%! score = [54 23 45 54 45 43 34 65 77 46 65]';
%! gender = {'male' 'male' 'male' 'male' 'male' 'female' 'female' 'female' ...
%!           'female' 'female' 'female'}';
%!
%! [P, T, STATS] = anovan (score,gender,'display','off');
%! assert (P(1), 0.2612876773271042,  1e-09);              # compared to p calculated by MATLAB anovan
%! assert (sqrt(T{2,6}), abs(1.198608733288208),  1e-09);  # compared to abs(t) calculated from sqrt(F) by MATLAB anovan
%! assert (P(1), 0.2612876773271047,  1e-09);              # compared to p calculated by MATLAB ttest2
%! assert (sqrt(T{2,6}), abs(-1.198608733288208),  1e-09); # compared to abs(t) calculated by MATLAB ttest2

## Test 2 for anovan example 2
## Test compares anovan to results from MATLAB's anovan and ttest functions
%!test
%! score = [4.5 5.6; 3.7 6.4; 5.3 6.4; 5.4 6.0; 3.9 5.7]';
%! treatment = {'before' 'after'; 'before' 'after'; 'before' 'after';
%!              'before' 'after'; 'before' 'after'}';
%! subject = {'GS' 'GS'; 'JM' 'JM'; 'HM' 'HM'; 'JW' 'JW'; 'PS' 'PS'}';
%!
%! [P, T, STATS] = anovan (score(:),{treatment(:),subject(:)},'display','off','sstype',2);
%! assert (P(1), 0.016004356735364,  1e-09);              # compared to p calculated by MATLAB anovan
%! assert (sqrt(T{2,6}), abs(4.00941576558195),  1e-09);  # compared to abs(t) calculated from sqrt(F) by MATLAB anovan
%! assert (P(1), 0.016004356735364,  1e-09);              # compared to p calculated by MATLAB ttest2
%! assert (sqrt(T{2,6}), abs(-4.00941576558195),  1e-09); # compared to abs(t) calculated by MATLAB ttest2

## Test 3 for anovan example 3
## Test compares anovan to results from MATLAB's anovan and anova1 functions
%!test
%! strength = [82 86 79 83 84 85 86 87 74 82 ...
%!            78 75 76 77 79 79 77 78 82 79]';
%! alloy = {'st','st','st','st','st','st','st','st',...
%!          'al1','al1','al1','al1','al1','al1',...
%!          'al2','al2','al2','al2','al2','al2'}';
%!
%! [P, T, STATS] = anovan (strength,{alloy},'display','off');
%! assert (P(1), 0.000152643638830491,  1e-09);
%! assert (T{2,6}, 15.4,  1e-09);

## Test 4 for anovan example 4
## Test compares anovan to results from MATLAB's anovan function
%!test
%! words = [10 13 13; 6 8 8; 11 14 14; 22 23 25; 16 18 20; ...
%!          15 17 17; 1 1 4; 12 15 17;  9 12 12;  8 9 12];
%! subject = [ 1  1  1;  2  2  2;  3  3  3;  4  4  4;  5  5  5; ...
%!             6  6  6;  7  7  7;  8  8  8;  9  9  9; 10 10 10];
%! seconds = [1 2 5; 1 2 5; 1 2 5; 1 2 5; 1 2 5; ...
%!            1 2 5; 1 2 5; 1 2 5; 1 2 5; 1 2 5;];
%! [P, T, STATS] = anovan (words(:),{seconds(:),subject(:)},'model','full','random',2,'sstype',2,'display','off');
%! assert (P(1), 1.51865926758752e-07,  1e-09);
%! assert (P(2), 1.49150337808586e-15,  1e-09);
%! assert (T{2,2}, 52.2666666666667,  1e-09);
%! assert (T{3,2}, 942.533333333333,  1e-09);
%! assert (T{4,2}, 11.0666666666667,  1e-09);

## Test 5 for anovan example 5
## Test compares anovan to results from MATLAB's anovan function
%!test
%! popcorn = [5.5, 4.5, 3.5; 5.5, 4.5, 4.0; 6.0, 4.0, 3.0; ...
%!            6.5, 5.0, 4.0; 7.0, 5.5, 5.0; 7.0, 5.0, 4.5];
%! brands = {'Gourmet', 'National', 'Generic'; ...
%!           'Gourmet', 'National', 'Generic'; ...
%!           'Gourmet', 'National', 'Generic'; ...
%!           'Gourmet', 'National', 'Generic'; ...
%!           'Gourmet', 'National', 'Generic'; ...
%!           'Gourmet', 'National', 'Generic'};
%! popper = {'oil', 'oil', 'oil'; 'oil', 'oil', 'oil'; 'oil', 'oil', 'oil'; ...
%!           'air', 'air', 'air'; 'air', 'air', 'air'; 'air', 'air', 'air'};
%!
%! [P, T, STATS] = anovan (popcorn(:),{brands(:),popper(:)},'display','off','model','full');
%! assert (P(1), 7.67895738278171e-07,  1e-09);
%! assert (P(2), 0.000100373896304998,  1e-09);
%! assert (P(3), 0.746215396636649,  1e-09);
%! assert (T{2,6}, 56.7,  1e-09);
%! assert (T{3,6}, 32.4,  1e-09);
%! assert (T{4,6}, 0.29999999999997,  1e-09);

## Test 6 for anovan example 6
## Test compares anovan to results from MATLAB's anovan function
%!test
%! salary = [24 26 25 24 27 24 27 23 15 17 20 16, ...
%!           25 29 27 19 18 21 20 21 22 19]';
%! gender = {'f' 'f' 'f' 'f' 'f' 'f' 'f' 'f' 'f' 'f' 'f' 'f'...
%!           'm' 'm' 'm' 'm' 'm' 'm' 'm' 'm' 'm' 'm'}';
%! degree = [1 1 1 1 1 1 1 1 0 0 0 0 1 1 1 0 0 0 0 0 0 0]';
%! [P, T, STATS] = anovan (salary,{gender,degree},'model','full','sstype',1,'display','off');
%! assert (P(1), 0.747462549227232,  1e-09);
%! assert (P(2), 1.03809316857694e-08,  1e-09);
%! assert (P(3), 0.523689833702691,  1e-09);
%! assert (T{2,2}, 0.296969696969699,  1e-09);
%! assert (T{3,2}, 272.391841491841,  1e-09);
%! assert (T{4,2}, 1.17482517482512,  1e-09);
%! assert (T{5,2}, 50.0000000000001,  1e-09);
%! [P, T, STATS] = anovan (salary,{degree,gender},'model','full','sstype',1,'display','off');
%! assert (P(1), 2.53445097305047e-08,  1e-09);
%! assert (P(2), 0.00388133678528749,  1e-09);
%! assert (P(3), 0.523689833702671,  1e-09);
%! assert (T{2,2}, 242.227272727273,  1e-09);
%! assert (T{3,2}, 30.4615384615384,  1e-09);
%! assert (T{4,2}, 1.17482517482523,  1e-09);
%! assert (T{5,2}, 50.0000000000001,  1e-09);
%! [P, T, STATS] = anovan (salary,{gender,degree},'model','full','sstype',2,'display','off');
%! assert (P(1), 0.00388133678528743,  1e-09);
%! assert (P(2), 1.03809316857694e-08,  1e-09);
%! assert (P(3), 0.523689833702691,  1e-09);
%! assert (T{2,2}, 30.4615384615385,  1e-09);
%! assert (T{3,2}, 272.391841491841,  1e-09);
%! assert (T{4,2}, 1.17482517482512,  1e-09);
%! assert (T{5,2}, 50.0000000000001,  1e-09);
%! [P, T, STATS] = anovan (salary,{gender,degree},'model','full','sstype',3,'display','off');
%! assert (P(1), 0.00442898146583742,  1e-09);
%! assert (P(2), 1.30634252053587e-08,  1e-09);
%! assert (P(3), 0.523689833702691,  1e-09);
%! assert (T{2,2}, 29.3706293706294,  1e-09);
%! assert (T{3,2}, 264.335664335664,  1e-09);
%! assert (T{4,2}, 1.17482517482512,  1e-09);
%! assert (T{5,2}, 50.0000000000001,  1e-09);

## Test 7 for anovan example 7
## Test compares anovan to results from MATLAB's anovan function
%!test
%! sugar = {'real' 'fake' 'fake' 'real' 'real' 'real' 'none' 'none' 'none' ...
%!          'fake' 'fake' 'fake' 'real' 'real' 'real' 'none' 'none' 'fake'}';
%! milk = {'yes' 'no' 'no' 'yes' 'yes' 'no' 'yes' 'yes' 'yes' ...
%!         'no' 'no' 'yes' 'no' 'no' 'no' 'no' 'no' 'yes'}';
%! babble = [4.6 4.4 3.9 5.6 5.1 5.5 3.9 3.5 3.7...
%!           5.6 4.7 5.9 6.0 5.4 6.6 5.8 5.3 5.7]';
%! [P, T, STATS] = anovan (babble,{sugar,milk},'model','full','sstype',1,'display','off');
%! assert (P(1), 0.0108632139833963,  1e-09);
%! assert (P(2), 0.0810606976703546,  1e-09);
%! assert (P(3), 0.00175433329935627,  1e-09);
%! assert (T{2,2}, 3.55752380952381,  1e-09);
%! assert (T{3,2}, 0.956108477471702,  1e-09);
%! assert (T{4,2}, 5.94386771300448,  1e-09);
%! assert (T{5,2}, 3.1625,  1e-09);
%! [P, T, STATS] = anovan (babble,{milk,sugar},'model','full','sstype',1,'display','off');
%! assert (P(1), 0.0373333189297505,  1e-09);
%! assert (P(2), 0.017075098787169,  1e-09);
%! assert (P(3), 0.00175433329935627,  1e-09);
%! assert (T{2,2}, 1.444,  1e-09);
%! assert (T{3,2}, 3.06963228699552,  1e-09);
%! assert (T{4,2}, 5.94386771300448,  1e-09);
%! assert (T{5,2}, 3.1625,  1e-09);
%! [P, T, STATS] = anovan (babble,{sugar,milk},'model','full','sstype',2,'display','off');
%! assert (P(1), 0.017075098787169,  1e-09);
%! assert (P(2), 0.0810606976703546,  1e-09);
%! assert (P(3), 0.00175433329935627,  1e-09);
%! assert (T{2,2}, 3.06963228699552,  1e-09);
%! assert (T{3,2}, 0.956108477471702,  1e-09);
%! assert (T{4,2}, 5.94386771300448,  1e-09);
%! assert (T{5,2}, 3.1625,  1e-09);
%! [P, T, STATS] = anovan (babble,{sugar,milk},'model','full','sstype',3,'display','off');
%! assert (P(1), 0.0454263063473954,  1e-09);
%! assert (P(2), 0.0746719907091438,  1e-09);
%! assert (P(3), 0.00175433329935627,  1e-09);
%! assert (T{2,2}, 2.13184977578476,  1e-09);
%! assert (T{3,2}, 1.00413461538462,  1e-09);
%! assert (T{4,2}, 5.94386771300448,  1e-09);
%! assert (T{5,2}, 3.1625,  1e-09);

## Test 8 for anovan example 8
## Test compares anovan to results from MATLAB's anovan function
%!test
%! drug = {'X' 'X' 'X' 'X' 'X' 'X' 'X' 'X' 'X' 'X' 'X' 'X' ...
%!         'X' 'X' 'X' 'X' 'X' 'X' 'X' 'X' 'X' 'X' 'X' 'X';
%!         'Y' 'Y' 'Y' 'Y' 'Y' 'Y' 'Y' 'Y' 'Y' 'Y' 'Y' 'Y' ...
%!         'Y' 'Y' 'Y' 'Y' 'Y' 'Y' 'Y' 'Y' 'Y' 'Y' 'Y' 'Y';
%!         'Z' 'Z' 'Z' 'Z' 'Z' 'Z' 'Z' 'Z' 'Z' 'Z' 'Z' 'Z' ...
%!         'Z' 'Z' 'Z' 'Z' 'Z' 'Z' 'Z' 'Z' 'Z' 'Z' 'Z' 'Z'};
%! feedback = [1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0;
%!             1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0;
%!             1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0];
%! diet = [0 0 0 0 0 0 1 1 1 1 1 1 0 0 0 0 0 0 1 1 1 1 1 1;
%!         0 0 0 0 0 0 1 1 1 1 1 1 0 0 0 0 0 0 1 1 1 1 1 1;
%!         0 0 0 0 0 0 1 1 1 1 1 1 0 0 0 0 0 0 1 1 1 1 1 1];
%! BP = [170 175 165 180 160 158 161 173 157 152 181 190 ...
%!       173 194 197 190 176 198 164 190 169 164 176 175;
%!       186 194 201 215 219 209 164 166 159 182 187 174 ...
%!       189 194 217 206 199 195 171 173 196 199 180 NaN;
%!       180 187 199 170 204 194 162 184 183 156 180 173 ...
%!       202 228 190 206 224 204 205 199 170 160 NaN NaN];
%! [P, T, STATS] = anovan (BP(:),{drug(:),feedback(:),diet(:)},'model','full','sstype', 1,'display','off');
%! assert (P(1), 7.02561843825325e-05,  1e-09);
%! assert (P(2), 0.000425806013389362,  1e-09);
%! assert (P(3), 6.16780773446401e-07,  1e-09);
%! assert (P(4), 0.261347622678438,  1e-09);
%! assert (P(5), 0.0542278432357043,  1e-09);
%! assert (P(6), 0.590353225626655,  1e-09);
%! assert (P(7), 0.0861628249564267,  1e-09);
%! assert (T{2,2}, 3614.70355731226,  1e-09);
%! assert (T{3,2}, 2227.46639771024,  1e-09);
%! assert (T{4,2}, 5008.25614451819,  1e-09);
%! assert (T{5,2}, 437.066007908781,  1e-09);
%! assert (T{6,2}, 976.180770397332,  1e-09);
%! assert (T{7,2}, 46.616653365254,  1e-09);
%! assert (T{8,2}, 814.345251396648,  1e-09);
%! assert (T{9,2}, 9065.8,  1e-09);
%! [P, T, STATS] = anovan (BP(:),{drug(:),feedback(:),diet(:)},'model','full','sstype',2,'display','off');
%! assert (P(1), 9.4879638470754e-05,  1e-09);
%! assert (P(2), 0.00124177666315809,  1e-09);
%! assert (P(3), 6.86162012732911e-07,  1e-09);
%! assert (P(4), 0.260856132341256,  1e-09);
%! assert (P(5), 0.0523758623892078,  1e-09);
%! assert (P(6), 0.590353225626655,  1e-09);
%! assert (P(7), 0.0861628249564267,  1e-09);
%! assert (T{2,2}, 3481.72176560122,  1e-09);
%! assert (T{3,2}, 1837.08812970469,  1e-09);
%! assert (T{4,2}, 4957.20277938622,  1e-09);
%! assert (T{5,2}, 437.693674777847,  1e-09);
%! assert (T{6,2}, 988.431929811402,  1e-09);
%! assert (T{7,2}, 46.616653365254,  1e-09);
%! assert (T{8,2}, 814.345251396648,  1e-09);
%! assert (T{9,2}, 9065.8,  1e-09);
%! [P, T, STATS] = anovan (BP(:),{drug(:),feedback(:),diet(:)},'model','full','sstype', 3,'display','off');
%! assert (P(1), 0.000106518678028207,  1e-09);
%! assert (P(2), 0.00125371366571508,  1e-09);
%! assert (P(3), 5.30813260778464e-07,  1e-09);
%! assert (P(4), 0.308353667232981,  1e-09);
%! assert (P(5), 0.0562901327343161,  1e-09);
%! assert (P(6), 0.599091042141092,  1e-09);
%! assert (P(7), 0.0861628249564267,  1e-09);
%! assert (T{2,2}, 3430.88156424581,  1e-09);
%! assert (T{3,2}, 1833.68031496063,  1e-09);
%! assert (T{4,2}, 5080.48346456693,  1e-09);
%! assert (T{5,2}, 382.07709497207,  1e-09);
%! assert (T{6,2}, 963.037988826813,  1e-09);
%! assert (T{7,2}, 44.4519685039322,  1e-09);
%! assert (T{8,2}, 814.345251396648,  1e-09);
%! assert (T{9,2}, 9065.8,  1e-09);

## Test 9 for anovan example 9
## Test compares anovan to results from MATLAB's anovan function
%!test
%! measurement = [444 614 423 625 408  856 447 719 ...
%!                764 831 586 782 609 1002 606 766]';
%! block = [1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2]';
%! strain= {'NIH','NIH','BALB/C','BALB/C','A/J','A/J','129/Ola','129/Ola',...
%!          'NIH','NIH','BALB/C','BALB/C','A/J','A/J','129/Ola','129/Ola'}';
%! treatment={'C' 'T' 'C' 'T' 'C' 'T' 'C' 'T' 'C' 'T' 'C' 'T' 'C' 'T' 'C' 'T'}';
%! [P, T, STATS] = anovan (measurement/10,{block,strain,treatment},'model','full','random',1,'display','off');
%! assert (P(1), 0.000339814602130731,  1e-09);
%! assert (P(2), 0.0914352969909372,  1e-09);
%! assert (P(3), 5.04077373924908e-05,  1e-09);
%! assert (P(4), 0.0283196918836667,  1e-09);
%! assert (T{2,2}, 1242.5625,  1e-09);
%! assert (T{3,2}, 286.132500000002,  1e-09);
%! assert (T{4,2}, 2275.29,  1e-09);
%! assert (T{5,2}, 495.905000000001,  1e-09);
%! assert (T{6,2}, 207.007499999999,  1e-09);

## Test 10 for anovan example 10
## Test compares anovan to results from MATLAB's anovan function
%!test
%! pulse = [67.9 65.1 77.3 78.7 79.4 80.4 85.8 86.6 87.5 89.1 ...
%!          98.6 100.8 99.3 101.7 44.3 47.2 47.6 49.6 50.3 51.8 ...
%!          60 58.5 58.9 60.7 69.8 70.9 76.2 76.1 77 77.7 84.7]';
%! temp = [20.8 20.8 24 24 24 24 26.2 26.2 26.2 26.2 28.4 ...
%!         29 30.4 30.4 17.2 18.3 18.3 18.3 18.9 18.9 20.4 ...
%!         21 21 22.1 23.5 24.2 25.9 26.5 26.5 26.5 28.6]';
%! species = {'ex' 'ex' 'ex' 'ex' 'ex' 'ex' 'ex' 'ex' 'ex' 'ex' 'ex' ...
%!            'ex' 'ex' 'ex' 'niv' 'niv' 'niv' 'niv' 'niv' 'niv' 'niv' ...
%!            'niv' 'niv' 'niv' 'niv' 'niv' 'niv' 'niv' 'niv' 'niv' 'niv'};
%! [P, T, STATS] = anovan (pulse,{species,temp},'model','linear','continuous',2,'sstype','h','display','off');
%! assert (P(1), 6.27153318786007e-14,  1e-09);
%! assert (P(2), 2.48773241196644e-25,  1e-09);
%! assert (T{2,2}, 598.003953318404,  1e-09);
%! assert (T{3,2}, 4376.08256843712,  1e-09);
%! assert (T{4,2}, 89.3498685376726,  1e-09);
%! assert (T{2,6}, 187.399388123951,  1e-09);
%! assert (T{3,6}, 1371.35413763454,  1e-09);

## Test 11 for anovan example 11
## Test compares anovan to results from MATLAB's anovan function
%!test
%! score = [95.6 82.2 97.2 96.4 81.4 83.6 89.4 83.8 83.3 85.7 ...
%!          97.2 78.2 78.9 91.8 86.9 84.1 88.6 89.8 87.3 85.4 ...
%!          81.8 65.8 68.1 70.0 69.9 75.1 72.3 70.9 71.5 72.5 ...
%!          84.9 96.1 94.6 82.5 90.7 87.0 86.8 93.3 87.6 92.4 ...
%!          100. 80.5 92.9 84.0 88.4 91.1 85.7 91.3 92.3 87.9 ...
%!          91.7 88.6 75.8 75.7 75.3 82.4 80.1 86.0 81.8 82.5]';
%! treatment = {'yes' 'yes' 'yes' 'yes' 'yes' 'yes' 'yes' 'yes' 'yes' 'yes' ...
%!              'yes' 'yes' 'yes' 'yes' 'yes' 'yes' 'yes' 'yes' 'yes' 'yes' ...
%!              'yes' 'yes' 'yes' 'yes' 'yes' 'yes' 'yes' 'yes' 'yes' 'yes' ...
%!              'no'  'no'  'no'  'no'  'no'  'no'  'no'  'no'  'no'  'no'  ...
%!              'no'  'no'  'no'  'no'  'no'  'no'  'no'  'no'  'no'  'no'  ...
%!              'no'  'no'  'no'  'no'  'no'  'no'  'no'  'no'  'no'  'no'}';
%! exercise = {'lo'  'lo'  'lo'  'lo'  'lo'  'lo'  'lo'  'lo'  'lo'  'lo'  ...
%!             'mid' 'mid' 'mid' 'mid' 'mid' 'mid' 'mid' 'mid' 'mid' 'mid' ...
%!             'hi'  'hi'  'hi'  'hi'  'hi'  'hi'  'hi'  'hi'  'hi'  'hi'  ...
%!             'lo'  'lo'  'lo'  'lo'  'lo'  'lo'  'lo'  'lo'  'lo'  'lo'  ...
%!             'mid' 'mid' 'mid' 'mid' 'mid' 'mid' 'mid' 'mid' 'mid' 'mid' ...
%!             'hi'  'hi'  'hi'  'hi'  'hi'  'hi'  'hi'  'hi'  'hi'  'hi'}';
%! age = [59 65 70 66 61 65 57 61 58 55 62 61 60 59 55 57 60 63 62 57 ...
%!        58 56 57 59 59 60 55 53 55 58 68 62 61 54 59 63 60 67 60 67 ...
%!        75 54 57 62 65 60 58 61 65 57 56 58 58 58 52 53 60 62 61 61]';
%! [P, T, STATS] = anovan (score,{treatment,exercise,age},'model','full','continuous',3,'sstype','h','display','off');
%! assert (P(5), 0.9245630968248468,  1e-09);
%! assert (P(6), 0.791115159521822,  1e-09);
%! assert (P(7), 0.9296668751457956,  1e-09);
%! [P, T, STATS] = anovan (score,{treatment,exercise,age},'model',[1 0 0; 0 1 0; 0 0 1; 1 1 0],'continuous',3,'sstype','h','display','off');
%! assert (P(1), 0.00158132928938933,  1e-09);
%! assert (P(2), 2.12537505039986e-07,  1e-09);
%! assert (P(3), 0.00390292555160047,  1e-09);
%! assert (P(4), 0.0164086580775543,  1e-09);
%! assert (T{2,6}, 11.0956027650549,  1e-09);
%! assert (T{3,6}, 20.8195665467178,  1e-09);
%! assert (T{4,6}, 9.10966630720186,  1e-09);
%! assert (T{5,6}, 4.4457923698584,  1e-09);

## Test 12 for anovan example 12
## Test compares anovan regression coefficients to R:
## https://www.uvm.edu/~statdhtx/StatPages/Unequal-ns/Unequal_n%27s_contrasts.html
%!test
%! dv =  [ 8.706 10.362 11.552  6.941 10.983 10.092  6.421 14.943 15.931 ...
%!        22.968 18.590 16.567 15.944 21.637 14.492 17.965 18.851 22.891 ...
%!        22.028 16.884 17.252 18.325 25.435 19.141 21.238 22.196 18.038 ...
%!        22.628 31.163 26.053 24.419 32.145 28.966 30.207 29.142 33.212 ...
%!        25.694 ]';
%! g = [1 1 1 1 1 1 1 1 2 2 2 2 2 3 3 3 3 3 3 3 3 4 4 4 4 4 4 4 5 5 5 5 5 5 5 5 5]';
%! C = [ 0.4001601  0.3333333  0.5  0.0
%!       0.4001601  0.3333333 -0.5  0.0
%!       0.4001601 -0.6666667  0.0  0.0
%!      -0.6002401  0.0000000  0.0  0.5
%!      -0.6002401  0.0000000  0.0 -0.5];
%! [P,T,STATS] = anovan (dv,g,'contrasts',{C},'display','off');
%! assert (STATS.coeffs(1,1), 19.4001,  1e-04);
%! assert (STATS.coeffs(2,1), -9.3297,  1e-04);
%! assert (STATS.coeffs(3,1), -5.0000,  1e-04);
%! assert (STATS.coeffs(4,1), -8.0000,  1e-04);
%! assert (STATS.coeffs(5,1), -8.0000,  1e-04);
%! assert (STATS.coeffs(1,2), 0.4831,  1e-04);
%! assert (STATS.coeffs(2,2), 0.9694,  1e-04);
%! assert (STATS.coeffs(3,2), 1.3073,  1e-04);
%! assert (STATS.coeffs(4,2), 1.6411,  1e-04);
%! assert (STATS.coeffs(5,2), 1.4507,  1e-04);
%! assert (STATS.coeffs(1,3), 40.161,  1e-03);
%! assert (STATS.coeffs(2,3), -9.624,  1e-03);
%! assert (STATS.coeffs(3,3), -3.825,  1e-03);
%! assert (STATS.coeffs(4,3), -4.875,  1e-03);
%! assert (STATS.coeffs(5,3), -5.515,  1e-03);
%! assert (STATS.coeffs(2,4), 5.74e-11,  1e-12);
%! assert (STATS.coeffs(3,4), 0.000572,  1e-06);
%! assert (STATS.coeffs(4,4), 2.86e-05,  1e-07);
%! assert (STATS.coeffs(5,4), 4.44e-06,  1e-08);
