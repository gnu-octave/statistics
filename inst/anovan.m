## Copyright (C) 2022 Andrew Penn <A.C.Penn@sussex.ac.uk>
## Copyright (C) 2022 Andreas Bertsatos <abertsatos@biol.uoa.gr>
## Copyright (C) 2021 Christian Scholz
## Copyright (C) 2003-2005 Andy Adler <adler@ncf.ca>
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
## Perform a multi (N)-way analysis of variance (ANOVA) to evaluate the effect  
## of a categorical predictor, or multiple predictors, on a continuous outcome. 
## The algorithms used make @code{anovan} suitable for balanced or unbalanced designs. 
## A range of experimental designs can be analysed using @code{anovan}. Examples
## of function usage can be found by entering the command @code{demo anovan}. 
##
## 
## Data is a single vector @var{Y} with groups specified by a corresponding 
## matrix or cell array of group labels @var{GROUP}, where each column of 
## @var{GROUP} has the same number of rows as @var{Y}. For example, if 
## @var{Y} = [1.1;1.2]; @var{GROUP} = [1,2,1; 1,5,2]; then observation 1.1 was 
## measured under conditions 1,2,1 and observation 1.2 was measured under 
## conditions 1,5,2. Note that groups do not need to be sequentially numbered.
##
## @code{anovan} can take a number of optional parameters as name-value pairs.
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
## 1 : Type I sequential sums-of-squares
##
## @item
## 2 : Type II hierarchical (or partially sequential) sums-of-squares
##
## @item
## 3 (default): Type III constrained, marginal or orthogonal sums-of-squares
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
## @var{dispopt} can be either "on" (default) or "off" and switches the display
## of the ANOVA table.
##
## @code{anovan} can return up to four output arguments:
##
## @var{P} = anovan (@dots{}) returns a vector of p-values, one for each term.
##
## [@var{P}, @var{T}] = anovan (@dots{}) returns a cell array containing the
## ANOVA table.
##
## [@var{P}, @var{T}, @var{STATS}] = anovan (@dots{}) returns a structure 
## containing additional statistics, including coefficients of the linear model,
## the model residuals, and the number of levels in each factor.
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
      
    if nargin <= 1
      error ("anovan usage: ""anovan (Y, GROUP)""; atleast 2 input arguments required");
    endif

    # Check supplied parameters
    modeltype = "linear";
    display = "on";
    sstype = 3;
    varnames = [];
    for idx = 3:2:nargin
      name = varargin{idx-2};
      value = varargin{idx-1};
      switch lower (name)
        case "model"
          modeltype = value;
        case "varnames"
          varnames = value;
        case "display"
          display = value;   
        case "sstype"
          sstype = value;
        otherwise 
          error (sprintf("anovan: parameter %s is not supported", name));
      endswitch
    endfor
    
    # Accomodate for different formats for GROUP 
    # GROUP can be a matrix of numeric identifiers of a cell arrays of strings or numeric idenitiers
    N = size (GROUP,2); # number of anova "ways"
    n = numel (Y);      # total number of observations
    if (prod (size (Y)) ~= n)
      error ("anovan: for ""anovan (Y, GROUP)"", Y must be a vector");
    endif
    if iscell(GROUP)
      if (size(GROUP, 1) == 1)
        for k = 1:N
          if isnumeric (GROUP{k})
            tmp(:,k) = cellstr (num2str (GROUP{k}));
          else
            tmp(:,k) = GROUP{k};
          endif
        endfor
        GROUP = tmp;
      else
        for k = 1:N
          tmp(:,k) = cellstr (char (GROUP{k}));
        endfor
      endif
    endif
    if (size (GROUP,1) ~= n)
      error ("anovan: GROUP must be a matrix of the same number of rows as Y");
    endif
    if ~isempty (varnames) 
      if iscell (varnames)
        if all (cellfun (@ischar, varnames))
          nvarnames = numel(varnames);
        else
          error ("anovan: all variable names must be character or character arrays");
        endif
      elseif ischar (varnames)
        nvarnames = 1;
        varnames = {varnames};
      elseif isstring (varnames)
        nvarnames = 1;
        varnames = {char(varnames)};
      else
        error ("anovan: varnames is not of a valid type. Must be cell array of character arrays, character array or string");
      endif
    else
      nvarnames = N;
      varnames = arrayfun(@(x) ["X",num2str(x)], 1:N, "UniformOutput", 0);
    endif
    if (nvarnames ~= N)
      error ("anovan: number of variable names is not equal to number of grouping variables");
    endif
    
    # Remove NaN or non-finite observations
    excl = logical (isnan(Y) + isinf(Y));
    Y(excl) = [];
    GROUP(excl,:) = [];
    if (size (Y, 1) == 1)
      Y = Y.';         # if Y is a row vector, make it a column vector
    endif
    n = numel (Y);     # recalculate total number of observations

    # Evaluate model type input argument and create terms matrix if not provided
    msg = "anovan: the number of columns in the term definitions cannot exceed the number of columns of GROUP";
    if ischar (modeltype)
      switch lower(modeltype)
        case "linear"
          modeltype = 1;
        case "interaction"
          modeltype = 2;
        case "full"
          modeltype = N;
        otherwise
          error ("anovan: model type not recognised");
      endswitch
    endif
    if isscalar (modeltype)
      TERMS = cell (modeltype,1);
      v = false (1, N);
      switch lower (modeltype)
        case 1
          # Create term definitions for an additive linear model
          TERMS = eye (N);
        case 2
          # Create term definitions for a model with two factor interactions
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
          if modeltype > N
            error (msg);
          endif
          # Create term definitions for a full model
          Nx = zeros (1, N-1);
          Nx = 0;
          for k = 1:N
            Nx = Nx + nchoosek(N,k);
          endfor
          for j = 1:modeltype
            v(1:j) = 1;
            TERMS(j) = flipud (unique (perms (v), "rows"));
          endfor
          TERMS = cell2mat (TERMS);
      endswitch
    else
      # Assume that the user provided a suitable matrix of term definitions
      if (size (modeltype, 2) > N)
        error (msg);
      endif
      TERMS = logical (modeltype);
    endif
    # Evaluate terms matrix
    Ng = sum (TERMS, 2); 
    if any (diff (Ng) < 0)
      error ("anovan: the model terms matrix must list main effects above/before interactions");
    endif
    Nm = sum (Ng == 1);
    Nx = sum (Ng > 1);
    Nt = Nm + Nx;
    if any (any (TERMS(1:Nm,:)) ~= any (TERMS))
      error ("anovan: all factors involved in interactions must have a main effect");  
    end
    
    # Calculate total sum-of-squares
    ct  = sum (Y)^2 / n;   % correction term
    sst = sum (Y.^2) - ct;
    dft = n - 1;
    
    # Fit linear models, and calculate sums-of-squares for ANOVA
    switch lower (sstype)
      case 1
        # Type I sequential sums-of-squares (sstype = 1)
        R = sst;
        ss = zeros (Nt,1);
        [X, grpnames, nlevels, df, termcols] = make_design_matrix (GROUP, TERMS, n, Nm, Nx, Ng);
        for j = 1:Nt
          XS = cell2mat (X(1:j+1));
          [b, sse, resid] = lmfit (XS, Y);
          ss(j) = R - sse;
          R = sse;
        endfor
        sstype_char = "I";
      case 2
        # Type II (hierarchical, or partially sequential) sums of squares
        ss = zeros (Nt,1);
        [X, grpnames, nlevels, df, termcols]  = make_design_matrix (GROUP, TERMS, n, Nm, Nx, Ng);
        for j = 1:Nt
          i = find (TERMS(j,:)); 
          k = cat (1, 1, 1 + find (any (~TERMS(:,i),2)));
          XS = cell2mat (X(k)); 
          [jnk, R1] = lmfit (XS, Y);
          k = cat (1, j+1, k);
          XS = cell2mat (X(k));
          [jnk, R2] = lmfit (XS, Y);
          ss(j) = R1 - R2;
        end
        [b, sse, resid] = lmfit (cell2mat (X), Y);
        sstype_char = "II";
      case 3
        # Type III (constrained, marginal or orthogonal) sums of squares
        ss = zeros (Nt, 1);
        [X, grpnames, nlevels, df, termcols] = make_design_matrix (GROUP, TERMS, n, Nm, Nx, Ng);
        [b, sse, resid] = lmfit (cell2mat (X), Y);
        for j = 1:Nt
          XS = cell2mat (X(1:Nt+1 ~= j+1));
          [jnk, R, resid] = lmfit (XS, Y);
          ss(j) = R - sse;
        endfor
        sstype_char = "III";
      otherwise
        # sstype "h" not supported
        error ("anovan: only sstype 1, 2 and 3 are supported");
    endswitch
    dfe = dft - sum (df);
    ms = ss ./ df;
    mse = sse / dfe;
    eta_sq = ss ./ sst;
    partial_eta_sq = ss ./ (ss + sse);
    F = ms / mse;
    P = 1 - fcdf (F, df, dfe);

    # Prepare stats output structure
    # Note that the information provided by STATS is not sufficient for MATLAB's multcompare function
    STATS = struct ("source","anovan", ...
                    "resid", resid, ...
                    "coeffs", b, ...
                    "Rtr", [], ...           # Not used in Octave
                    "rowbasis", [], ...      # Not used in Octave
                    "dfe", dfe, ...
                    "mse", mse, ...
                    "nullproject", [], ...   # Not used in Octave
                    "terms", TERMS, ...
                    "nlevels", nlevels, ...  
                    "continuous", [], ...
                    "vmeans", [], ...        # Not used since "continuous" argument name not supported
                    "termcols", termcols, ...
                    "coeffnames", [], ...    # Not used in Octave
                    "vars", [], ...          # Not used in Octave
                    "varnames", {varnames}, ...
                    "grpnames", {grpnames}, ...
                    "vnested", [], ...       # Not used since "nested" argument name not supported
                    "ems", [], ...           # Not used since "nested" argument name not supported
                    "denom", [], ...         # Not used since "random" argument name not supported
                    "dfdenom", [], ...       # Not used since "random" argument name not supported
                    "msdenom", [], ...       # Not used since "random" argument name not supported
                    "varest", [], ...        # Not used since "random" argument name not supported
                    "varci", [], ...         # Not used since "random" argument name not supported
                    "txtdenom", [], ...      # Not used since "random" argument name not supported
                    "txtems", [], ...        # Not used since "random" argument name not supported
                    "rtnames", [],...        # Not used since "random" argument name not supported
                    "eta_squared", eta_sq,...
                    "partial_eta_squared", partial_eta_sq);
    
    # Prepare cell array containing the ANOVA table (T)
    T = cell (Nt + 3, 7);
    T(1,:) = {"Source","Sum Sq.","d.f.","Mean Sq.","Eta Sq.","F","Prob>F"};
    T(2:Nt+1,2:7) = num2cell([ss df ms partial_eta_sq F P]);
    T(end-1,1:4) = {"Error",sse,dfe,mse};
    T(end,1:3) = {"Total",sst,dft};
    for i = 1:Nt
      str = sprintf("%s*",varnames{find(TERMS(i,:))});
      T(i+1,1) = str(1:end-1);
    endfor
    
    # Print ANOVA table 
    switch lower (display)
      case "on"
        # Get dimensions of the ANOVA table
        [nrows, ncols] = size (T);
        # Print table
        fprintf("\n%d-way ANOVA table (Type %s sums of squares):\n\n", Nm, sstype_char);
        fprintf("Source                   Sum Sq.    d.f.    Mean Sq.  R Sq.            F  Prob>F\n");
        fprintf("********************************************************************************\n");  
        for i = 1:Nt
          str = T{i+1,1};
          l = numel(str);  # Needed to truncate source term name at 18 characters
          # Format and print the statistics for each model term
          # Format F statistics and p-values in APA style
          if (P(i) < 0.001)
            fprintf ("%-20s  %10.5g  %6d  %10.5g  %4.3f  %11.2f   <.001 \n", str(1:min(18,l)), T{i+1,2:end-1});
          elseif (P(i) < 1.0)
            fprintf ("%-20s  %10.5g  %6d  %10.5g  %4.3f  %11.2f    .%03u \n", str(1:min(18,l)), T{i+1,2:end-1}, round (P(i) * 1e+03));
          else
            fprintf ("%-20s  %10.5g  %6d  %10.5g  %4.3f  %11.2f   1.000 \n", str(1:min(18,l)), T{i+1,2:end-1}); 
          endif
        endfor
        fprintf("Error                 %10.5g  %6d  %10.5g\n", T{end-1,2:4});               
        fprintf("Total                 %10.5g  %6d \n", T{end,2:3});  
        fprintf("\n");
      case "off"
        # do nothing
      otherwise
        error ("anovan: wrong value for ""display"" parameter.");    
    endswitch
  
endfunction


function [X, levels, nlevels, df, termcols] = make_design_matrix (GROUP, TERMS, n, Nm, Nx, Ng)
  
  # Returns a cell array of the design matrix for each term in the model
  
  # Fetch factor levels from each column (i.e. factor) in GROUP
  levels = cell (Nm, 1);
  gid = zeros (n, Nm);
  nlevels = zeros (Nm, 1);
  df = zeros (Nm + Nx, 1);
  termcols = ones (1 + Nm + Nx, 1);
  for j = 1:Nm
    m = find (TERMS(j,:));
    [levels{j}, jnk, gid(:,j)] = unique (GROUP (:,m), "legacy");
    nlevels(j) = numel (levels{j});
    termcols(j+1) = nlevels(j);
    df(j) = nlevels(j) - 1;
  endfor
 
  # Create contrast matrix C and design matrix X
  # Prepare design matrix columns for the main effects
  X = cell (1, 1 + Nm + Nx);
  X(1) = ones (n, 1);
  for j = 1:Nm
    C = contr_sum (nlevels(j));
    func = @(x) x(gid(:,j));
    X(1+j) = cell2mat (cellfun (func, num2cell (C, 1), "UniformOutput", false));
  endfor
  # If applicable, prepare design matrix columns for all the interaction terms
  if (Nx > 0)
    row = TERMS((Ng > 1),:);
    for i = 1:Nx
      I = 1 + find (row(i,:));
      df(Nm+i) = prod (df(I-1));
      termcols(1+Nm+i) = prod (df(I-1) + 1);
      X{1+Nm+i} = X{1};
      for k = 1:numel(I)
        X(1+Nm+i) = bsxfun (@times, X{1+Nm+i}, X{I(k)});
      endfor
    endfor
  endif

endfunction


function C = contr_sum (N)

  # Create contrast matrix (of doubles) using deviation coding 
  # These contrasts sum to 0
  C =  cat (1, eye (N-1), - (ones (1,N-1)));
  
endfunction


function [b, sse, resid] = lmfit (X, Y)
  
  # Get model coefficients by solving the linear equation by QR decomposition 
  # (this achieves the same thing as b = X \ Y)
  # The number of free parameters (i.e. intercept + coefficients) is equal to n - dfe
  [Q, R] = qr (X, 0);
  b = R \ Q' * Y;
 
  # Get fitted values 
  fit = X * b;
  # Get residuals from the fit
  resid = Y - fit;
  # Calculate residual sums-of-squares
  sse = sum ((resid).^2);
  
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
%! # taking the square root of the reported F statistic.
%!
%! score = [4.5 5.6; 3.7 6.4; 5.3 6.4; 5.4 6.0; 3.9 5.7]';
%! treatment = {"before" "after"; "before" "after"; "before" "after"; 
%!              "before" "after"; "before" "after"}';
%! subject = {"GS" "GS"; "JM" "JM"; "HM" "HM"; "JW" "JW"; "PS" "PS"}';
%!
%! [P, T] = anovan (score(:), {treatment(:), subject(:)}, "display", "on",...
%!                  "sstype", 2, "varnames", {"treatment", "subject"});

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
%! # (1994) Psychon Bull Rev. 1(4):476-490, Table 2
%!
%! words = [10 13 13; 6 8 8; 11 14 14; 22 23 25; 16 18 20; ... 
%!          15 17 17; 1 1 4; 12 15 17;  9 12 12;  8 9 12];
%! seconds = [1 2 5; 1 2 5; 1 2 5; 1 2 5; 1 2 5; ...
%!            1 2 5; 1 2 5; 1 2 5; 1 2 5; 1 2 5;];
%! subject = [ 1  1  1;  2  2  2;  3  3  3;  4  4  4;  5  5  5; ...
%!             6  6  6;  7  7  7;  8  8  8;  9  9  9; 10 10 10];
%!
%! [P, T] = anovan (words(:), {seconds(:), subject(:)}, "display", "on",...
%!                  "sstype", 2, "varnames", {"seconds", "subject"});

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
%! # different mouse strains, from Festing (2014), ILAR Journal, 55(3):427-476
%!
%! measurement = [444 614 423 625 408  856 447 719 ...
%!                764 831 586 782 609 1002 606 766]';
%! block = [1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2]';
%! strain= {"NIH","NIH","BALB/C","BALB/C","A/J","A/J","129/Ola","129/Ola",...
%!          "NIH","NIH","BALB/C","BALB/C","A/J","A/J","129/Ola","129/Ola"}';
%! treatment={"C" "T" "C" "T" "C" "T" "C" "T" "C" "T" "C" "T" "C" "T" "C" "T"}';
%!
%! [P, T] = anovan (measurement/10, {block, strain, treatment}, "sstype", 2, ...
%!                  "model", [1 0 0;0 1 0;0 0 1;0 1 1], "display", "on", ...
%!                  "varnames", {"block", "strain", "treatment"});

## Test 1 for anovan example 1 
## Test compares anovan to results from MATLAB's anovan and ttest2 functions
%!test
%! score = [54 23 45 54 45 43 34 65 77 46 65]';
%! gender = {'male' 'male' 'male' 'male' 'male' 'female' 'female' 'female' ...
%!           'female' 'female' 'female'}'; 
%!
%! [P, T] = anovan (score, gender, 'display', 'off', 'varnames', 'gender');
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
%! [P, T] = anovan (score(:), {treatment(:), subject(:)}, 'display', 'off',...
%!                  'sstype', 2, 'varnames', {'treatment', 'subject'});
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
%! [P, T] = anovan (strength, {alloy}, 'display', 'off', 'varnames', {'alloy'});
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
%! [P, T] = anovan (words(:),{seconds(:),subject(:)},'model','linear','sstype',2,'display','off','varnames',{'seconds','subject'});
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
%! [P, T] = anovan (popcorn(:), {brands(:), popper(:)}, 'display', 'off',...
%!                  'model', 'full', 'varnames', {'brands', 'popper'});
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
%! [P, T] = anovan (salary,{gender, degree}, 'model', 'full', 'sstype', 1, 'display','off');
%! assert (P(1), 0.747462549227232,  1e-09);
%! assert (P(2), 1.03809316857694e-08,  1e-09);
%! assert (P(3), 0.523689833702691,  1e-09);
%! assert (T{2,2}, 0.296969696969699,  1e-09);
%! assert (T{3,2}, 272.391841491841,  1e-09);
%! assert (T{4,2}, 1.17482517482512,  1e-09);
%! assert (T{5,2}, 50.0000000000001,  1e-09);
%! [P, T] = anovan (salary,{degree, gender}, 'model', 'full', 'sstype', 1, 'display','off');
%! assert (P(1), 2.53445097305047e-08,  1e-09);
%! assert (P(2), 0.00388133678528749,  1e-09);
%! assert (P(3), 0.523689833702671,  1e-09);
%! assert (T{2,2}, 242.227272727273,  1e-09);
%! assert (T{3,2}, 30.4615384615384,  1e-09);
%! assert (T{4,2}, 1.17482517482523,  1e-09);
%! assert (T{5,2}, 50.0000000000001,  1e-09);
%! [P, T] = anovan (salary,{gender, degree}, 'model', 'full', 'sstype', 2, 'display','off');
%! assert (P(1), 0.00388133678528743,  1e-09);
%! assert (P(2), 1.03809316857694e-08,  1e-09);
%! assert (P(3), 0.523689833702691,  1e-09);
%! assert (T{2,2}, 30.4615384615385,  1e-09);
%! assert (T{3,2}, 272.391841491841,  1e-09);
%! assert (T{4,2}, 1.17482517482512,  1e-09);
%! assert (T{5,2}, 50.0000000000001,  1e-09);
%! [P, T] = anovan (salary,{gender, degree}, 'model', 'full', 'sstype', 3, 'display','off');
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
%! [P, T] = anovan (babble, {sugar, milk}, 'model', 'full', 'sstype', 1, 'display','off');
%! assert (P(1), 0.0108632139833963,  1e-09);
%! assert (P(2), 0.0810606976703546,  1e-09);
%! assert (P(3), 0.00175433329935627,  1e-09);
%! assert (T{2,2}, 3.55752380952381,  1e-09);
%! assert (T{3,2}, 0.956108477471702,  1e-09);
%! assert (T{4,2}, 5.94386771300448,  1e-09);
%! assert (T{5,2}, 3.1625,  1e-09);
%! [P, T] = anovan (babble, {milk, sugar}, 'model', 'full', 'sstype', 1, 'display','off');
%! assert (P(1), 0.0373333189297505,  1e-09);
%! assert (P(2), 0.017075098787169,  1e-09);
%! assert (P(3), 0.00175433329935627,  1e-09);
%! assert (T{2,2}, 1.444,  1e-09);
%! assert (T{3,2}, 3.06963228699552,  1e-09);
%! assert (T{4,2}, 5.94386771300448,  1e-09);
%! assert (T{5,2}, 3.1625,  1e-09);
%! [P, T] = anovan (babble, {sugar, milk}, 'model', 'full', 'sstype', 2, 'display','off');
%! assert (P(1), 0.017075098787169,  1e-09);
%! assert (P(2), 0.0810606976703546,  1e-09);
%! assert (P(3), 0.00175433329935627,  1e-09);
%! assert (T{2,2}, 3.06963228699552,  1e-09);
%! assert (T{3,2}, 0.956108477471702,  1e-09);
%! assert (T{4,2}, 5.94386771300448,  1e-09);
%! assert (T{5,2}, 3.1625,  1e-09);
%! [P, T] = anovan (babble, {sugar, milk}, 'model', 'full', 'sstype', 3, 'display','off');
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
%! [P, T] = anovan (BP(:), {drug(:), feedback(:), diet(:)}, 'model', 'full', 'sstype', 1, 'display','off');
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
%! [P, T] = anovan (BP(:), {drug(:), feedback(:), diet(:)}, 'model', 'full', 'sstype', 2, 'display','off');
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
%! [P, T] = anovan (BP(:), {drug(:), feedback(:), diet(:)}, 'model', 'full', 'sstype', 3, 'display','off');
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
%! [P, T] = anovan (measurement/10,{block,strain,treatment},'model', [1 0 0;0 1 0;0 0 1;0 1 1], 'display','off');
%! assert (P(1), 0.000339814602130731,  1e-09);
%! assert (P(2), 0.0914352969909372,  1e-09);
%! assert (P(3), 5.04077373924908e-05,  1e-09);
%! assert (P(4), 0.0283196918836667,  1e-09);
%! assert (T{2,2}, 1242.5625,  1e-09);
%! assert (T{3,2}, 286.132500000002,  1e-09);
%! assert (T{4,2}, 2275.29,  1e-09);
%! assert (T{5,2}, 495.905000000001,  1e-09);
%! assert (T{6,2}, 207.007499999999,  1e-09);


