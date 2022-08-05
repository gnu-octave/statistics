## Copyright (C) 2022 Andrew Penn <A.C.Penn@sussex.ac.uk>
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
## @deftypefn {Function File} {[@var{P}, @var{T}, @var{STATS}, @var{TERMS}] =} anovan (@var{Y}, @var{GROUP})
## @deftypefnx {Function File} {[@var{P}, @var{T}, @var{STATS}, @var{TERMS}] =}  anovan (@var{Y}, @var{GROUP}, 'name', @var{value})
##
##  Perform a multi-way analysis of variance (ANOVA) with categorical predictors.
## 
##  Data is a single vector @var{Y} with groups specified by a corresponding matrix or 
##  cell array of group labels @var{GROUP}, where each column of @var{GROUP} has the same 
##  number of rows as @var{Y}. For example, if @var{Y} = [1.1;1.2]; @var{GROUP} = [1,2,1; 1,5,2]; 
##  then observation 1.1 was measured under conditions 1,2,1 and observation 1.2 
##  was measured under conditions 1,5,2. Note that groups do not need to be sequentially 
##  numbered.
## 
##  By default, a 'linear' model is used, computing the main effects with no
##  interactions. 
##
##  The settings of anovan can be configured with the following name-value pairs.
##  
##  @var{P} = anovan (@var{Y}, @var{GROUP}, 'model', modeltype)
##  The model to use (modeltype) can specified as one of the following:
##    modeltype = 'linear' to compute n main effects,
##    modeltype = 'interaction' to compute n effects and n*(n-1) two-factor interactions,
##    modeltype = 'full' to compute the n main effects and interactions at all levels,
##    an integer representing the maximum interaction order, or
##    a matrix of term definitions, where each row is a term and each column is a factor.
##    For example, a two-way ANOVA with interaction would be: [1 0; 0 1; 1 1]
##
##  @var{P} = anovan (@var{Y}, @var{GROUP}, 'sstype', sstype)
##  The type of sum-of-squares 1, 2 or 3 (default = 3)
## 
##  @var{P} = anovan (@var{Y}, @var{GROUP}, 'varnames', varnames)
##  - optionally, a factor name for each column of GROUP can be provided in the 
##    input argument. varnames should be a cell array of strings. By default, 
##    varnames are 'X1','X2','X3', etc. 
##  
##  @var{P} = anovan (@var{Y}, @var{GROUP}, 'display', 'on')
##  - 'on' (default) | 'off': switch display of ANOVA table on/off
##
##  [@var{P}, @var{T}] = anovan (...) returns a cell array containing the ANOVA table
##
##  [@var{P}, @var{T}, @var{STATS}] = anovan (...) returns a structure containing additional
##  statistics, including coefficients of the linear model, the model residuals, 
##  and the number of levels in each factor.
## 
##  [@var{P}, @var{T}, @var{STATS}, @var{TERMS}] = anovan (...) returns the model term definitions
##
## @end deftypefn

##  Author: Andrew Penn <a.c.penn@sussex.ac.uk>
##  Includes some code by: Andy Adler <adler@ncf.ca> and Christian Scholz
## 

function [P, T, STATS, TERMS] = anovan (Y, GROUP, varargin)
      
    if nargin <= 1
      error ('anovan usage: ''anovan (Y, GROUP)''; atleast 2 input arguments required');
    endif

    # Check supplied parameters
    modeltype = 'linear';
    display = 'on';
    sstype = 3;
    varnames = [];
    for idx = 3:2:nargin
      name = varargin{idx-2};
      value = varargin{idx-1};
      if strcmpi (name, 'model')
        modeltype = value;
      elseif strcmpi (name, 'varnames')
        varnames = value;
      elseif strcmpi (name, 'display')
        display = value;   
      elseif strcmpi (name, 'sstype') 
        sstype = value;
      else 
        error (sprintf('anovan: parameter %s is not supported', name));
      endif
    endfor
    
    # Remove NaN or non-finite observations
    excl = logical (isnan(Y) + isinf(Y));
    Y(excl) = [];
    GROUP(excl,:) = [];
    n = numel (Y);
    if (prod (size (Y)) ~= n)
      error ('anovan: for ''anovan (Y, GROUP)'', Y must be a vector');
    endif
    if (size (Y, 2) > 1)
      Y = Y(:);
    endif
    N = size (GROUP,2); # number of anova "ways"
    # Accomodate for different formats for GROUP 
    # GROUP can be a matrix of numeric identifiers of a cell arrays of strings or numeric idenitiers
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
      error ('anovan: GROUP must be a matrix of the same number of rows as Y');
    endif
    if ~isempty (varnames) 
      if iscell (varnames)
        if all (cellfun (@ischar, varnames))
          nvarnames = numel(varnames);
        else
          error ('anovan: all variable names must be character or character arrays');
        endif
      elseif ischar (varnames)
        nvarnames = 1;
        varnames = {varnames};
      elseif isstring (varnames)
        nvarnames = 1;
        varnames = {char(varnames)};
      else
        error ('anovan: varnames is not of a valid type. Must be cell array of character arrays, character array or string');
      endif
    else
      nvarnames = N;
      varnames = arrayfun(@(x) ['X',num2str(x)], 1:N, 'UniformOutput', 0);
    endif
    if (nvarnames ~= N)
      error ('anovan: number of variable names is not equal to number of grouping variables');
    endif

    # Evaluate model type input argument and create terms matrix if not provided
    if ischar (modeltype)
      switch lower(modeltype)
        case 'linear'
          modeltype = 1;
        case 'interaction'
          modeltype = 2;
        case 'full'
          modeltype = N;
      endswitch
    endif
    if isscalar (modeltype)
      TERMS = cell (modeltype,1);
      v = false (1,N);
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
          # Create term definitions for a full model
          Nx = zeros (1, N-1);
          Nx = 0;
          for k = 1:N
            Nx = Nx + nchoosek(N,k);
          endfor
          for j = 1:modeltype
            v(1:j) = 1;
            TERMS(j) = flipud (unique (perms (v), 'rows'));
          endfor
          TERMS = cell2mat (TERMS);
      endswitch
    else
      # Assume that the user provided a suitable matrix of term definitions
      if (size (modeltype, 2) > N)
        error ('anovan: the number of columns in the term definitions cannot exceed the number of columns of GROUP')
      endif
      TERMS = logical (modeltype);
    endif
    # Evaluate terms matrix
    Ng = sum (TERMS, 2); 
    if any (diff (Ng) < 0)
      error ('anovan: the model terms matrix must list main effects above interactions')
    endif
    Nm = sum (Ng == 1);
    Nx = sum (Ng > 1);
    Nt = Nm + Nx;
    if any (any (TERMS(1:Nm,:)) ~= any (TERMS))
      error ('anovan: all factors involved in interactions must have a main effect')  
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
        sstype_char = 'I';
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
        sstype_char = 'II';
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
        sstype_char = 'III';
      otherwise
        # sstype 'h' not supported
        error ('anovan: only sstype 1, 2 and 3 are supported')
    endswitch
    dfe = dft - sum (df);
    ms = ss ./ df;
    mse = sse / dfe;
    F = ms / mse;
    P = 1 - fcdf (F, df, dfe);

    # Prepare stats output structure
    # Note that the information provided by STATS is not sufficient for MATLAB's multcompare function
    STATS = struct ('source','anovan', ...
                    'resid', resid, ...
                    'coeffs', b, ...
                    'Rtr', [], ...           # Not used in Octave
                    'rowbasis', [], ...      # Not used in Octave
                    'dfe', dfe, ...
                    'mse', mse, ...
                    'nullproject', [], ...   # Not used in Octave
                    'terms', TERMS, ...
                    'nlevels', nlevels, ...  
                    'continuous', [], ...
                    'vmeans', [], ...        # Not used since 'continuous' argument name not supported
                    'termcols', termcols, ...
                    'coeffnames', [], ...    # Not used in Octave
                    'vars', [], ...          # Not used in Octave
                    'varnames', {varnames}, ...
                    'grpnames', {grpnames}, ...
                    'vnested', [], ...       # Not used since 'nested' argument name not supported
                    'ems', [], ...           # Not used since 'nested' argument name not supported
                    'denom', [], ...         # Not used since 'random' argument name not supported
                    'dfdenom', [], ...       # Not used since 'random' argument name not supported
                    'msdenom', [], ...       # Not used since 'random' argument name not supported
                    'varest', [], ...        # Not used since 'random' argument name not supported
                    'varci', [], ...         # Not used since 'random' argument name not supported
                    'txtdenom', [], ...      # Not used since 'random' argument name not supported
                    'txtems', [], ...        # Not used since 'random' argument name not supported
                    'rtnames', []);          # Not used since 'random' argument name not supported
    
    # Prepare cell array containing the ANOVA table (T)
    T = cell (Nt + 3, 6);
    T(1,:) = {'Source','Sum Sq.','d.f.','Mean Sq.','F','Prob>F'};
    T(2:Nt+1,2:6) = num2cell([ss df ms F P]);
    T(end-1,1:4) = {'Error',sse,dfe,mse};
    T(end,1:3) = {'Total',sst,dft};
    for i = 1:Nt
      str = sprintf('%s*',varnames{find(TERMS(i,:))});
      T(i+1,1) = str(1:end-1);
    endfor
    
    # Print ANOVA table 
    if strcmpi(display,'on')
      # Get dimensions of the ANOVA table
      [nrows, ncols] = size (T);
      # Print table
      fprintf('\n%d-way ANOVA table (Type %s sums of squares):\n\n', Nm, sstype_char);
      fprintf('Source                    Sum Sq.    d.f.    Mean Sq.           F  Prob>F\n');
      fprintf('*************************************************************************\n');  
      for i = 1:Nt
        str = T{i+1,1};
        l = numel(str);  # Needed to truncate source term name at 21 characters
        # Format and print the statistics for each model term
        # Format F statistics and p-values in APA style
        if (P(i) < 0.001)
          fprintf ('%-21s  %10.5g  %6d  %10.5g %11.2f   <.001 \n', str(1:min(21,l)), T{i+1,2:end-1});
        elseif (P(i) < 1.0)
          fprintf ('%-21s  %10.5g  %6d  %10.5g %11.2f    .%03u \n', str(1:min(21,l)), T{i+1,2:end-1}, round (P(i) * 1e+03));
        else
          fprintf ('%-21s  %10.5g  %6d  %10.5g %11.2f   1.000 \n', str(1:min(21,l)), T{i+1,2:end-1}); 
        endif
      endfor
      fprintf('Error                  %10.5g  %6d  %10.5g\n', T{end-1,2:4});               
      fprintf('Total                  %10.5g  %6d \n', T{end,2:3});  
      fprintf('\n');
    elseif strcmp(display,'off')
      # do nothing
    else
      error ('anovan: unknown display option');    
    endif
  
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
    [levels{j}, jnk, gid(:,j)] = unique (GROUP (:,m), 'legacy');
    nlevels(j) = numel (levels{j});
    termcols(j+1) = nlevels(j);
    df(j) = nlevels(j) - 1;
  endfor
 
  # Create contrast matrix C and dummy variables X
  # Prepare design matrix columns for the main effects
  X = cell (1, 1 + Nm + Nx);
  X(1) = ones (n, 1);
  for j = 1:Nm
    C = contr_sum (nlevels(j));
    func = @(x) x(gid(:,j));
    X(1+j) = cell2mat (cellfun (func, num2cell (C, 1), 'UniformOutput', false));
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


## Test for unbalanced two-way ANOVA (2x2) from Maxwell, Delaney and Kelly (2018): Chapter 7, Table 15)
## https://designingexperiments.com/csv-chapter-data/
## Test compares to results in matlab 
%!test
%! salary = [24 26 25 24 27 24 27 23 15 17 20 16, ...
%!           25 29 27 19 18 21 20 21 22 19]';
%! gender = {'f' 'f' 'f' 'f' 'f' 'f' 'f' 'f' 'f' 'f' 'f' 'f'...
%!           'm' 'm' 'm' 'm' 'm' 'm' 'm' 'm' 'm' 'm'}';
%! degree = [1 1 1 1 1 1 1 1 0 0 0 0 1 1 1 0 0 0 0 0 0 0]';
%! [P, T] = anovan (salary,{gender, degree}, 'model', 'full', 'sstype', 1, 'display','off');
%! assert (P(1), 0.747462549227232, 1e-12);
%! assert (P(2), 1.03809316857694e-08, 1e-12);
%! assert (P(3), 0.523689833702691, 1e-12);
%! assert (T{2,2}, 0.296969696969699, 1e-12);
%! assert (T{3,2}, 272.391841491841, 1e-12);
%! assert (T{4,2}, 1.17482517482512, 1e-12);
%! assert (T{5,2}, 50.0000000000001, 1e-12);
%! [P, T] = anovan (salary,{degree, gender}, 'model', 'full', 'sstype', 1, 'display','off');
%! assert (P(1), 2.53445097305047e-08, 1e-12);
%! assert (P(2), 0.00388133678528749, 1e-12);
%! assert (P(3), 0.523689833702671, 1e-12);
%! assert (T{2,2}, 242.227272727273, 1e-12);
%! assert (T{3,2}, 30.4615384615384, 1e-12);
%! assert (T{4,2}, 1.17482517482523, 1e-12);
%! assert (T{5,2}, 50.0000000000001, 1e-12);
%! [P, T] = anovan (salary,{gender, degree}, 'model', 'full', 'sstype', 2, 'display','off');
%! assert (P(1), 0.00388133678528743, 1e-12);
%! assert (P(2), 1.03809316857694e-08, 1e-12);
%! assert (P(3), 0.523689833702691, 1e-12);
%! assert (T{2,2}, 30.4615384615385, 1e-12);
%! assert (T{3,2}, 272.391841491841, 1e-12);
%! assert (T{4,2}, 1.17482517482512, 1e-12);
%! assert (T{5,2}, 50.0000000000001, 1e-12);
%! [P, T] = anovan (salary,{gender, degree}, 'model', 'full', 'sstype', 3, 'display','off');
%! assert (P(1), 0.00442898146583742, 1e-12);
%! assert (P(2), 1.30634252053587e-08, 1e-12);
%! assert (P(3), 0.523689833702691, 1e-12);
%! assert (T{2,2}, 29.3706293706294, 1e-12);
%! assert (T{3,2}, 264.335664335664, 1e-12);
%! assert (T{4,2}, 1.17482517482512, 1e-12);
%! assert (T{5,2}, 50.0000000000001, 1e-12);


