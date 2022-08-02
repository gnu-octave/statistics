## Copyright (C) 2022 Andrew Penn <A.C.Penn@sussex.ac.uk>
## Copyright (C) 2022 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## Perform a multi-way analysis of variance (ANOVA) with categorical predictors.
## When interpreting the results from an analysis of an unbalanced design, be  
## aware that this function calculates sum-of-squared residuals sequentially 
## from means that are weighted by sample size (i.e. 'sstype' = 1) and that
## the order of the factors (columns in @var{GROUP}) matters.
## 
## Data is a single vector @var{Y} with groups specified by a corresponding
## matrix or cell array of group labels @var{GROUP}, where each column of
## @var{GROUP} has the same number of rows as @var{Y}.  For example, if 
## @code{@var{Y} = [1.1;1.2]} and @code{@var{GROUP} = [1,2,1; 1,5,2]}, then
## observation 1.1 was measured under conditions 1,2,1 and observation 1.2 was
## measured under conditions 1,5,2.  Note that groups do not need to be
## sequentially numbered.
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
## @code{[@dots{}] = anovan (@var{Y}, @var{GROUP}, "varnames", @var{varnames})}
##
## @itemize
## @item
## @var{varnames} must be a cell array of strings with each element containing a
## factor name for each column of GROUP.  By default (if not parsed as optional
## argument), @var{varnames} are 'X1','X2','X3', etc. 
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
## @end deftypefn

##  Author: Andrew Penn <a.c.penn@sussex.ac.uk>
##  Includes some code by: Andy Adler <adler@ncf.ca> and Christian Scholz

function [P, T, STATS, TERMS] = anovan (Y, GROUP, varargin)
      
  if nargin < 2
    error ("anovan: at least 2 input arguments required.");
  endif

  ## Check supplied parameters
  modeltype = "linear";
  display = "on";
  sstype = 1; # default: this is not optional (yet)
  varnames = [];
  for idx = 3:2:nargin
    name = varargin{idx-2};
    value = varargin{idx-1};
    if strcmpi (name, "model")
      modeltype = value;
    elseif strcmpi (name, "varnames")
      varnames = value;
    elseif strcmpi (name, "display")
      if strcmpi (value, "on") || strcmpi (value, "off")
        display = value;
      else
        error ("anovan: wrong value for 'display' parameter.");
      endif
    elseif strcmpi (name, "sstype") 
      sstype = value;
    else 
      error (sprintf ("anovan: parameter %s is not supported.", name));
    endif
  endfor
  
  ## Remove NaN or non-finite observations
  excl = logical (isnan (Y) + isinf (Y));
  Y(excl) = [];
  GROUP(excl,:) = [];
  N = numel (Y);
  if prod (size (Y)) != N
    error ("anovan: Y must be a vector.");
  endif
  if (size (Y, 2) > 1)
    Y = Y(:);
  endif
  n = size (GROUP, 2); # number of anova "ways"
  ## Accomodate for different formats for GROUP 
  ## GROUP can be a matrix of numeric identifiers of a cell arrays of strings
  ## or numeric idenitiers
  if iscell(GROUP)
    if (size(GROUP, 1) == 1)
      for k = 1:n
        if isnumeric (GROUP{k})
          tmp(:,k) = cellstr (num2str (GROUP{k}));
        else
          tmp(:,k) = GROUP{k};
        endif
      endfor
      GROUP = tmp;
    else
      for k = 1:n
        tmp(:,k) = cellstr (char (GROUP{k}));
      endfor
    endif
  endif
  if (size (GROUP,1) != N)
    error ("anovan: GROUP must be a matrix of the same number of rows as Y.");
  end
  if ~isempty (varnames) 
    if iscell (varnames)
      if all (cellfun (@ischar, varnames))
        ng = numel (varnames);
      else
        error ("anovan: variable names must be character or character arrays.");
      endif
    elseif ischar (varnames)
      ng = 1;
      varnames = {varnames};
    elseif isstring (varnames)
      ng = 1;
      varnames = {char(varnames)};
    else
      error ("anovan: varnames is not of a valid type.");
    endif
  else
    ng = n;
    varnames = arrayfun(@(x) ['X',num2str(x)], 1:n, "UniformOutput", 0);
  endif
  if (ng ~= n)
    error ("anovan: variable names is not equal to groups.");
  endif

  ## Evaluate model type input argument and create terms matrix if not provided
  msg1 = "anovan: interactions of more than two factors are not supported.";
  msg2 = "anovan: model terms matrix must list main effects above interactions.";
  if (isnumeric (modeltype) && (numel (modeltype) == 1))
    if (modeltype == 1)
      modeltype = "linear";
    elseif (modeltype == 2)
      modeltype = "interaction";
    else
      error (msg);
    endif
  endif
  if ! isnumeric (modeltype)
    if strcmpi (modeltype, "full") && (n < 3)
      modeltype = "interaction";
    endif
    switch lower(modeltype)
      case "linear"
        nx = 0;
        TERMS = zeros (n + nx, n);
        TERMS(1:n,:) = diag (ones (n, 1));
      case "interaction"
        nx = nchoosek (n, 2);
        TERMS = zeros (n + nx, n);
        TERMS(1:n,:) = diag (ones (n, 1));
        for j = 1:n
          for i = j:n-1
            TERMS(n+j+i-1,j) = 1;
            TERMS(n+j+i-1,i+1) = 1;
          endfor
        endfor
      case "full"
        error (msg);
    endswitch
  else
    ## Assume that the user provided a terms matrix
    TERMS = modeltype;
  endif
  TERMS = logical(TERMS);
  ## Evaluate terms matrix
  ng = sum (TERMS, 2); 
  if any(diff(ng) < 0)
    error (msg2);
  endif
  main = (ng == 1);
  inte = (ng == 2);
  if any (ng > 2)
    error (msg1);
  endif
  nm = sum (main);
  nx = sum (inte);
  nt = nm + nx;

  ## Calculate total sum-of-squares
  ct  = sum (Y)^2 / N;   % correction term
  sst = sum (Y.^2) - ct;
  dft = N - 1;

  ## Fit linear models, and calculate sums-of-squares for ANOVA
  switch lower (sstype)
    case 1
      ## Type I sequential sums-of-squares (sstype = 1)
      R = sst;
      ss = zeros (nt,1);
      df = zeros (nt,1);
      [X, grpnames, nlevels, df, termcols] = mk_design_matrix (GROUP, ...
                                             TERMS, main, inte, N, nm, nx);
      for j = 1:nt
        XS = cell2mat (X(1:j+1));
        [b, sse, resid] = lmfit (XS, Y);
        ss(j) = R - sse;
        R = sse;
      endfor
    otherwise
      ## sstype 2, 3, or 'h' not supported (yet)
      error ("anovan: only sstype = 1 is currently supported.")
  endswitch
  dfe = dft - sum (df);
  ms = ss ./ df;
  mse = sse / dfe;
  F = ms / mse;
  P = 1 - fcdf (F, df, dfe);

  ## Prepare stats output structure
  ## Note that the information provided by STATS is not sufficient
  ## for MATLAB's multcompare function
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
                  'vmeans', [], ...        # 'continuous' argument not supported
                  'termcols', termcols, ...
                  'coeffnames', [], ...    # Not used in Octave
                  'vars', [], ...          # Not used in Octave
                  'varnames', {varnames}, ...
                  'grpnames', {grpnames}, ...
                  'vnested', [], ...       # 'nested' argument not supported
                  'ems', [], ...           # 'nested' argument not supported
                  'denom', [], ...         # 'random' argument not supported
                  'dfdenom', [], ...       # 'random' argument not supported
                  'msdenom', [], ...       # 'random' argument not supported
                  'varest', [], ...        # 'random' argument not supported
                  'varci', [], ...         # 'random' argument not supported
                  'txtdenom', [], ...      # 'random' argument not supported
                  'txtems', [], ...        # 'random' argument not supported
                  'rtnames', []);          # 'random' argument not supported
  
  ## Prepare cell array containing the ANOVA table (T)
  T = cell (nt+3, 6);
  T(1,:) = {'Source','Sum Sq.','d.f.','Mean Sq.','F','Prob>F'};
  T(2:nt+1,2:6) = num2cell([ss df ms F P]);
  T(end-1,1:4) = {'Error',sse,dfe,mse};
  T(end,1:3) = {'Total',sst,dft};
  for i=1:nt
    str = sprintf('%s*',varnames{find(TERMS(i,:))});
    T(i+1,1) = str(1:end-1);
  endfor
  
  ## Print ANOVA table 
  if strcmpi (display, "on")
    ## Get dimensions of the ANOVA table
    [nrows, ncols] = size (T);
    ## Print table
    fprintf("\n%d-way ANOVA Table:\n\n", nm);
    fprintf("Source         Sum Sq.      df     Mean Sq.       F     Prob>F\n");
    fprintf("**************************************************************\n");  
    for i = 1:nt
      str = T{i+1,1};
      l = numel(str);  # Needed to truncate source term name at 10 characters
      ## Format and print tehe statistics for each model term
      ## Format F statistics and p-values in APA style
      if (P(i) < 0.001)
        fprintf ("%-10s  %10.5g  %6d  %10.5g %10.2f    <.001 \n", ...
                 str(1:min(10,l)), T{i+1,2:end-1});
      elseif (P(i) < 1.0)
        fprintf ("%-10s  %10.5g  %6d  %10.5g %10.2f     .%03u \n", ...
                 str(1:min(10,l)), T{i+1,2:end-1}, round (P(i) * 1e+03));
      else
        fprintf ("%-10s  %10.5g  %6d  %10.5g %10.2f    1.000 \n", ...
                 str(1:min(10,l)), T{i+1,2:end-1}); 
      endif
    endfor
    fprintf("Error       %10.5g  %6d  %10.5g\n", T{end-1,2:4});               
    fprintf("Total       %10.5g  %6d \n", T{end,2:3});  
    fprintf("\n"); 
  endif
endfunction


function [X, levels, nlevels, df, termcols] = mk_design_matrix (GROUP, ...
                                              TERMS, main, inte, N, nm, nx)
  ## Returns a cell array of dummy-coded levels for each term in the linear model
  ## Fetch factor levels from each column (i.e. factor) in GROUP
  levels = cell (nm, 1);
  gid = zeros (N, nm);
  nlevels = zeros (nm, 1);
  df = zeros (nm + nx, 1);
  termcols = ones (1 + nm + nx, 1);
  for j = 1:nm
    m = find (TERMS(j,:));
    [levels{j}, jnk, gid(:,j)] = unique (GROUP (:,m), 'legacy');
    nlevels(j) = numel (levels{j});
    termcols(j+1) = nlevels(j);
    df(j) = nlevels(j) - 1;
  endfor
  ## Create contrast matrix C and dummy variables X
  ## Prepare dummy variables for main effects
  X = cell (1, 1 + nm + nx);
  X(1) = ones (N, 1);
  for j = 1:nm
    C = contr_sum (nlevels(j));
    func = @(x) x(gid(:,j));
    X(1+j) = cell2mat (cellfun (func, num2cell (C, 1), 'UniformOutput', false));
  endfor
  ## If applicable, prepare dummy variables for all two-factor interactions
  if (nx > 0)
    pairs = TERMS(inte,:);
    for i = 1:nx
      I = 1 + find (pairs(i,:));
      X(1+nm+i) = bsxfun (@times, X{I});
      df(nm+i) = prod (df(I-1));
      termcols(1+nm+i) = prod (df(I-1) + 1);
    endfor
  endif
endfunction


function C = contr_sum (n)
  ## Create contrast matrix (of doubles) using deviation coding 
  ## These contrasts sum to 0
  C =  cat (1, diag (ones (n-1, 1)), - (ones (1,n-1)));
endfunction


function [b, sse, resid] = lmfit (X, Y)
  ## Get model coefficients by solving the linear equation by QR decomposition 
  ## (this achieves the same thing as b = X \ Y)
  ## The number of free parameters (i.e. intercept + coefficients)
  ## is equal to N - dfe
  [Q, R] = qr (X, 0);
  b = R \ Q' * Y;
  ## Get fitted values 
  fit = X * b;
  ## Get residuals from the fit
  resid = Y - fit;
  ## Calculate residual sums-of-squares
  sse = sum ((resid).^2);
endfunction

%!demo
%! salary = [24 26 25 24 27 24 27 23 15 17 20 16, ...
%!           25 29 27 19 18 21 20 21 22 19]';
%! gender = {"f" "f" "f" "f" "f" "f" "f" "f" "f" "f" "f" "f"...
%!           "m" "m" "m" "m" "m" "m" "m" "m" "m" "m"}';
%! degree = [1 1 1 1 1 1 1 1 0 0 0 0 1 1 1 0 0 0 0 0 0 0]';
%! [P, T] = anovan (salary,{gender, degree}, "model", "interaction", ...
%!                  "sstype", 1, "display", "on");

## Test for unbalanced two-way ANOVA (2x2) from Maxwell, Delaney and Kelly (2018
## Chapter 7, Table 15) https://designingexperiments.com/csv-chapter-data/
## Test compares to results in matlab 
%!test
%! salary = [24 26 25 24 27 24 27 23 15 17 20 16, ...
%!           25 29 27 19 18 21 20 21 22 19]';
%! gender = {"f" "f" "f" "f" "f" "f" "f" "f" "f" "f" "f" "f"...
%!           "m" "m" "m" "m" "m" "m" "m" "m" "m" "m"}';
%! degree = [1 1 1 1 1 1 1 1 0 0 0 0 1 1 1 0 0 0 0 0 0 0]';
%! [P, T] = anovan (salary,{gender, degree}, "model", "interaction", ...
%!                  "sstype", 1, "display", "off");
%! assert (P(1), 0.747462549227232, 1e-12);
%! assert (P(2), 1.03809316857694e-08, 1e-12);
%! assert (P(3), 0.523689833702691, 1e-12);
%! assert (T{2,2}, 0.296969696969699, 1e-12);
%! assert (T{3,2}, 272.391841491841, 1e-12);
%! assert (T{4,2}, 1.17482517482512, 1e-12);
%! [P, T] = anovan (salary,{degree, gender}, "model", "interaction", ...
%!                  "sstype", 1, "display", "off");
%! assert (P(1), 2.53445097305047e-08, 1e-12);
%! assert (P(2), 0.00388133678528749, 1e-12);
%! assert (P(3), 0.523689833702671, 1e-12);
%! assert (T{2,2}, 242.227272727273, 1e-12);
%! assert (T{3,2}, 30.4615384615384, 1e-12);
%! assert (T{4,2}, 1.17482517482523, 1e-12);



