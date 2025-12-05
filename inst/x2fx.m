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
## @deftypefn  {statistics} {[@var{d}, @var{model}, @var{termstart}, @var{termend}] =} x2fx (@var{x})
## @deftypefnx {statistics} {[@var{d}, @var{model}, @var{termstart}, @var{termend}] =} x2fx (@var{x}, @var{model})
## @deftypefnx {statistics} {[@var{d}, @var{model}, @var{termstart}, @var{termend}] =} x2fx (@var{x}, @var{model}, @var{categ})
## @deftypefnx {statistics} {[@var{d}, @var{model}, @var{termstart}, @var{termend}] =} x2fx (@var{x}, @var{model}, @var{categ}, @var{catlevels})
##
## Convert predictors to design matrix.
##
## @code{@var{d} = x2fx (@var{x}, @var{model})} converts a matrix of predictors
## @var{x} to a design matrix @var{d} for regression analysis.  Distinct
## predictor variables should appear in different columns of @var{x}.  The input
## @var{x} can be a numeric matrix, or a mixed-type array including categorical,
## string, character, and numeric columns.
##
## The optional input @var{model} controls the regression model.  By default,
## @code{x2fx} returns the design matrix for a linear additive model with a
## constant term.  @var{model} can be any one of the following strings:
##
## @multitable @columnfractions 0.05 0.2 0.75
## @item @tab "linear" @tab Constant and linear terms (the default)
## @item @tab "interaction" @tab Constant, linear, and interaction terms
## @item @tab "quadratic" @tab Constant, linear, interaction, and squared terms
## @item @tab "purequadratic" @tab Constant, linear, and squared terms
## @end multitable
##
## If @var{x} has n columns, the order of the columns of @var{d} for a full
## quadratic model is:
##
## @itemize
## @item
## The constant term.
## @item
## The linear terms (the columns of X, in order 1,2,...,n).
## @item
## The interaction terms (pairwise products of columns of @var{x}, in order
## (1,2), (1,3), ..., (1,n), (2,3), ..., (n-1,n).
## @item
## The squared terms (in the order 1,2,...,n).
## @end itemize
##
## Other models use a subset of these terms, in the same order.
##
## Alternatively, MODEL can be a matrix specifying polynomial terms of arbitrary
## order.  In this case, MODEL should have one column for each column in X and
## one row for each term in the model.  The entries in any row of MODEL are powers
## for the corresponding columns of @var{x}.  For example, if @var{x} has
## columns X1, X2, and X3, then a row [0 1 2] in @var{model} would specify the
## term (X1.^0).*(X2.^1).*(X3.^2).  A row of all zeros in @var{model} specifies
## a constant term, which you can omit.
##
## MODEL can also be a function handle or a string naming a function that
## will be applied to @var{x} to generate the design matrix.
##
## @code{@var{d} = x2fx (@var{x}, @var{model}, @var{categ})} treats columns with
## numbers listed in the vector @var{categ} as categorical variables.  Terms
## involving categorical variables produce dummy variable columns in @var{d}.
## Dummy variables are computed under the assumption that possible categorical
## levels are completely enumerated by the unique values that appear in the
## corresponding column of @var{x}.
##
## @code{@var{d} = x2fx (@var{x}, @var{model}, @var{categ}, @var{catlevels})}
## accepts a vector @var{catlevels} the same length as @var{categ}, specifying
## the number of levels in each categorical variable.  In this case, values in
## the corresponding column of @var{x} must be integers in the range from 1 to
## the specified number of levels.  Not all of the levels need to appear in
## @var{x}.
##
## @end deftypefn

function [D, model, termstart, termend] = x2fx (x, model, categ, catlevels)
  original_x = x;
  [m, n] = size (x);
  ## Determine data class for output
  if (isa (x, 'single'))
    data_class = 'single';
  else
    data_class = 'double';
  endif

  auto_detected_categ = [];
  x_converted = zeros (m, n, data_class);
  for col_idx = 1:n
    col = x(:, col_idx);
    ## Check for categorical type
    if (isa (col, 'categorical'))
      x_converted(:, col_idx) = double (col);
      auto_detected_categ = [auto_detected_categ, col_idx];
    ## Check for string type
    elseif (isa (col, 'string'))
      [~, ~, idx] = unique (col);
      x_converted(:, col_idx) = double (idx);
      auto_detected_categ = [auto_detected_categ, col_idx];
    ## Check for cellstr
    elseif (iscellstr (col))
      [~, ~, idx] = unique (col);
      x_converted(:, col_idx) = double (idx);
      auto_detected_categ = [auto_detected_categ, col_idx];
    ## Check for char matrix
    elseif (ischar (col) && ndims (col) == 2)
      col_cellstr = cellstr (col);
      [~, ~, idx] = unique (col_cellstr);
      x_converted(:, col_idx) = double (idx);
      auto_detected_categ = [auto_detected_categ, col_idx];
    ## Numeric or logical columns
    elseif (isnumeric (col) || islogical (col))
      x_converted(:, col_idx) = double (col);
    else
      error ("x2fx: unsupported data type in column %d.", col_idx);
    endif
  endfor
  x = x_converted;

  ## Check for input arguments
  if (nargin < 2 || isempty (model))
    model = 'linear';
  endif
  if (nargin < 3)
    categ = [];
  else
    if (! all (ismember (categ, 1:n)))
      error ("x2fx: category index exceeds number of columns in X.");
    endif
  endif

  if (! isempty (auto_detected_categ))
    categ = unique ([categ(:)', auto_detected_categ]);
  endif

  if (nargin < 4)
    catlevels = [];
  endif

  ## Handle function handles or string function names
  if (isa (model, 'function_handle'))
    D = feval (model, original_x);
    termstart = [];
    termend = [];
    return
  endif

  ## Convert model string to model matrix
  if (ischar (model) || isstring (model))
    model = char (model);
    if (strcmpi (model, "linear") || strcmpi (model, "additive"))
      interactions = false;
      quadratic = false;
    elseif (strcmpi (model, "interaction"))
      interactions = true;
      quadratic = false;
    elseif (strcmpi (model, "quadratic"))
      interactions = true;
      quadratic = true;
    elseif (strcmpi (model, "purequadratic"))
      interactions = false;
      quadratic = true;
    else
      try
        D = feval (model, original_x);
      catch
        error ("x2fx: unrecognized function '%s'.", model);
      end_try_catch
      termstart = [];
      termend = [];
      return
    endif
    I = eye (n);
    ## Construct interactions part
    if (interactions && n > 1)
      [r, c] = find (tril (ones (n), -1));
      nt = length (r);
      intpart = zeros (nt, n);
      intpart(sub2ind (size (intpart), (1:nt)', r)) = 1;
      intpart(sub2ind (size (intpart), (1:nt)', c)) = 1;
    else
      intpart = zeros (0, n);
    endif
    ## Construct quadratic part
    if (quadratic)
      quadpart = 2 * I;
      quadpart(categ,:) = [];
    else
      quadpart = zeros (0, n);
    endif
    model = [zeros(1, n); I];
    model = [model; intpart; quadpart];
  endif

  ## Process categorical variables and convert to numeric indices
  catmember = ismember (1:n, categ);
  var_DF = ones (1, n);
  if (isempty (catlevels))
    ## Get unique values for each categorical variable
    for idx = 1:length (categ)
      categ_idx = categ(idx);
      [unique_vals, ~, J] = unique (x(:,categ_idx), 'rows');
      var_DF(categ_idx) = length (unique_vals) - 1;
      x(:,categ_idx) = J;
    endfor
  else
    ## Use specified category levels
    var_DF(categ) = catlevels - 1;
    for idx = 1:length (categ)
      categ_idx = categ(idx);
      if (any (! ismember (x(:,categ_idx), 1:catlevels(idx))))
        error ("x2fx: wrong value %f in category %d.", ...
               catlevels(idx), categ_idx);
      endif
    endfor
  endif

  ## Get size of model matrix
  [r, c] = size (model);
  ## Check for equal number of columns between x and model
  if (c != n)
    error ("x2fx: wrong number of columns between X and MODEL.");
  endif

  ## Allocate space for the design matrix
  termdf = prod (max (1, (model > 0) .* repmat (var_DF, r, 1)), 2);
  termend = cumsum (termdf);
  termstart = termend - termdf + 1;
  D = zeros (m, termend(end), data_class);
  allrows = (1:m)';
  for idx = 1:r
    cols = termstart(idx):termend(idx);
    pwrs = model(idx,:);
    t = pwrs > 0;
    C = 1;
    if (any (t))
      if (any (pwrs(! catmember)))
        pwrs_cont = pwrs .* (! catmember);
        C = ones (size (x, 1), 1);
        collist = find (pwrs_cont > 0);
        for j = 1:length (collist)
          col_idx = collist(j);
          C = C .* (x(:,col_idx) .^ pwrs_cont(col_idx));
        endfor
      endif
      if (any (pwrs(catmember) > 0))
        Z = zeros (m, termdf(idx));
        collist = find (pwrs > 0 & catmember);
        xcol = x(:,collist(1));
        keep = (xcol <= var_DF(collist(1)));
        colnum = xcol;
        cumdf = 1;
        for j = 2:length (collist)
          cumdf = cumdf * var_DF(collist(j-1));
          xcol = x(:,collist(j));
          keep = keep & (xcol <= var_DF(collist(j)));
          colnum = colnum + cumdf * (xcol - 1);
        endfor
        if (length (C) > 1)
          C = C(keep);
        endif
        Z(sub2ind(size(Z), allrows(keep), colnum(keep))) = C;
        C = Z;
      endif
    endif
    D(:,cols) = C;
  endfor
endfunction

%!test
%! X = [1, 10; 2, 20; 3, 10; 4, 20; 5, 15; 6, 15];
%! D = x2fx(X,'quadratic');
%! assert (D(1,:), [1, 1, 10, 10, 1, 100]);
%! assert (D(2,:), [1, 2, 20, 40, 4, 400]);

%!test
%! X = [1, 10; 2, 20; 3, 10; 4, 20; 5, 15; 6, 15];
%! model = [0, 0; 1, 0; 0, 1; 1, 1; 2, 0];
%! D = x2fx(X,model);
%! assert (D(1,:), [1, 1, 10, 10, 1]);
%! assert (D(2,:), [1, 2, 20, 40, 4]);
%! assert (D(4,:), [1, 4, 20, 80, 16]);

%!test
%! x = [1, 2, 3; 2, 3, 4; 3, 4, 5];
%! D = x2fx (x, 'linear');
%! assert (D, [1, 1, 2, 3; 1, 2, 3, 4;, 1, 3, 4, 5]);
%! D = x2fx (x, 'interaction');
%! assert (D(1,:), [1, 1, 2, 3, 2, 3, 6]);
%! assert (D(2,:), [1, 2, 3, 4, 6, 8, 12]);
%! assert (D(3,:), [1, 3, 4, 5, 12, 15, 20]);
%! D = x2fx (x, 'quadratic');
%! assert (D(1,:), [1, 1, 2, 3, 2, 3, 6, 1, 4, 9]);
%! assert (D(2,:), [1, 2, 3, 4, 6, 8, 12, 4, 9, 16]);
%! assert (D(3,:), [1, 3, 4, 5, 12, 15, 20, 9, 16, 25]);
%! D = x2fx (x, 'purequadratic');
%! assert (D(1,:), [1, 1, 2, 3, 1, 4, 9]);
%! assert (D(2,:), [1, 2, 3, 4, 4, 9, 16]);
%! assert (D(3,:), [1, 3, 4, 5, 9, 16, 25]);

%!test
%! x = [1, 2, 3; 2, 3, 4; 3, 4, 5];
%! D = x2fx (x, [0, 0, 1; 1, 0, 2]);
%! assert (D, [3, 9; 4, 32; 5, 75]);

%!test
%! x = [1, 2, 3; 2, 3, 4; 3, 4, 5];
%! D = x2fx (x, 'linear', [1, 3]);
%! assert (D, [1, 1, 0, 2, 1, 0; 1, 0, 1, 3, 0, 1; 1, 0, 0, 4, 0, 0]);

%!test
%! x = [1, 2, 3; 2, 3, 4; 3, 4, 5];
%! D = x2fx (x, 'quadratic', [1, 3]);
%! assert (D(1,:), [1, 1, 0, 2, 1, 0, 2, 0, 1, 0, 0, 0, 2, 0, 4]);
%! assert (D(2,:), [1, 0, 1, 3, 0, 1, 0, 3, 0, 0, 0, 1, 0, 3, 9]);
%! assert (D(3,:), [1, 0, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 16]);

%!test
%! x = [1, 2, 3; 2, 3, 4; 3, 4, 5];
%! D = x2fx (x, 'cos');
%! assert (D(1,:), [0.5403, -0.4161, -0.9900], 1e-4);
%! assert (D(2,:), [-0.4161, -0.9900, -0.6536], 1e-4);
%! assert (D(3,:), [-0.9900, -0.6536, 0.2837], 1e-4);

## Tests for auto-detection of categorical data
%!test
%! x = [1, 2, 1; 2, 3, 2; 3, 4, 1];
%! D = x2fx (x, 'linear');
%! assert (size (D, 1), 3);
%! assert (isnumeric (D));

%!test
%! x = [1, 2, 1; 2, 3, 2; 3, 4, 1];
%! D = x2fx (x, 'linear');
%! assert (size (D), [3, 4]);

%!test
%! x = [1, 2, 1; 2, 3, 2; 3, 4, 1];
%! D = x2fx (x, 'linear', [3]);
%! assert (size (D, 1), 3);
%! assert (size (D, 2), 4);

%!test
%! x = [1, 2, 1; 2, 3, 2; 3, 4, 1];
%! D = x2fx (x, 'interaction', [3]);
%! assert (size (D, 1), 3);
%! assert (size (D, 2) > 5);

%!test
%! x = [1, 2, 1; 2, 3, 2; 3, 4, 1];
%! D1 = x2fx (x, 'linear', [3]);
%! assert (size (D1, 1), 3);
%! assert (size (D1, 2), 4);

%!test
%! x = [1, 2, 1; 2, 3, 2; 3, 4, 1];
%! D = x2fx (x, 'quadratic', [3]);
%! assert (size (D, 1), 3);
%! assert (size (D, 2) > 5);

%!test
%! x = [1, 2, 1; 2, 3, 2; 3, 4, 1];
%! D = x2fx (x, 'purequadratic', [3]);
%! assert (size (D, 1), 3);

%!test
%! x = [1, 2, 3; 2, 3, 4; 3, 4, 5];
%! D = x2fx (x, 'linear');
%! assert (size (D), [3, 4]);

%!test
%! x = [1, 2, 3; 2, 3, 4; 3, 4, 5];
%! D = x2fx (x, @cos);
%! assert (size (D), [3, 3]);

%!test
%! x = [1, 2; 2, 3; 3, 4];
%! D = x2fx (x, [0 0; 1 0; 0 1]);
%! assert (size (D, 1), 3);

%!test
%! x = [1, 2, 3; 2, 3, 4; 3, 4, 5];
%! [D, model, termstart, termend] = x2fx (x, 'linear');
%! assert (!isempty (model));
%! assert (!isempty (termstart));
%! assert (!isempty (termend));
%! assert (length (termstart) == length (termend));

%!error <x2fx: category index exceeds number of columns in X.> ...
%! x2fx ([1, 2, 3; 2, 3, 4], 'quadratic', [1, 4])
%!error <x2fx: unrecognized function 'cosine'.> ...
%! D = x2fx ([1, 2, 3; 2, 3, 4; 3, 4, 5], 'cosine')
%!error <x2fx: wrong number of columns between X and MODEL.> ...
%! x2fx ([1, 10; 2, 20; 3, 10], [0; 1]);
%!error <x2fx: wrong number of columns between X and MODEL.> ...
%! x2fx ([1, 10, 15; 2, 20, 40; 3, 10, 25], [0, 0; 1, 0; 0, 1; 1, 1; 2, 0]);
