## Copyright (C) 2026 Jayant Chauhan <0001jayant@gmail.com>
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
## @deftypefn {statistics} {[@var{intMat}, @var{response}, @var{hasIntercept}] =} FormulaParser (@var{formula}, @var{predictorNames})
## @deftypefnx {statistics} {@var{intMat} =} FormulaParser (@var{interactions}, @var{numPredictors})
##
## Parse a Wilkinson notation formula string or interaction specification into a binary term matrix.
##
## This utility is the central parsing engine for the statistics package. It converts
## high-level model specifications (strings like @qcode{"y ~ A * B"} or interaction
## flags) into the numeric design matrix logic required by regression and
## classification algorithms (e.g., @code{fitlm}, @code{fitrgam}).
##
## @strong{Modes of Operation:}
##
## @strong{1. Formula Parsing Mode}
## @*
## Usage: @code{FormulaParser ("y ~ A + B", @{"A", "B"@})}
## @*
## Parses a string containing a tilde (@qcode{~}) to separate the response variable
## from predictor terms. Supports complex expansion rules including factorials,
## nesting, and term deletion.
##
## @strong{2. Interaction Parsing Mode}
## @*
## Usage: @code{FormulaParser ("all", 3)}
## @*
## Generates an interaction matrix based on a numeric or string descriptor. Useful
## for automated model building where a specific formula is not provided.
##
## @strong{Input Arguments:}
##
## @multitable @columnfractions 0.2 0.8
## @item @var{formula}
## @tab A character string specifying the model relationship. The string must contain
## a tilde (@qcode{~}) separating the response (LHS) from the terms (RHS).
## @item @var{predictorNames}
## @tab A cell array of strings containing the names of valid predictor variables
## in the dataset. Used to map string terms to column indices.
## @item @var{interactions}
## @tab Alternative input. Can be a numeric scalar (max interaction degree),
## a logical matrix (manual specification), or the string @qcode{"all"}.
## @item @var{numPredictors}
## @tab Integer number of total predictor variables (used only in Interaction Mode).
## @end multitable
##
## @strong{Output Arguments:}
##
## @multitable @columnfractions 0.2 0.8
## @item @var{intMat}
## @tab A logical matrix of size @math{K \times P}, where @math{K} is the number
## of terms and @math{P} is the number of predictors. A value of @code{true} (1)
## at @code{(i,j)} indicates that predictor @math{j} is included in term @math{i}.
## @item @var{response}
## @tab The name of the response variable extracted from the Left-Hand Side (LHS)
## of the formula. Returns an empty string if LHS is missing (invalid syntax).
## @item @var{hasIntercept}
## @tab A boolean flag indicating whether an intercept term should be included.
## Default is @code{true} unless explicitly removed via @qcode{"- 1"}.
## @end multitable
##
## @strong{Supported Wilkinson Syntax:}
##
## @table @code
## @item +
## Additive term. @qcode{"A + B"} creates two separate terms: Main Effect A and
## Main Effect B.
## @item -
## Deletion term. @qcode{"A * B - B"} removes the main effect B from the expansion,
## leaving only A and the A:B interaction.
## @item :
## Interaction term. @qcode{"A:B"} creates a single term representing the
## interaction between A and B (no main effects included).
## @item *
## Factorial expansion. @qcode{"A * B"} is shorthand for @qcode{"A + B + A:B"}.
## It includes both main effects and their interaction.
## @item /
## Nesting. @qcode{"A / B"} is shorthand for @qcode{"A + A:B"}. It represents
## B nested within A.
## @item - 1
## Intercept removal. @qcode{"y ~ x - 1"} forces the model to pass through the
## origin by setting @var{hasIntercept} to @code{false}.
## @item ()
## Grouping. Parentheses can be used to group terms for operators. For example,
## @qcode{"(A + B):C"} distributes the interaction to expand to @qcode{"A:C + B:C"}.
## @end table
##
## @strong{Note:}
## The Power operator (@qcode{^}) and recursive deletion operators (@qcode{-*}, @qcode{-/})
## are currently not implemented.
##
## @strong{Term Ordering:}
## The resulting @var{intMat} is sorted to ensure consistency with standard model
## hierarchies:
## @enumerate
## @item Main effects appear first (lower complexity).
## @item Interaction terms appear subsequently (higher complexity).
## @item Within the same complexity level, terms are sorted by the order of
## variables in @var{predictorNames}.
## @end enumerate
##
## @seealso{fitcgam, fitrgam, fitlm}
## @end deftypefn

function [intMat, response, hasIntercept] = FormulaParser (varargin)

  intMat = [];
  response = "";
  hasIntercept = true;

  if (nargin != 2)
    error ("FormulaParser: invalid input arguments.");
  endif

  arg1 = varargin{1};
  arg2 = varargin{2};

  ## --- Mode 1: Formula Parsing ---
  ## Condition: arg2 is a cell array (predictor names)
  if (iscellstr (arg2))
    formula = arg1;
    predictorNames = arg2;

    if (! ischar (formula))
       error ("FormulaParser: Formula must be a string.");
    endif

    ## 1. Syntax Check
    if (isempty (strfind (formula, '~')))
      error ("FormulaParser: invalid syntax. Formula must contain '~'.");
    endif

    ## 2. Split LHS (Response) and RHS (Predictors)
    parts = strsplit (formula, '~');
    if (numel (parts) < 2)
      error ("FormulaParser: no predictor terms found.");
    endif

    response = strtrim (parts{1});
    rhs = strtrim (parts{2});

    if (isempty (rhs))
      error ("FormulaParser: no predictor terms found.");
    endif

    ## 3. Handle Intercepts (+1 / -1)
    if (! isempty (strfind (rhs, '-1')) || ! isempty (strfind (rhs, '- 1')))
      hasIntercept = false;
      rhs = strrep (rhs, '- 1', '');
      rhs = strrep (rhs, '-1', '');
    endif

    rhs = strrep (rhs, '+ 1', '');
    rhs = strrep (rhs, '+1', '');

    ## 4. Recursive Expansion
    termsList = parse_expression (rhs, predictorNames);

    ## Convert list to matrix (cell2mat handles empty list correctly)
    intMat = cell2mat (termsList');

    ## 5. Sort and Unique
    if (! isempty (intMat))
      intMat = unique (intMat, 'rows');
      termComplexity = sum (intMat, 2);
      [~, idx] = sortrows ([termComplexity, -intMat]);
      intMat = intMat(idx, :);
    else
      ## This error block is reachable if termsList was empty (e.g. "y ~ -1")
      error ("FormulaParser: formula resulted in no valid terms.");
    endif

  ## --- Mode 2: Interaction Parsing ---
  ## Condition: arg2 is a numeric scalar (number of predictors)
  elseif (isnumeric (arg2) && isscalar (arg2))
    interactions = arg1;
    numPredictors = arg2;

    if (islogical (interactions))
      if (numPredictors != columns (interactions))
        error ("FormulaParser: columns in Interactions logical matrix must equal the number of predictors.");
      endif
      intMat = interactions;

    elseif (isnumeric (interactions) && isscalar (interactions))
      if (interactions < 0)
        error ("FormulaParser: Invalid interaction argument.");
      elseif (interactions > numPredictors * (numPredictors - 1) / 2)
        error ("FormulaParser: number of interaction terms requested is larger than all possible combinations.");
      endif

      allMat = ff2n (numPredictors);
      iterms = find (sum (allMat, 2) > 1);

      ## Cast to logical to match Mode 1 output
      intMat = logical (allMat(iterms, :));

    elseif (ischar (interactions) && strcmpi (interactions, "all"))
      allMat = ff2n (numPredictors);
      iterms = find (sum (allMat, 2) > 1);

      ## Cast to logical to match Mode 1 output
      intMat = logical (allMat(iterms, :));
    else
      error ("FormulaParser: Invalid interaction argument.");
    endif

    ## Sort Mode 2 results
    if (! isempty (intMat))
       [~, idx] = sortrows ([sum(intMat, 2), -intMat]);
       intMat = intMat(idx, :);
    endif

  else
    error ("FormulaParser: invalid input arguments. Second argument must be PredictorNames (cellstr) or NumPredictors (scalar).");
  endif

endfunction

## --- Helper Functions for Recursive Parsing ---

function terms = parse_expression (str, pNames)
  terms = {};
  current_op = '+'; # Implicit start with addition

  ## We must scan for '+' and '-' at the top level (depth 0)
  depth = 0;
  start_idx = 1;
  len = length (str);

  for i = 1:len
    c = str(i);
    if (c == '(')
      depth += 1;
    elseif (c == ')')
      depth -= 1;
    elseif (depth == 0 && (c == '+' || c == '-'))
      ## Found an operator. Process the chunk before it.
      chunk = strtrim (str(start_idx:i-1));
      if (! isempty (chunk))
        subTerms = parse_term (chunk, pNames);
        if (current_op == '+')
          terms = [terms, subTerms];
        else
          terms = remove_terms (terms, subTerms);
        endif
      endif
      ## Update operator and reset start index
      current_op = c;
      start_idx = i + 1;
    endif
  endfor

  ## Process the final chunk
  chunk = strtrim (str(start_idx:end));
  if (! isempty (chunk))
    subTerms = parse_term (chunk, pNames);
    if (current_op == '+')
      terms = [terms, subTerms];
    else
      terms = remove_terms (terms, subTerms);
    endif
  endif
endfunction

function terms = remove_terms (current_terms, terms_to_remove)
  ## Set difference for cell array of term vectors
  if (isempty (current_terms) || isempty (terms_to_remove))
    terms = current_terms;
    return;
  endif

  ## Convert to matrices for efficient row comparison
  mat_curr = cell2mat (current_terms');
  mat_rem  = cell2mat (terms_to_remove');

  ## Use setdiff with 'rows' to find remaining terms
  remaining_rows = setdiff (mat_curr, mat_rem, 'rows');

  ## Convert back to cell array
  if (isempty (remaining_rows))
    terms = {};
  else
    terms = num2cell (remaining_rows, 2)';
  endif
endfunction

function terms = parse_term (str, pNames)

  ## Case 1: Nesting (/)
  ## Rule: A / B -> A + A:B
  ## Rule: A / B / C -> A + A:B + A:B:C
  factors = paren_split (str, '/');
  if (numel (factors) > 1)
    parsed_factors = {};
    for i = 1:numel (factors)
      parsed_factors{end+1} = parse_term (factors{i}, pNames);
    endfor

    terms = {};
    ## Start with the first factor (the hierarchy root)
    cumul_inter = parsed_factors{1};
    terms = [terms, cumul_inter];

    ## For subsequent factors, interact them with the cumulative chain
    for k = 2:numel (parsed_factors)
      next_fac = parsed_factors{k};
      new_cumul = {};
      ## Distribute interaction: cumul_inter : next_fac
      for t1 = 1:numel (cumul_inter)
        for t2 = 1:numel (next_fac)
          new_cumul{end+1} = cumul_inter{t1} | next_fac{t2};
        endfor
      endfor
      cumul_inter = new_cumul;
      terms = [terms, cumul_inter];
    endfor
    return;
  endif

  ## Case 2: Factorial Expansion (*)
  factors = paren_split (str, '*');
  if (numel (factors) > 1)
    expanded_factors = {};
    for i = 1:numel (factors)
      expanded_factors{end+1} = parse_term (factors{i}, pNames);
    endfor

    n_fac = numel (factors);
    combo_map = ff2n (n_fac);
    combo_map(1, :) = [];

    terms = {};
    for r = 1:rows (combo_map)
      active_indices = find (combo_map(r, :));
      current_terms = expanded_factors{active_indices(1)};

      for k = 2:numel (active_indices)
        next_terms = expanded_factors{active_indices(k)};
        new_combo = {};
        for t1 = 1:numel (current_terms)
          for t2 = 1:numel (next_terms)
            new_combo{end+1} = current_terms{t1} | next_terms{t2};
          endfor
        endfor
        current_terms = new_combo;
      endfor
      terms = [terms, current_terms];
    endfor
    return;
  endif

  ## Case 3: Interaction (:)
  factors = paren_split (str, ':');
  if (numel (factors) > 1)
    component_terms = {};
    for i = 1:numel (factors)
      component_terms{end+1} = parse_term (factors{i}, pNames);
    endfor

    terms = component_terms{1};
    for k = 2:numel (component_terms)
      next_set = component_terms{k};
      product_set = {};
      for t1 = 1:numel (terms)
        for t2 = 1:numel (next_set)
          product_set{end+1} = terms{t1} | next_set{t2};
        endfor
      endfor
      terms = product_set;
    endfor
    return;
  endif

  ## Case 4: Parentheses
  if (length (str) >= 2 && str(1) == '(' && str(end) == ')')
    terms = parse_expression (str(2:end-1), pNames);
    return;
  endif

  ## Case 5: Base Predictor
  idx = find (strcmp (pNames, str));
  if (isempty (idx))
    if (isempty (str))
       error ("FormulaParser: invalid syntax (empty term derived).");
    endif
    error ("FormulaParser: predictor '%s' not found in PredictorNames.", str);
  endif

  termVec = logical (zeros (1, numel (pNames)));
  termVec(idx) = 1;
  terms = {termVec};

endfunction

function parts = paren_split (str, delim)
  parts = {};
  depth = 0;
  start_idx = 1;
  len = length (str);

  for i = 1:len
    c = str(i);
    if (c == '(')
      depth += 1;
    elseif (c == ')')
      depth -= 1;
    elseif (c == delim && depth == 0)
      parts{end+1} = strtrim (str(start_idx:i-1));
      start_idx = i + 1;
    endif
  endfor

  parts{end+1} = strtrim (str(start_idx:end));
endfunction

%!demo
%! ## 1. Basic Additive Model
%! ## Defines a model with two main effects (Age, Weight) and an intercept.
%! ## The response 'BP' is extracted automatically.
%! pnames = {"Age", "Height", "Weight"};
%! formula = "BP ~ Age + Weight";
%! [intMat, resp, hasInt] = FormulaParser (formula, pnames)

%!demo
%! ## 2. Interaction Term
%! ## Defines a model with main effects A, B and their specific interaction A:B.
%! ## Note: The term 'A:B' in the matrix will have 1s at both A and B's indices.
%! pnames = {"A", "B", "C"};
%! formula = "y ~ A + B + A:B";
%! [intMat, resp, hasInt] = FormulaParser (formula, pnames)

%!demo
%! ## 3. Factorial Expansion (*)
%! ## The * operator is shorthand for Main Effects + Interaction.
%! ## "A * B" automatically expands to "A + B + A:B".
%! pnames = {"A", "B"};
%! formula = "y ~ A * B";
%! [intMat, resp, ~] = FormulaParser (formula, pnames)

%!demo
%! ## 4. Intercept Removal
%! ## Use "- 1" to exclude the constant term from the model.
%! ## 'hasInt' will return false (0).
%! pnames = {"x1", "x2"};
%! formula = "y ~ x1 + x2 - 1";
%! [intMat, resp, hasInt] = FormulaParser (formula, pnames)

%!demo
%! ## 5. Nested Grouping and Distribution
%! ## Parentheses can be used to group terms. "(A + B):C" distributes the
%! ## interaction to expand to "A:C + B:C".
%! pnames = {"A", "B", "C"};
%! formula = "Y ~ (A + B):C";
%! [intMat, ~, ~] = FormulaParser (formula, pnames)

%!demo
%! ## 6. Generating 'All' Interactions
%! ## Generates all possible pairwise and higher-order interactions.
%! ## Useful for generating a hypothesis space without a formula.
%! numPreds = 3;
%! intMat = FormulaParser ("all", numPreds)

%!test
%! ## Test 1: Simple additive formula (Standard Regression)
%! pnames = {"x1", "x2", "x3"};
%! [intMat, resp, hasInt] = FormulaParser ("y ~ x1 + x2", pnames);
%! expected = logical ([1, 0, 0; 0, 1, 0]);
%! assert (intMat, expected);
%! assert (resp, "y");
%! assert (hasInt, true);

%!test
%! ## Test 2: Factorial Expansion and Ordering (A * B)
%! ## Checks if A*B expands correctly to Main Effects + Interaction
%! pnames = {"A", "B"};
%! [intMat, resp, ~] = FormulaParser ("y ~ A * B", pnames);
%! expected = logical ([1, 0; 0, 1; 1, 1]);
%! assert (intMat, expected);

%!test
%! ## Test 3: Intercept Removal (-1)
%! pnames = {"A", "B"};
%! [~, ~, hasInt] = FormulaParser ("Cost ~ A + B - 1", pnames);
%! assert (hasInt, false);

%!test
%! ## Test 4: Grouping Distribution ((A + B):C)
%! ## Verifies distributive property: A:C + B:C
%! pnames = {"A", "B", "C"};
%! [intMat, ~, ~] = FormulaParser ("y ~ (A + B):C", pnames);
%! expected = logical ([1, 0, 1; 0, 1, 1]);
%! assert (intMat, expected);

%!test
%! ## Test 5: Complex Nested Factorial ((A + B) * C)
%! ## Expands to Main(A,B,C) + Inter(A:C, B:C)
%! pnames = {"A", "B", "C"};
%! [intMat, ~, ~] = FormulaParser ("y ~ (A + B) * C", pnames);
%! expected = logical ([1,0,0; 0,1,0; 0,0,1; 1,0,1; 0,1,1]);
%! assert (intMat, expected);

%!test
%! ## Test 6: Whitespace Robustness
%! pnames = {"x1", "x2"};
%! [intMat, ~, ~] = FormulaParser ("y ~ x1 + x1 : x2", pnames);
%! assert (intMat, logical ([1, 0; 1, 1]));

%!test
%! ## Test 7: Interaction Argument 'all' (Package Feature)
%! ## 3 predictors -> A:B, A:C, B:C, A:B:C
%! intMat = FormulaParser ("all", 3);
%! expected = logical ([1, 1, 0; 1, 0, 1; 0, 1, 1; 1, 1, 1]);
%! assert (intMat, expected);

%!test
%! ## Test 8: Error Handling (Trailing Syntax)
%! pnames = {"x"};
%! fail ("FormulaParser ('y ~ x:', pnames)", "invalid syntax");

%!test
%! ## Test 9: Triple Factorial Expansion (A * B * C)
%! ## Checks full model expansion: Mains + Pairs + Triple
%! pnames = {"A", "B", "C"};
%! [intMat, ~, ~] = FormulaParser ("y ~ A * B * C", pnames);
%! expected = logical ([1,0,0; 0,1,0; 0,0,1; 1,1,0; 1,0,1; 0,1,1; 1,1,1]);
%! assert (intMat, expected);

%!test
%! ## Test 10: Mixed Operators and Precedence (A + B:C + D)
%! ## Ensures mixing + and : works without grouping
%! pnames = {"A", "B", "C", "D"};
%! [intMat, ~, ~] = FormulaParser ("y ~ A + B : C + D", pnames);
%! expected = logical ([1,0,0,0; 0,0,0,1; 0,1,1,0]);
%! assert (intMat, expected);

%!test
%! ## Test 11: Real-World Model Specification
%! ## y ~ A + B + C + A:B + B:C (Specific manual interactions)
%! pnames = {"A", "B", "C"};
%! [intMat, ~, ~] = FormulaParser ("y ~ A + B + C + A:B + B:C", pnames);
%! expected = logical ([1,0,0; 0,1,0; 0,0,1; 1,1,0; 0,1,1]);
%! assert (intMat, expected);

%!test
%! ## Test 12: High-Order Interaction Chain (A:B:C)
%! ## Verifies pure 3-way interaction term
%! pnames = {"A", "B", "C"};
%! [intMat, ~, ~] = FormulaParser ("y ~ A : B : C", pnames);
%! expected = logical ([1, 1, 1]);
%! assert (intMat, expected);

%!test
%! ## Test 13: Substring Variable Names (Partial Match Check)
%! ## Ensures 'Age' does not incorrectly match 'AgeGroup'
%! pnames = {"Age", "AgeGroup"};
%! [intMat, ~, ~] = FormulaParser ("y ~ Age + AgeGroup", pnames);
%! expected = logical ([1, 0; 0, 1]);
%! assert (intMat, expected);

%!test
%! ## Test 14: Associativity of Factorial
%! ## A * (B * C) should equal A * B * C
%! pnames = {"A", "B", "C"};
%! mat1 = FormulaParser ("y ~ A * (B * C)", pnames);
%! mat2 = FormulaParser ("y ~ A * B * C", pnames);
%! assert (isequal (mat1, mat2));

%!test
%! ## Test 15: Variable Names with Special Characters
%! ## Robustness check for underscores/numbers in variable names
%! pnames = {"Var_1", "Var_2"};
%! [intMat, ~, ~] = FormulaParser ("y ~ Var_1 * Var_2", pnames);
%! expected = logical ([1, 0; 0, 1; 1, 1]);
%! assert (intMat, expected);

%!test
%! ## Test 16: Deletion Operator (-)
%! ## Wilkinson (1973): A*B - B should yield A + A:B
%! pnames = {"A", "B"};
%! [intMat, ~, ~] = FormulaParser ("y ~ A * B - B", pnames);
%! expected = logical ([1, 0; 1, 1]); %% A(1,0) and A:B(1,1). B(0,1) is removed.
%! assert (intMat, expected);

%!test
%! ## Test 17: Nesting Operator (/)
%! ## Wilkinson (1973): A/B expands to A + A:B
%! pnames = {"A", "B"};
%! [intMat, ~, ~] = FormulaParser ("y ~ A / B", pnames);
%! expected = logical ([1, 0; 1, 1]);
%! assert (intMat, expected);

%!test
%! ## Test 18: Nested Chain (A/B/C)
%! ## Expands to A + A:B + A:B:C
%! pnames = {"A", "B", "C"};
%! [intMat, ~, ~] = FormulaParser ("y ~ A / B / C", pnames);
%! expected = logical ([1,0,0; 1,1,0; 1,1,1]);
%! assert (intMat, expected);

%!error <FormulaParser: invalid syntax. Formula must contain '~'.> FormulaParser ("y x1 + x2", {"x1"})
%!error <FormulaParser: no predictor terms found.> FormulaParser ("y ~ ", {"x1"})
%!error <FormulaParser: predictor 'z' not found in PredictorNames.> FormulaParser ("y ~ x1 + z", {"x1"})
%!error <FormulaParser: invalid syntax> FormulaParser ("y ~ x1:", {"x1"})
%!error <FormulaParser: formula resulted in no valid terms.> FormulaParser ("y ~ -1", {"x1"})
