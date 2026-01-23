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
## from predictor terms. Supports complex expansion rules including factorials and
## nesting.
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
## @item :
## Interaction term. @qcode{"A:B"} creates a single term representing the
## interaction between A and B (no main effects included).
## @item *
## Factorial expansion. @qcode{"A * B"} is shorthand for @qcode{"A + B + A:B"}.
## It includes both main effects and their interaction.
## @item - 1
## Intercept removal. @qcode{"y ~ x - 1"} forces the model to pass through the
## origin by setting @var{hasIntercept} to @code{false}.
## @item ()
## Grouping. Parentheses can be used to group terms for operators. For example,
## @qcode{"(A + B):C"} distributes the interaction to expand to @qcode{"A:C + B:C"}.
## @end table
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
  blocks = paren_split (str, '+');

  for i = 1:numel (blocks)
    block = strtrim (blocks{i});
    if (isempty (block))
      continue;
    endif
    subTerms = parse_term (block, pNames);
    terms = [terms, subTerms];
  endfor
endfunction

function terms = parse_term (str, pNames)

  ## Case 1: Factorial Expansion (*)
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

  ## Case 2: Interaction (:)
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

  ## Case 3: Parentheses
  if (length (str) >= 2 && str(1) == '(' && str(end) == ')')
    terms = parse_expression (str(2:end-1), pNames);
    return;
  endif

  ## Case 4: Base Predictor
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
%! ## Test 1: Simple additive formula
%! pnames = {"x1", "x2", "x3"};
%! [intMat, resp, hasInt] = FormulaParser ("y ~ x1 + x2", pnames);
%! expected = logical ([1, 0, 0; 0, 1, 0]);
%! assert (intMat, expected);
%! assert (resp, "y");
%! assert (hasInt, true);

%!test
%! ## Test 2: Factorial Expansion and Ordering
%! pnames = {"A", "B"};
%! ## A * B -> A + B + A:B
%! ## Sorted: A (1,0), B (0,1), AB (1,1)
%! [intMat, resp, ~] = FormulaParser ("y ~ A * B", pnames);
%! expected = logical ([1, 0; 0, 1; 1, 1]);
%! assert (intMat, expected);

%!test
%! ## Test 3: Intercept Removal
%! pnames = {"A", "B"};
%! [~, ~, hasInt] = FormulaParser ("Cost ~ A + B - 1", pnames);
%! assert (hasInt, false);

%!test
%! ## Test 4: Nested Grouping Distribution
%! pnames = {"A", "B", "C"};
%! ## (A + B):C -> A:C + B:C
%! ## A:C is [1 0 1], B:C is [0 1 1]
%! [intMat, ~, ~] = FormulaParser ("y ~ (A + B):C", pnames);
%! expected = logical ([1, 0, 1; 0, 1, 1]);
%! assert (intMat, expected);

%!test
%! ## Test 5: Complex Nested Factorial
%! pnames = {"A", "B", "C"};
%! ## (A+B)*C -> A, B, C, A:C, B:C
%! [intMat, ~, ~] = FormulaParser ("y ~ (A + B) * C", pnames);
%! expected = logical ([1,0,0; 0,1,0; 0,0,1; 1,0,1; 0,1,1]);
%! assert (intMat, expected);

%!test
%! ## Test 6: Whitespace handling
%! pnames = {"x1", "x2"};
%! [intMat, ~, ~] = FormulaParser ("y ~ x1 + x1 : x2", pnames);
%! assert (intMat, logical ([1, 0; 1, 1]));

%!test
%! ## Test 7: Interaction Argument 'all' (Sorted)
%! ## 3 predictors -> A:B, A:C, B:C, A:B:C
%! intMat = FormulaParser ("all", 3);
%! expected = logical ([1, 1, 0; 1, 0, 1; 0, 1, 1; 1, 1, 1]);
%! assert (intMat, expected);

%!test
%! ## Test 8: Trailing Syntax Error
%! pnames = {"x"};
%! fail ("FormulaParser ('y ~ x:', pnames)", "invalid syntax");

%!test
%! ## Test 9: Deep Nesting Distribution
%! ## ((A + B) : C) : D -> (A:C + B:C) : D -> A:C:D + B:C:D
%! pnames = {"A", "B", "C", "D"};
%! [intMat, ~, ~] = FormulaParser ("y ~ ((A + B) : C) : D", pnames);
%! ## A:C:D [1 0 1 1], B:C:D [0 1 1 1]
%! expected = logical ([1, 0, 1, 1; 0, 1, 1, 1]);
%! assert (intMat, expected);

%!test
%! ## Test 10: FOIL Distribution (Two Groups)
%! ## (A + B) : (C + D) -> A:C + A:D + B:C + B:D
%! pnames = {"A", "B", "C", "D"};
%! [intMat, ~, ~] = FormulaParser ("y ~ (A + B) : (C + D)", pnames);
%! expected = logical ([1,0,1,0; 1,0,0,1; 0,1,1,0; 0,1,0,1]);
%! assert (intMat, expected);

%!test
%! ## Test 11: Mixed Factorial and Interaction
%! ## (A * B) : C -> (A + B + A:B) : C -> A:C + B:C + A:B:C
%! pnames = {"A", "B", "C"};
%! [intMat, ~, ~] = FormulaParser ("y ~ (A * B) : C", pnames);
%! expected = logical ([1, 0, 1; 0, 1, 1; 1, 1, 1]);
%! assert (intMat, expected);

%!test
%! ## Test 12: Triple Factorial Expansion
%! ## A * B * C -> A + B + C + AB + AC + BC + ABC
%! pnames = {"A", "B", "C"};
%! [intMat, ~, ~] = FormulaParser ("y ~ A * B * C", pnames);
%! ## Order: Main(A,B,C), Pairs(AB,AC,BC), Trip(ABC)
%! expected = logical ([1,0,0; 0,1,0; 0,0,1; 1,1,0; 1,0,1; 0,1,1; 1,1,1]);
%! assert (intMat, expected);

%!test
%! ## Test 13: Redundant Terms Handling
%! ## A + A + (A:B) + (B:A) -> Collapses to unique terms: A, A:B
%! pnames = {"A", "B"};
%! [intMat, ~, ~] = FormulaParser ("y ~ A + A + (A:B) + (B:A)", pnames);
%! expected = logical ([1, 0; 1, 1]);
%! assert (intMat, expected);

%!test
%! ## Test 14: Parentheses Stripping
%! ## (A) + ((B)) -> A + B
%! pnames = {"A", "B"};
%! [intMat, ~, ~] = FormulaParser ("y ~ (A) + ((B))", pnames);
%! expected = logical ([1, 0; 0, 1]);
%! assert (intMat, expected);

%!test
%! ## Test 15: Interaction of Factorials (Stress Test)
%! ## (A * B) : (C * D) -> (A+B+AB):(C+D+CD) -> 9 terms (all combos of left x right)
%! pnames = {"A", "B", "C", "D"};
%! [intMat, ~, ~] = FormulaParser ("y ~ (A * B) : (C * D)", pnames);
%! ## Should produce 9 rows. We just check size to be safe.
%! assert (size (intMat), [9, 4]);
%! ## Check one complex term: A:B:C:D (1 1 1 1) should be present
%! assert (ismember ([1, 1, 1, 1], intMat, "rows"));

%!test
%! ## Test 16: Intercept inside grouping (Should treat '-1' as token removal)
%! ## y ~ (A + B) - 1 -> A + B, no intercept
%! pnames = {"A", "B"};
%! [intMat, ~, hasInt] = FormulaParser ("y ~ (A + B) - 1", pnames);
%! assert (hasInt, false);
%! assert (intMat, logical ([1, 0; 0, 1]));

%!test
%! ## Test 17: Complex Real-World (Main + Specific Interactions)
%! ## y ~ A + B + C + A:B + B:C
%! pnames = {"A", "B", "C"};
%! [intMat, ~, ~] = FormulaParser ("y ~ A + B + C + A:B + B:C", pnames);
%! expected = logical ([1,0,0; 0,1,0; 0,0,1; 1,1,0; 0,1,1]);
%! assert (intMat, expected);

%!test
%! ## Test 18: Empty Grouping / Edge Case
%! ## y ~ A + () -> A (Parser should skip empty group)
%! pnames = {"A"};
%! [intMat, ~, ~] = FormulaParser ("y ~ A + ()", pnames);
%! assert (intMat, logical ([1]));

%!test
%! ## Test 19: High-Arity Distribution (FOIL on Steroids)
%! ## (A + B + C) : (D + E) -> 3x2 = 6 interaction terms
%! pnames = {"A", "B", "C", "D", "E"};
%! [intMat, ~, ~] = FormulaParser ("y ~ (A + B + C) : (D + E)", pnames);
%! assert (rows (intMat), 6);
%! ## Check A:D term [1 0 0 1 0]
%! assert (ismember ([1, 0, 0, 1, 0], intMat, "rows"));

%!test
%! ## Test 20: Operator Precedence (Standard R/Wilkinson)
%! ## A : B * C -> A : (B * C) [because * is expansion, : is interaction]
%! ## Expands to: A : (B + C + B:C) -> A:B + A:C + A:B:C
%! pnames = {"A", "B", "C"};
%! [intMat, ~, ~] = FormulaParser ("y ~ A : B * C", pnames);
%! ## A:B (110), C (001), A:B:C (111)
%! ## Sorted: C (001), A:B (110), A:B:C (111)
%! expected = logical ([0, 0, 1; 1, 1, 0; 1, 1, 1]);
%! assert (intMat, expected);

%!test
%! ## Test 21: Extreme Parenthesis Depth
%! ## ((((A)))) + B -> A + B
%! pnames = {"A", "B"};
%! [intMat, ~, ~] = FormulaParser ("y ~ ((((A)))) + B", pnames);
%! assert (intMat, logical ([1, 0; 0, 1]));

%!test
%! ## Test 22: Self-Interaction Redundancy
%! ## A : A should behave as A (Logical OR 1|1 = 1)
%! pnames = {"A"};
%! [intMat, ~, ~] = FormulaParser ("y ~ A : A", pnames);
%! assert (intMat, logical ([1]));

%!test
%! ## Test 23: Multi-Level Factorial
%! ## (A*B)*(C*D) -> (A+B+AB) * (C+D+CD) -> Full cross product
%! ## Result should have 3*3 = 9 terms, then their interaction expansions...
%! ## Actually A*B*C*D generates 15 terms (2^4 - 1).
%! ## (A*B)*(C*D) is structurally A*B*C*D.
%! pnames = {"A", "B", "C", "D"};
%! [intMat, ~, ~] = FormulaParser ("y ~ (A * B) * (C * D)", pnames);
%! assert (rows (intMat), 15);
%! assert (ismember ([1, 1, 1, 1], intMat, "rows"));

%!test
%! ## Test 24: Consistency Check ("all" vs formula)
%! ## "all" with 2 vars should equal "A * B"
%! pnames = {"A", "B"};
%! mat1 = FormulaParser ("y ~ A * B", pnames);
%! mat2 = FormulaParser ("all", 2);
%! ## Note: "all" generates interactions, A*B generates main+inter
%! ## Wait, mode 2 "all" in FormulaParser returns INTERACTIONS ONLY (sum > 1).
%! ## "A*B" returns Main + Inter.
%! ## So we expect mat1 to contain mat2.
%! assert (all (ismember (mat2, mat1, "rows")));

%!test
%! ## Test 25: Variable Names with Numbers and Underscores
%! pnames = {"Var_1", "Var_2"};
%! [intMat, ~, ~] = FormulaParser ("y ~ Var_1 * Var_2", pnames);
%! expected = logical ([1, 0; 0, 1; 1, 1]);
%! assert (intMat, expected);

%!test
%! ## Test 26: Interleaved Operators
%! ## A + B : C + D -> A, B:C, D
%! ## Sorted: A (1000), D (0001), B:C (0110)
%! pnames = {"A", "B", "C", "D"};
%! [intMat, ~, ~] = FormulaParser ("y ~ A + B : C + D", pnames);
%! expected = logical ([1,0,0,0; 0,0,0,1; 0,1,1,0]);
%! assert (intMat, expected);

%!test
%! ## Test 27: The "Empty Set" Interaction / Redundancy
%! ## A : (B + B) -> A : B
%! pnames = {"A", "B"};
%! [intMat, ~, ~] = FormulaParser ("y ~ A : (B + B)", pnames);
%! assert (intMat, logical ([1, 1]));

%!test
%! ## Test 28: Complex Whitespace Abuse
%! ## y ~  A   * B
%! pnames = {"A", "B"};
%! [intMat, ~, ~] = FormulaParser ("y ~   A   * B   ", pnames);
%! assert (intMat, logical ([1, 0; 0, 1; 1, 1]));

%!test
%! ## Test 29: Triple Interaction Chain
%! ## A : B : C
%! pnames = {"A", "B", "C"};
%! [intMat, ~, ~] = FormulaParser ("y ~ A : B : C", pnames);
%! expected = logical ([1, 1, 1]);
%! assert (intMat, expected);

%!test
%! ## Test 30: Factorial inside Interaction inside Factorial (Torture Test)
%! ## ((A*B):C)*D -> ( (A+B+AB):C ) * D -> (AC + BC + ABC) * D
%! ## -> AC + BC + ABC + D + ACD + BCD + ABCD
%! pnames = {"A", "B", "C", "D"};
%! [intMat, ~, ~] = FormulaParser ("y ~ ((A*B):C)*D", pnames);
%! ## AC(1010), BC(0110), ABC(1110)
%! ## D(0001)
%! ## ACD(1011), BCD(0111), ABCD(1111)
%! expected = logical ([1,0,1,0; 0,1,1,0; 1,1,1,0; 0,0,0,1; 1,0,1,1; 0,1,1,1; 1,1,1,1]);
%! ## Ordering might differ, use unique/sort check implied by function
%! assert (rows (intMat), 7);
%! assert (all (ismember (expected, intMat, "rows")));

%!test
%! ## Test 31: Substring Variable Names
%! ## Predictor names where one is a substring of another.
%! pnames = {"Age", "AgeGroup"};
%! ## If parser is not exact, "Age" might match "AgeGroup".
%! [intMat, ~, ~] = FormulaParser ("y ~ Age + AgeGroup", pnames);
%! expected = logical ([1, 0; 0, 1]);
%! assert (intMat, expected);

%!test
%! ## Test 32: Associativity of Factorial
%! ## A * (B * C) should equal A * B * C
%! pnames = {"A", "B", "C"};
%! mat1 = FormulaParser ("y ~ A * (B * C)", pnames);
%! mat2 = FormulaParser ("y ~ A * B * C", pnames);
%! assert (isequal (mat1, mat2));

%!test
%! ## Test 33: Three-Way Group Distribution (3-Way FOIL)
%! ## (A+B) : (C+D) : (E+F) -> 2*2*2 = 8 terms
%! pnames = {"A", "B", "C", "D", "E", "F"};
%! [intMat, ~, ~] = FormulaParser ("y ~ (A+B):(C+D):(E+F)", pnames);
%! assert (rows (intMat), 8);
%! ## Check A:C:E (1 0 1 0 1 0)
%! assert (ismember ([1, 0, 1, 0, 1, 0], intMat, "rows"));

%!test
%! ## Test 34: Self-Factorial
%! ## A * A -> A + A + A:A -> A
%! pnames = {"A"};
%! [intMat, ~, ~] = FormulaParser ("y ~ A * A", pnames);
%! assert (intMat, logical ([1]));

%!test
%! ## Test 35: Complex Mixed Structure
%! ## (A + B:C) * D -> (A + BC) * D -> A + BC + D + AD + BCD
%! pnames = {"A", "B", "C", "D"};
%! [intMat, ~, ~] = FormulaParser ("y ~ (A + B:C) * D", pnames);
%! ## A(1000), BC(0110), D(0001), AD(1001), BCD(0111)
%! expected = logical ([1,0,0,0; 0,1,1,0; 0,0,0,1; 1,0,0,1; 0,1,1,1]);
%! assert (rows (intMat), 5);
%! assert (all (ismember (expected, intMat, "rows")));

%!test
%! ## Test 36: Duplicate Group Factorial
%! ## (A+B) * (A+B) -> Equivalent to (A+B) * (A+B)
%! ## (A+B) * (A+B) -> A+B + A+B + (A+B):(A+B)
%! ## -> A + B + (AA + AB + BA + BB) -> A + B + A + AB + AB + B
%! ## -> A + B + AB.
%! pnames = {"A", "B"};
%! [intMat, ~, ~] = FormulaParser ("y ~ (A+B) * (A+B)", pnames);
%! expected = logical ([1, 0; 0, 1; 1, 1]);
%! assert (intMat, expected);

%!test
%! ## Test 37: Leading Operator Resilience
%! ## + A + B should parse as A + B
%! pnames = {"A", "B"};
%! [intMat, ~, ~] = FormulaParser ("y ~ + A + B", pnames);
%! assert (intMat, logical ([1, 0; 0, 1]));

%!test
%! ## Test 38: Deeply Nested Interaction Chain
%! ## A : (B : (C : D)) -> ABCD
%! pnames = {"A", "B", "C", "D"};
%! [intMat, ~, ~] = FormulaParser ("y ~ A : (B : (C : D))", pnames);
%! assert (intMat, logical ([1, 1, 1, 1]));

%!test
%! ## Test 39: Interaction of High-Order Factorial
%! ## (A*B*C) : D -> (A+B+C+AB+AC+BC+ABC) : D
%! ## -> AD + BD + CD + ABD + ACD + BCD + ABCD
%! pnames = {"A", "B", "C", "D"};
%! [intMat, ~, ~] = FormulaParser ("y ~ (A*B*C):D", pnames);
%! assert (rows (intMat), 7);
%! ## Check ABCD
%! assert (ismember ([1, 1, 1, 1], intMat, "rows"));

%!test
%! ## Test 40: The Kitchen Sink
%! ## (A*B) + (C:D) + ((E+F):G)
%! ## -> A, B, AB
%! ## -> CD
%! ## -> EG, FG
%! pnames = {"A", "B", "C", "D", "E", "F", "G"};
%! [intMat, ~, ~] = FormulaParser ("y ~ (A*B) + (C:D) + ((E+F):G)", pnames);
%! assert (rows (intMat), 6);
%! ## Check AB (1100000)
%! assert (ismember ([1,1,0,0,0,0,0], intMat, "rows"));
%! ## Check EG (0000101)
%! assert (ismember ([0,0,0,0,1,0,1], intMat, "rows"));

%!error <FormulaParser: invalid syntax. Formula must contain '~'.> FormulaParser ("y x1 + x2", {"x1"})
%!error <FormulaParser: no predictor terms found.> FormulaParser ("y ~ ", {"x1"})
%!error <FormulaParser: predictor 'z' not found in PredictorNames.> FormulaParser ("y ~ x1 + z", {"x1"})
%!error <FormulaParser: invalid syntax> FormulaParser ("y ~ x1:", {"x1"})
%!error <FormulaParser: formula resulted in no valid terms.> FormulaParser ("y ~ -1", {"x1"})
