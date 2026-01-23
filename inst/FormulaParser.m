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

classdef FormulaParser
  ## -*- texinfo -*-
  ## @deftp {statistics} FormulaParser
  ##
  ## Utility class for parsing Wilkinson notation formulas and interaction specifications.
  ##
  ## This class provides static methods to convert formula strings (e.g., "y ~ x1 + x2")
  ## or interaction specifications into binary terms matrices required by regression
  ## and classification models.
  ##
  ## @end deftp

  methods (Static)

    ## -*- texinfo -*-
    ## @deftypefn {FormulaParser} {[@var{intMat}, @var{response}, @var{hasIntercept}] =} parseFormula (@var{formula}, @var{predictorNames})
    ## Parse a Wilkinson notation formula string into components.
    ##
    ## @var{formula} is a string (e.g., "y ~ A * B").
    ## @var{predictorNames} is a cell array of strings defining valid predictors.
    ##
    ## Output:
    ## @var{intMat}: Binary matrix (Terms x Predictors). 1 indicates variable presence.
    ## @var{response}: String name of the response variable.
    ## @var{hasIntercept}: Boolean. True unless "-1" is specified.
    ## @end deftypefn
    function [intMat, response, hasIntercept] = parseFormula (formula, predictorNames)
      intMat = [];
      response = "";
      hasIntercept = true; ## Default for Wilkinson notation

      ## 1. Syntax Check
      if (isempty (strfind (formula, '~')))
        error ("FormulaParser: invalid syntax. Formula must contain '~'.");
      endif

      ## 2. Split LHS (Response) and RHS (Predictors)
      parts = strsplit (formula, '~');
      if (numel (parts) < 2)
        error ("FormulaParser: no predictor terms in Formula.");
      endif

      response = strtrim (parts{1});
      rhs = strtrim (parts{2});

      if (isempty (rhs))
        error ("FormulaParser: no predictor terms in Formula.");
      endif

      ## 3. Handle Intercepts (+1 / -1)
      ## Check for explicit removal
      if (! isempty (strfind (rhs, '-1')) || ! isempty (strfind (rhs, '- 1')))
        hasIntercept = false;
        rhs = strrep (rhs, '- 1', '');
        rhs = strrep (rhs, '-1', '');
      endif

      ## Remove explicit inclusion (redundant but valid syntax)
      rhs = strrep (rhs, '+ 1', '');
      rhs = strrep (rhs, '+1', '');

      ## 4. Expand Terms (Handle + and *)
      ## First, split by '+' to get high-level blocks
      blocks = strtrim (strsplit (rhs, '+'));
      blocks(cellfun (@isempty, blocks)) = []; ## Clean up empty cells

      finalTerms = {};

      for i = 1:numel (blocks)
        block = blocks{i};

        ## Handle Cross Operator (*)
        ## A * B expands to: A + B + A:B
        ## FIX: Replaced 'contains' with '!isempty(strfind)' for Octave compatibility
        if (! isempty (strfind (block, '*')))
          factors = strtrim (strsplit (block, '*'));
          ## Generate all combinations (Power Set minus empty set)
          n_factors = numel (factors);
          combo_indices = ff2n (n_factors);
          ## Remove the 'all zeros' row
          combo_indices(1, :) = [];

          ## Construct terms from combinations
          for r = 1:rows (combo_indices)
            active_factors = factors(logical (combo_indices(r, :)));
            term_str = strjoin (active_factors, ':');
            finalTerms{end+1} = term_str;
          endfor
        else
          ## No *, just add the block (could be "A" or "A:B")
          finalTerms{end+1} = block;
        end
      endfor

      ## 5. Map Terms to Predictor Indices
      numPreds = numel (predictorNames);
      if (numPreds == 0)
         return;
      endif

      ## We use a map/hashtable approach or standard loop to build the matrix
      ## Rows = Terms, Cols = Predictors

      for i = 1:numel (finalTerms)
        currentTerm = finalTerms{i};

        ## Split by ':' to find interacting variables
        vars = strtrim (strsplit (currentTerm, ':'));

        ## Check for invalid syntax (e.g. "x1:") which creates empty strings
        if (any (cellfun (@isempty, vars)))
           error ("FormulaParser: invalid syntax in interaction term '%s'.", currentTerm);
        endif

        termRow = zeros (1, numPreds);

        for v = 1:numel (vars)
          idx = find (strcmp (predictorNames, vars{v}));
          if (isempty (idx))
             error ("FormulaParser: predictor '%s' not found in PredictorNames.", vars{v});
          endif
          termRow(idx) = 1;
        endfor

        ## Append to matrix
        intMat = [intMat; termRow];
      endfor

      ## 6. Cleanup (Unique rows)
      if (! isempty (intMat))
        ## Use 'stable' to preserve the order terms appeared in the formula (MATLAB compat)
        intMat = unique (intMat, 'rows', 'stable');
      else
        error ("FormulaParser: formula resulted in no valid terms.");
      endif

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {FormulaParser} {@var{intMat} =} parseInteractions (@var{interactions}, @var{numPredictors})
    ## Generate an interaction matrix from an interactions argument.
    ##
    ## @var{interactions} can be "all", a logical matrix, or a numeric scalar.
    ## @end deftypefn
    function intMat = parseInteractions (interactions, numPredictors)
      if (islogical (interactions))
        if (numPredictors != columns (interactions))
          error ("FormulaParser: columns in Interactions logical matrix must equal the number of predictors.");
        endif
        intMat = interactions;

      elseif (isnumeric (interactions) && isscalar (interactions))
        p = numPredictors;
        if (interactions < 0)
          error ("FormulaParser: Invalid interaction argument.");
        elseif (interactions > p * (p - 1) / 2)
          error ("FormulaParser: number of interaction terms requested is larger than all possible combinations.");
        endif

        allMat = ff2n (p);
        iterms = find (sum (allMat, 2) > 1);
        intMat = allMat(iterms, :);

      elseif (ischar (interactions) && strcmpi (interactions, "all"))
        p = numPredictors;
        allMat = ff2n (p);
        iterms = find (sum (allMat, 2) > 1);
        intMat = allMat(iterms, :);
      else
        error ("FormulaParser: Invalid interaction argument.");
      endif
    endfunction

  endmethods
endclassdef

%!test
%! pnames = {"x1", "x2", "x3"};
%! ## Test 1: Simple additive formula
%! [intMat, resp, hasInt] = FormulaParser.parseFormula ("y ~ x1 + x2", pnames);
%! expected = [1, 0, 0; 0, 1, 0];
%! assert (intMat, expected);
%! assert (resp, "y");
%! assert (hasInt, true);

%!test
%! pnames = {"A", "B"};
%! ## Test 2: Interaction term with colon and Intercept Removal
%! [intMat, resp, hasInt] = FormulaParser.parseFormula ("Cost ~ A + A:B - 1", pnames);
%! expected = [1, 0; 1, 1];
%! assert (intMat, expected);
%! assert (resp, "Cost");
%! assert (hasInt, false);

%!test
%! pnames = {"x1", "x2"};
%! ## Test 3: Whitespace handling
%! [intMat, ~, ~] = FormulaParser.parseFormula ("y ~ x1 + x1 : x2", pnames);
%! expected = [1, 0; 1, 1];
%! assert (intMat, expected);

%!test
%! ## Test 4: Parse Interactions "all"
%! intMat = FormulaParser.parseInteractions ("all", 3);
%! assert (size (intMat, 1), 4);
%! assert (columns (intMat), 3);

%!error <FormulaParser: invalid syntax. Formula must contain '~'.> FormulaParser.parseFormula ("y x1 + x2", {"x1"})
%!error <FormulaParser: no predictor terms in Formula.> FormulaParser.parseFormula ("y ~ ", {"x1"})
%!error <FormulaParser: predictor 'z' not found in PredictorNames.> FormulaParser.parseFormula ("y ~ x1 + z", {"x1"})
%!error <FormulaParser: invalid syntax in interaction term 'x1:'> FormulaParser.parseFormula ("y ~ x1:", {"x1"})
