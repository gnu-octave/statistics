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
    ## @deftypefn {FormulaParser} {@var{intMat} =} parseFormula (@var{formula}, @var{predictorNames})
    ## Parse a Wilkinson notation formula string into an interaction matrix.
    ##
    ## @var{formula} is a string like "y ~ x1 + x1:x2".
    ## @var{predictorNames} is a cell array of strings containing the names of valid predictors.
    ## @end deftypefn
    function intMat = parseFormula (formula, predictorNames)
      intMat = [];
      ## Check formula for syntax
      if (isempty (strfind (formula, '~')))
        error ("FormulaParser: invalid syntax in Formula. Missing '~'.");
      endif

      ## Split formula and keep predictor terms
      formulaParts = strsplit (formula, '~');

      ## Check there is some string after '~'
      if (numel (formulaParts) < 2)
        error ("FormulaParser: no predictor terms in Formula.");
      endif

      predictorString = strtrim (formulaParts{2});
      if (isempty (predictorString))
        error ("FormulaParser: no predictor terms in Formula.");
      endif

      ## Split additive terms (between + sign)
      aterms = strtrim (strsplit (predictorString, '+'));

      ## Process all terms
      for i = 1:numel (aterms)
        ## Find individual terms (string missing ':')
        if (isempty (strfind (aterms{i}, ':')))
          ## Search PredictorNames to associate with column in X
          sterms = strcmp (predictorNames, aterms(i));
          ## Append to interactions matrix
          intMat = [intMat; sterms];
        else
          ## Split interaction terms (string contains ':')
          ## Trim whitespace to handle "x1 : x2" correctly
          mterms = strtrim (strsplit (aterms{i}, ':'));

          ## Add each individual predictor to interaction term vector
          iterms = logical (zeros (1, numel (predictorNames)));
          for t = 1:numel (mterms)
            iterms = iterms | strcmp (predictorNames, mterms(t));
          endfor

          ## Check that all predictors have been identified
          if (sum (iterms) != numel (mterms))
            error ("FormulaParser: some predictors in the formula were not found in PredictorNames.");
          endif
          ## Append to interactions matrix
          intMat = [intMat; iterms];
        endif
      endfor

      ## Check that all terms have been identified
      if (! all (sum (intMat, 2) > 0))
        error ("FormulaParser: some terms have not been identified.");
      endif
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {FormulaParser} {@var{intMat} =} parseInteractions (@var{interactions}, @var{numPredictors})
    ## Generate an interaction matrix from an interactions argument.
    ##
    ## @var{interactions} can be "all", a logical matrix, or a numeric scalar.
    ## @var{numPredictors} is the integer number of predictor variables.
    ## @end deftypefn
    function intMat = parseInteractions (interactions, numPredictors)
      if (islogical (interactions))
        ## Check that interaction matrix corresponds to predictors
        if (numPredictors != columns (interactions))
          error ("FormulaParser: columns in Interactions logical matrix must equal the number of predictors.");
        endif
        intMat = interactions;

      elseif (isnumeric (interactions))
        ## Just check that the given number is not higher than
        ## p*(p-1)/2, where p is the number of predictors.
        p = numPredictors;
        if (interactions > p * (p - 1) / 2)
          error ("FormulaParser: number of interaction terms requested is larger than all possible combinations.");
        endif

        ## Generate full factorial combinations (binary)
        ## Use ff2n (Full Factorial 2-levels) to get true binary matrix
        allMat = ff2n (p);

        ## Only keep interaction terms (sum of row > 1)
        iterms = find (sum (allMat, 2) > 1);
        intMat = allMat(iterms, :);

      elseif (ischar (interactions) && strcmpi (interactions, "all"))
        p = numPredictors;
        ## Generate full factorial combinations (binary)
        allMat = ff2n (p);

        ## Only keep interaction terms (sum of row > 1)
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
%! intMat = FormulaParser.parseFormula ("y ~ x1 + x2", pnames);
%! expected = [1, 0, 0; 0, 1, 0];
%! assert (intMat, expected);

%!test
%! pnames = {"A", "B"};
%! ## Test 2: Interaction term with colon
%! intMat = FormulaParser.parseFormula ("y ~ A + A:B", pnames);
%! expected = [1, 0; 1, 1];
%! assert (intMat, expected);

%!test
%! pnames = {"x1", "x2"};
%! ## Test 3: Whitespace handling (robustness check)
%! intMat = FormulaParser.parseFormula ("y ~ x1 + x1 : x2", pnames);
%! expected = [1, 0; 1, 1];
%! assert (intMat, expected);

%!test
%! ## Test 4: Parse Interactions "all"
%! ## 3 vars -> Pairs: (1,2), (1,3), (2,3) -> 3 rows of 2s.
%! ## Plus triplet (1,2,3) -> 1 row of 3s.
%! ## ff2n(3) gives 8 rows.
%! ## Row sums: 0(1), 1(3), 2(3), 3(1).
%! ## Our logic keeps sum > 1. So 3+1 = 4 rows.
%! intMat = FormulaParser.parseInteractions ("all", 3);
%! assert (size (intMat, 1), 4);
%! assert (columns (intMat), 3);

%!error <FormulaParser: invalid syntax in Formula> FormulaParser.parseFormula ("y x1 + x2", {"x1"})
%!error <FormulaParser: no predictor terms> FormulaParser.parseFormula ("y ~ ", {"x1"})
%!error <FormulaParser: some predictors in the formula were not found> FormulaParser.parseFormula ("y ~ x1 + z", {"x1"})