## Copyright (C) 2026 Jayant Chauhan <0001jayant@gmail.com>

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

classdef FormulaParser < handle
    properties
        FormulaString
        VariableNames
        Terms          % The binary matrix (Model Terms)
        ResponseName   % The name of the Y variable
        PredictorNames % The names of the X variables used
        HasIntercept   % Boolean
    end

    methods
        function obj = FormulaParser(formulaStr, varNames)
             obj.FormulaString = formulaStr;
             obj.VariableNames = varNames;
             obj.parse();
        end

        function parse(obj)
            % 1. Split by '~' to get Response and Predictors
            % 2. Parse ResponseName (trim whitespace)
            % 3. Check for '-1' or '0' to determine HasIntercept
            % 4. Split Predictors by '+'
            % 5. Handle ':' interactions
            % 6. Map strings to indices in VariableNames to build obj.Terms
        end
    end
end
