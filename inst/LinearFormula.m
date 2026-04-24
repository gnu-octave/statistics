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

classdef LinearFormula
## -*- texinfo -*-
## @deftp {statistics} LinearFormula
##
## Linear formula encapsulation class.
##
## A @code{LinearFormula} object encapsulates the internal formula
## representation for a @code{LinearModel} or @code{CompactLinearModel}.
## When displayed in the console, it shows only the human-readable
## formula string (e.g., @code{y ~ 1 + x1 + x2}).
##
## @strong{Properties}
##
## @multitable @columnfractions 0.30 0.70
## @headitem Property @tab Description
## @item @code{ResponseName} @tab Name of the response variable.
## @item @code{LinearPredictor} @tab Human-readable formula string.
## @item @code{Terms} @tab Binary terms matrix (n_terms x p).
## @item @code{HasIntercept} @tab Logical: true if model has intercept.
## @item @code{InModel} @tab Logical vector indicating which terms are
##   active.
## @item @code{CoefficientNames} @tab Cell array of coefficient names.
## @item @code{PredictorNames} @tab Cell array of predictor names.
## @item @code{VariableNames} @tab Cell array of all variable names.
## @end multitable
##
## @seealso{fitlm, LinearModel, CompactLinearModel}
## @end deftp

  properties (SetAccess = protected)
    ResponseName     = "";
    LinearPredictor  = "";
    Terms            = [];
    HasIntercept     = true;
    InModel          = [];
    CoefficientNames = {};
    PredictorNames   = {};
    VariableNames    = {};
  endproperties

  methods

    function obj = LinearFormula (s)
      ## LinearFormula  Constructor.
      ##
      ## obj = LinearFormula (s)
      ##
      ## S is a struct with formula fields.  If S is already a
      ## LinearFormula, it is returned as-is.

      if (nargin < 1)
        return;
      endif

      if (isa (s, "LinearFormula"))
        obj = s;
        return;
      endif

      if (! isstruct (s))
        error ("LinearFormula: constructor argument must be a struct.");
      endif

      if (isfield (s, "ResponseName"))
        obj.ResponseName = s.ResponseName;
      endif
      if (isfield (s, "LinearPredictor"))
        obj.LinearPredictor = s.LinearPredictor;
      endif
      if (isfield (s, "Terms"))
        obj.Terms = s.Terms;
      endif
      if (isfield (s, "HasIntercept"))
        obj.HasIntercept = s.HasIntercept;
      endif
      if (isfield (s, "InModel"))
        obj.InModel = s.InModel;
      endif
      if (isfield (s, "CoefficientNames"))
        obj.CoefficientNames = s.CoefficientNames;
      endif
      if (isfield (s, "PredictorNames"))
        obj.PredictorNames = s.PredictorNames;
      endif
      if (isfield (s, "VariableNames"))
        obj.VariableNames = s.VariableNames;
      endif
    endfunction

    function disp (obj)
      ## disp  Display only the formula string.
      ##
      ## Matches MATLAB's LinearFormula display behavior: shows only
      ## the human-readable formula (e.g., "y ~ 1 + x1 + x2").
      fprintf ("    %s ~ %s\n", obj.ResponseName, obj.LinearPredictor);
    endfunction

    function s = char (obj)
      ## char  Convert to character representation.
      s = sprintf ("%s ~ %s", obj.ResponseName, obj.LinearPredictor);
    endfunction

  endmethods

endclassdef
