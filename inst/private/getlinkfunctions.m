## Copyright (C) 2025 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn {Private Function} [@var{flink}, @var{dlink}, @var{ilink}, @var{errmsg}] = getlinkfunctions (@var{linkArg})
##
## Return the link function and its derivative and inverse base.
##
## @end deftypefn

function [flink, dlink, ilink, errmsg] = getlinkfunctions (linkArg)
  flink = dlink = ilink = [];
  errmsg = "";
  ## linkArg is a scalar structure
  if (isstruct (linkArg))
    if (! isscalar (linkArg))
      errmsg = "structure with custom link functions must be a scalar.";
      return;
    endif
    rf = {"Link", "Derivative", "Inverse"};
    if (! all (ismember (rf, fieldnames (linkArg))))
      errmsg = ["structure with custom link functions requires", ...
                " the fields 'Link', 'Derivative', and 'Inverse'."];
      return;
    endif
    if (ischar (linkArg.Link) && ! isempty (which (linkArg.Link)))
      flink = @(mu) feval (linkArg.Link, mu);
    elseif (isa (linkArg.Link, "function_handle"))
      flink = linkArg.Link;
    else
      errmsg = ["bad 'Link' function in custom link function structure."];
      return;
    endif
    errmsg = testLinkFunction (flink, "Link");
    if (! isempty (errmsg))
      return;
    endif
    if (ischar (linkArg.Derivative) && ! isempty (which (linkArg.Derivative)))
      dlink = @(mu) feval (linkArg.Derivative, mu);
    elseif (isa (linkArg.Derivative, "function_handle"))
      dlink = linkArg.Derivative;
    else
      errmsg = ["bad 'Derivative' function in custom link function structure."];
      return;
    endif
    errmsg = testLinkFunction (dlink, "Derivative");
    if (! isempty (errmsg))
      return;
    endif
    if (ischar (linkArg.Inverse) && ! isempty (which (linkArg.Inverse)))
      ilink = @(mu) feval (linkArg.Inverse, mu);
    elseif (isa (linkArg.Inverse, "function_handle"))
      ilink = linkArg.Inverse;
    else
      errmsg = ["bad 'Inverse' function in custom link function structure."];
      return;
    endif
    errmsg = testLinkFunction (ilink, "Inverse");
    if (! isempty (errmsg))
      return;
    endif

  ## linkArg is a cell array
  elseif (iscell (linkArg))
    if (numel (linkArg) != 3)
      errmsg = "cell array with custom link functions must have three elements.";
      return;
    endif
    if (isa (linkArg{1}, "function_handle"))
      flink = linkArg{1};
    else
      errmsg = ["bad 'Link' function in custom link function cell array."];
      return;
    endif
    errmsg = testLinkFunction (flink, "Link");
    if (! isempty (errmsg))
      return;
    endif
    if (isa (linkArg{2}, "function_handle"))
      dlink = linkArg{2};
    else
      errmsg = ["bad 'Derivative' function in custom link function cell array."];
      return;
    endif
    errmsg = testLinkFunction (dlink, "Derivative");
    if (! isempty (errmsg))
      return;
    endif
    if (isa (linkArg{3}, "function_handle"))
      ilink = linkArg{3};
    else
      errmsg = ["bad 'Inverse' function in custom link function cell array."];
      return;
    endif
    errmsg = testLinkFunction (ilink, "Inverse");
    if (! isempty (errmsg))
      return;
    endif

  ## linkArg is a scalar value
  elseif (isnumeric (linkArg))
    if (! (isscalar (linkArg)  && isfinite (linkArg) && isreal (linkArg)))
       errmsg = ["numeric input for custom link function", ...
                 " must be a finite real scalar value."];
      return;
    endif
    flink = @(mu) mu .^ linkArg;
    dlink = @(mu) linkArg .* mu .^ (linkArg - 1);
    ilink = @(eta) eta .^ (1 / linkArg);

  ## linkArg is character vector
  elseif (ischar (linkArg))
    if (! isvector (linkArg))
      errmsg = "canonical link function name must be a character vector.";
      return;
    endif
    supported_link_functions = {"identity", "log", "logit", "probit", ...
                                "loglog", "comploglog", "reciprocal"};
    if (! any (strcmpi (linkArg, supported_link_functions)))
      errmsg = sprintf ("canonical link function '%s' is not supported.", ...
                        linkArg);
      return;
    endif
    ## Select a canonical link function
    switch (linkArg)
      case "identity"
        flink = @(mu) mu;
        dlink = @(mu) 1;
        ilink = @(eta) eta;
      case "log"
        flink = @(mu) log (mu);
        dlink = @(mu)  1 ./ (mu);
        ilink = @(eta) exp (eta);
      case "logit"
        flink = @(mu) log (mu ./ (1 - mu));
        dlink = @(mu) 1 ./ (mu .* (1 - mu));
        ilink = @(eta) 1 ./ (1 + exp (- eta));
      case "probit"
        flink = @(mu) norminv (mu);
        dlink = @(mu) 1 ./ normpdf (norminv (mu));
        ilink = @(eta) normcdf (eta);
      case "loglog"
        flink = @(mu) log (- log (mu));
        dlink = @(mu)  1 ./ (mu .* log (mu));
        ilink = @(eta) exp (- exp (eta));
      case "comploglog"
        flink = @(mu) log (- log1p (- mu));
        dlink = @(mu) 1 ./ - ((1 - mu) .* log1p (- mu));
        ilink = @(eta) -expm1 (- exp (eta));
      case "reciprocal"
        flink = @(mu) 1 ./ mu;
        dlink = @(mu) -1 ./ (mu .^ 2);
        ilink = @(eta) 1 ./ (eta);
    endswitch
  else
    errmsg = "invalid value for custom link function.";
  endif
endfunction

function errmsg = testLinkFunction (flink, linkname);
  errmsg = "";
  testInput = [1; 2; 3; 4; 5];
  try
    testOutput = flink (testInput);
    if (! isequal (size (testInput), size (testOutput)))
      errmsg = sprintf (["custom '%s' function must return an output", ...
                         " of the same size as input."], linkname);
    endif
  catch
    errmsg = sprintf ("invalid custom '%s' function.", linkname);
  end_try_catch
endfunction
