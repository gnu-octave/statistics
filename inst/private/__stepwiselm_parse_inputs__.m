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
## @deftypefn {Private Function} {[@var{raw_X}, @var{raw_y}, @var{var_names}, @
##   @var{cat_flags}, @var{weights}, @var{excl_mask}, @var{modelspec}, @
##   @var{obs_names}, @var{resp_name}, @var{pred_names}, @var{has_intercept}, @
##   @var{raw_cols}, @var{lower_spec}, @var{upper_spec}, @var{p_enter}, @
##   @var{p_remove}, @var{criterion}, @var{n_steps}, @var{verbose}] =} @
##   __stepwiselm_parse_inputs__ (@var{varargin})
##
## Parse stepwiselm inputs.  Extracts the seven stepwise-specific name-value
## arguments (Lower, Upper, PEnter, PRemove, Criterion, NSteps, Verbose),
## validates them, then delegates the remaining arguments to
## @code{__fitlm_parse_inputs__} for standard fitlm parsing.
##
## Validation rules:
## @itemize
## @item @code{PEnter} must be strictly less than @code{PRemove}, regardless
## of criterion.  This matches MATLAB's behavior (even for AIC/BIC).
## @item @code{Criterion} must be one of @code{'sse'}, @code{'aic'},
## @code{'bic'}, @code{'rsquared'}, @code{'adjrsquared'}.
## @item @code{NSteps} must be a positive integer or @code{Inf}.
## @item @code{Verbose} must be 0, 1, or 2.
## @end itemize
##
## @seealso{stepwiselm, __fitlm_parse_inputs__}
## @end deftypefn

function [raw_X, raw_y, var_names, cat_flags, weights, excl_mask, ...
          modelspec, obs_names, resp_name, pred_names, has_intercept, ...
          raw_cols, lower_spec, upper_spec, p_enter, p_remove, ...
          criterion, n_steps, verbose] = ...
    __stepwiselm_parse_inputs__ (varargin)

  ## ── Initialize stepwise-specific defaults ─────────────────────────────────
  lower_spec = "constant";
  upper_spec = "interactions";
  p_enter    = 0.05;
  p_remove   = 0.10;
  criterion  = "sse";
  n_steps    = Inf;
  verbose    = 1;

  ## ── Extract stepwise NV args from varargin before passing to fitlm parser ─
  ## We walk through varargin and pull out stepwise-specific keys,
  ## leaving everything else for __fitlm_parse_inputs__.
  remaining = {};
  stepwise_keys = {"lower", "upper", "penter", "premove", ...
                   "criterion", "nsteps", "verbose"};

  i = 1;
  while (i <= numel (varargin))
    ## Non-char args (data matrices, table, numeric modelspec) pass through
    if (! ischar (varargin{i}))
      remaining{end+1} = varargin{i};
      i = i + 1;
      continue;
    endif

    key_lower = lower (varargin{i});

    ## Check if this is a stepwise-specific key
    if (any (strcmp (key_lower, stepwise_keys)))
      if (i + 1 > numel (varargin))
        error ("stepwiselm: name-value argument '%s' has no value.", ...
               varargin{i});
      endif
      val = varargin{i+1};
      switch (key_lower)
        case "lower"
          lower_spec = val;
        case "upper"
          upper_spec = val;
        case "penter"
          p_enter = double (val);
        case "premove"
          p_remove = double (val);
        case "criterion"
          criterion = lower (val);
        case "nsteps"
          n_steps = double (val);
        case "verbose"
          verbose = double (val);
      endswitch
      i = i + 2;
    else
      ## Not a stepwise key — pass through to fitlm parser
      remaining{end+1} = varargin{i};
      i = i + 1;
    endif
  endwhile

  ## ── Validate stepwise arguments ───────────────────────────────────────────

  ## PEnter < PRemove always required (matches MATLAB, even for AIC/BIC)
  if (p_enter >= p_remove)
    error ("stepwiselm: Cannot have PEnter (%g) >= PRemove (%g) for the %s criterion.", ...
           p_enter, p_remove, criterion);
  endif

  ## Criterion must be valid
  valid_criteria = {"sse", "aic", "bic", "rsquared", "adjrsquared"};
  if (! any (strcmp (criterion, valid_criteria)))
    error ("stepwiselm: Criterion must be one of: %s.", ...
           strjoin (valid_criteria, ", "));
  endif

  ## NSteps must be positive
  if (! (isinf (n_steps) || (n_steps > 0 && n_steps == round (n_steps))))
    error ("stepwiselm: NSteps must be a positive integer or Inf.");
  endif

  ## Verbose must be 0, 1, or 2
  if (! any (verbose == [0, 1, 2]))
    error ("stepwiselm: Verbose must be 0, 1, or 2.");
  endif

  ## ── Delegate to fitlm parser ──────────────────────────────────────────────
  [raw_X, raw_y, var_names, cat_flags, weights, excl_mask, ...
   modelspec, obs_names, resp_name, pred_names, has_intercept, raw_cols] = ...
    __fitlm_parse_inputs__ (remaining{:});

endfunction
