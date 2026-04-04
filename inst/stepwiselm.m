## Copyright (C) 2026 Jayant Chauhan <0001jayant@gmail.com>
##
## This file is part of the statistics package for GNU Octave.
##
## This program is free software; you can redistribute it and / or modify it under
##the terms of the GNU General Public License as published by the Free Software
##Foundation; either version 3 of the License, or (at your option) any later
##version.
##
## This program is distributed in the hope that it will be useful, but WITHOUT
## ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
## FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
## details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see < http: // www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn  {statistics} {@var{mdl} =} stepwiselm (@var{X}, @var{y})
## @deftypefnx {statistics} {@var{mdl} =} stepwiselm (@var{X}, @var{y}, @var{modelspec})
## @deftypefnx {statistics} {@var{mdl} =} stepwiselm (@dots{}, @var{Name}, @var{Value})
##
## Stepwise linear regression with automatic predictor selection.
##
## @code{stepwiselm} fits a linear regression model by iteratively adding
## and removing predictor terms.  Starting from a lower-bound model, the
## function adds the term that most improves the fit (while its addition
## p-value stays below @code{PEnter}) and removes the term whose p-value
## exceeds @code{PRemove}, repeating until the model stabilises or
## @code{NSteps} iterations are exhausted.
##
## @var{X} is an @var{n}-by-@var{p} numeric matrix of predictor variables.
## @var{y} is an @var{n}-by-1 numeric response vector.
##
## @var{modelspec} specifies the starting model and must be one of:
## @itemize
## @item @qcode{"constant"} -- intercept only (no predictors)
## @item @qcode{"linear"} (default) -- intercept plus all linear terms
## @item @qcode{"interactions"} -- linear terms plus all pairwise products
## @item @qcode{"full"} -- all main effects and higher-order interactions
## @end itemize
##
## @subheading Name-Value Options
##
## @table @asis
## @item @qcode{"Upper"}
## Upper-bound model: maximum complexity allowed.  Accepts the same
## keywords as @var{modelspec}.  Default is @qcode{"linear"}.
##
## @item @qcode{"Lower"}
## Lower-bound model: terms that are always kept in the model.
## Default is @qcode{"constant"} (intercept only).
##
## @item @qcode{"Criterion"}
## Selection criterion.  One of @qcode{"sse"} (p-value, default),
## @qcode{"aic"}, @qcode{"bic"}, @qcode{"rsquared"},
## @qcode{"adjrsquared"}.
##
## @item @qcode{"PEnter"}
## Threshold for adding a term.  For @qcode{"sse"}: maximum p-value for
## entry; default @code{0.05}, must be in @code{(0,1)}.
## For @qcode{"aic"}/@qcode{"bic"}: maximum AIC/BIC change allowed for
## entry; default @code{0} (add term only if criterion strictly decreases).
## For @qcode{"rsquared"}: minimum R@sup{2} increase required; default
## @code{0.1}.  For @qcode{"adjrsquared"}: default @code{0}.
## Must satisfy @code{PEnter < PRemove} for all criteria.
## @strong{Note:} for non-@qcode{"sse"} criteria, PEnter and PRemove are
## accepted but currently ignored internally; the greedy (best-move) rule
## is used instead.  This matches MATLAB's default-threshold behaviour.
##
## @item @qcode{"PRemove"}
## Threshold for removing a term.  For @qcode{"sse"}: minimum p-value to
## trigger removal; default @code{0.10}, must be in @code{(0,1)}, and
## must satisfy @code{PRemove > PEnter} strictly.
## For @qcode{"aic"}/@qcode{"bic"}: default @code{0.01}.
## For @qcode{"rsquared"}: default @code{0.05}.
## For @qcode{"adjrsquared"}: default @code{-0.05}.
## Must satisfy @code{PEnter < PRemove} for all criteria.
##
## @item @qcode{"NSteps"}
## Maximum number of stepwise iterations.  Default @code{Inf}.
##
## @item @qcode{"Verbose"}
## Display level: @code{0} (silent), @code{1} (final summary, default),
## @code{2} (per-step detail).
##