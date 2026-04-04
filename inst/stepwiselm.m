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