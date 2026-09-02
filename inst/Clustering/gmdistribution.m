## Copyright (C) 2015 Lachlan Andrew <lachlanbis@gmail.com>
## Copyright (C) 2018-2020 John Donoghue <john.donoghue@ieee.org>
## Copyright (C) 2022 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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

classdef gmdistribution
   ## -*- texinfo -*-
   ## @deftypefn  {statistics} {@var{GMdist} =} gmdistribution (@var{mu}, @var{Sigma})
   ## @deftypefnx {statistics} {@var{GMdist} =} gmdistribution (@var{mu}, @var{Sigma}, @var{p})
   ## @deftypefnx {statistics} {@var{GMdist} =} gmdistribution (@var{mu}, @var{Sigma}, @var{p}, @var{extra})
   ##
   ## Create an object of the gmdistribution class which represents a Gaussian
   ## mixture model with k components of n-dimensional Gaussians.
   ##
   ## Input @var{mu} is a k-by-n matrix specifying the n-dimensional mean of
   ## each of the k components of the distribution.
   ##
   ## Input @var{Sigma} is an array that specifies the variances of the
   ## distributions, in one of four forms depending on its dimension.
   ## @itemize
   ##   @item n-by-n-by-k: Slice @var{Sigma}(:,:,i) is the variance of the
   ##         i'th component
   ##   @item 1-by-n-by-k: Slice diag(@var{Sigma}(1,:,i)) is the variance of the
   ##         i'th component
   ##   @item n-by-n: @var{Sigma} is the variance of every component
   ##   @item 1-by-n-by-k: Slice diag(@var{Sigma}) is the variance of every
   ##         component
   ## @end itemize
   ##
   ## If @var{p} is specified, it is a vector of length k specifying the
   ## proportion of each component.  If it is omitted or empty, each component
   ## has an equal proportion.
   ##
   ## Input @var{extra} is used by fitgmdist to indicate the parameters of the
   ## fitting process.
   ## @seealso{fitgmdist}
   ## @end deftypefn

   properties

      ## -*- texinfo -*-
      ## @deftp {gmdistribution} {property} mu
      ##
      ## Component means
      ##
      ## A k-by-d matrix holding the mean of each of the k components, one
      ## per row, where d is the number of variables.  This property is
      ## read-only.
      ##
      ## @end deftp
      mu

      ## -*- texinfo -*-
      ## @deftp {gmdistribution} {property} Sigma
      ##
      ## Component covariances
      ##
      ## The covariances of the components, in the form they were given in.
      ## A full covariance per component is d-by-d-by-k, a diagonal one per
      ## component 1-by-d-by-k, a full covariance shared by every component
      ## d-by-d, and a shared diagonal one 1-by-d.  This property is
      ## read-only.
      ##
      ## @end deftp
      Sigma

      ## -*- texinfo -*-
      ## @deftp {gmdistribution} {property} ComponentProportion
      ##
      ## Mixing proportions
      ##
      ## A 1-by-k row vector holding the proportion of each component.  The
      ## proportions are scaled to sum to 1, and are all equal when none was
      ## given.  This property is read-only.
      ##
      ## @end deftp
      ComponentProportion

      ## -*- texinfo -*-
      ## @deftp {gmdistribution} {property} DistributionName
      ##
      ## Name of the distribution
      ##
      ## The character vector @qcode{'gaussian mixture distribution'}.  This
      ## property is read-only.
      ##
      ## @end deftp
      DistributionName

      ## -*- texinfo -*-
      ## @deftp {gmdistribution} {property} NumComponents
      ##
      ## Number of mixture components
      ##
      ## A positive integer, the number of rows of @code{mu}.  This property
      ## is read-only.
      ##
      ## @end deftp
      NumComponents

      ## -*- texinfo -*-
      ## @deftp {gmdistribution} {property} NumVariables
      ##
      ## Number of variables
      ##
      ## A positive integer, the dimension of each component, which is the
      ## number of columns of @code{mu}.  This property is read-only.
      ##
      ## @end deftp
      NumVariables

      ## -*- texinfo -*-
      ## @deftp {gmdistribution} {property} CovarianceType
      ##
      ## Form of the component covariances
      ##
      ## A character vector, @qcode{'diagonal'} when each covariance was
      ## given as a row of variances and @qcode{'full'} when it was given as
      ## a matrix.  This property is read-only.
      ##
      ## @end deftp
      CovarianceType

      ## -*- texinfo -*-
      ## @deftp {gmdistribution} {property} SharedCovariance
      ##
      ## Whether the components share one covariance
      ##
      ## A logical scalar, true when a single covariance was given for every
      ## component and false when one was given per component.  This property
      ## is read-only.
      ##
      ## @end deftp
      SharedCovariance

      ## -*- texinfo -*-
      ## @deftp {gmdistribution} {property} AIC
      ##
      ## Akaike information criterion
      ##
      ## A scalar, twice the negative log-likelihood plus twice the number of
      ## estimated parameters.  It is empty unless the object came from
      ## @code{fitgmdist}.  This property is read-only.
      ##
      ## @end deftp
      AIC

      ## -*- texinfo -*-
      ## @deftp {gmdistribution} {property} BIC
      ##
      ## Bayesian information criterion
      ##
      ## A scalar, twice the negative log-likelihood plus the number of
      ## estimated parameters times the log of the number of observations.
      ## It is empty unless the object came from @code{fitgmdist}.  This
      ## property is read-only.
      ##
      ## @end deftp
      BIC

      ## -*- texinfo -*-
      ## @deftp {gmdistribution} {property} Converged
      ##
      ## Whether the fit converged
      ##
      ## A logical scalar, true when the fit reached its tolerance within the
      ## iteration limit.  It is empty unless the object came from
      ## @code{fitgmdist}.  This property is read-only.
      ##
      ## @end deftp
      Converged

      ## -*- texinfo -*-
      ## @deftp {gmdistribution} {property} NegativeLogLikelihood
      ##
      ## Negative log-likelihood at the fit
      ##
      ## A scalar, the negative of the log-likelihood of the data under the
      ## fitted mixture.  It is empty unless the object came from
      ## @code{fitgmdist}.  This property is read-only.
      ##
      ## @end deftp
      NegativeLogLikelihood

      ## -*- texinfo -*-
      ## @deftp {gmdistribution} {property} NumIterations
      ##
      ## Iterations the fit took
      ##
      ## A positive integer.  It is empty unless the object came from
      ## @code{fitgmdist}.  This property is read-only.
      ##
      ## @end deftp
      NumIterations

      ## -*- texinfo -*-
      ## @deftp {gmdistribution} {property} RegularizationValue
      ##
      ## Regularization added to the covariance diagonal
      ##
      ## A nonnegative scalar added to the diagonal of each covariance to
      ## keep it positive definite.  It is empty unless the object came from
      ## @code{fitgmdist}.  This property is read-only.
      ##
      ## @end deftp
      RegularizationValue

   endproperties

   ## Kept out of the documented surface.  MATLAB carries the same alias and
   ## hides it, along with eight other legacy names this class does not have.
   properties (Hidden)

      ## A second name for NegativeLogLikelihood, holding the same value.
      NlogL

   endproperties

   properties(Access = private)
      DiagonalCovariance        ## bool summary of "CovarianceType"
   endproperties

   methods
      ## -*- texinfo -*-
      ## @deftypefn  {gmdistribution} {@var{obj} =} gmdistribution (@var{mu}, @var{Sigma})
      ## @deftypefnx {gmdistribution} {@var{obj} =} gmdistribution (@var{mu}, @var{Sigma}, @var{p})
      ## @deftypefnx {gmdistribution} {@var{obj} =} gmdistribution (@var{mu}, @var{Sigma}, @var{p}, @var{extra})
      ##
      ## Create a Gaussian mixture distribution.
      ##
      ## @var{mu} is a k-by-d matrix holding the mean of each of the k
      ## components, one per row, where d is the number of variables.
      ##
      ## @var{Sigma} holds the covariances in one of four forms.  A full
      ## covariance per component is d-by-d-by-k and a diagonal one per
      ## component is 1-by-d-by-k, while a single d-by-d matrix or a single
      ## 1-by-d row of variances is shared by every component.  The form given
      ## sets @code{CovarianceType} and @code{SharedCovariance}.
      ##
      ## @var{p} is a vector of k mixing proportions, scaled to sum to 1.  A
      ## proportion may not be negative and they may not all be zero.  When
      ## @var{p} is omitted or empty the components are equally weighted.
      ##
      ## @var{extra} carries the results of a fit and is passed by
      ## @code{fitgmdist}.  It fills @code{AIC}, @code{BIC},
      ## @code{Converged}, @code{NegativeLogLikelihood},
      ## @code{NumIterations} and @code{RegularizationValue}, which stay
      ## empty for an object built by hand.
      ##
      ## @end deftypefn
      function obj = gmdistribution (mu,sigma,p = [],extra = [])
        obj.DistributionName = 'gaussian mixture distribution';
        obj.mu = mu;
        obj.Sigma = sigma;
        obj.NumComponents = rows (mu);
        obj.NumVariables = columns (mu);
        if (isempty (p))
          obj.ComponentProportion = ones (1,obj.NumComponents) / ...
                                            obj.NumComponents;
        else
          if any (p < 0)
            error ("gmdistribution: component weights must be non-negative");
          endif
          s = sum (p);
          if (s == 0)
            error ("gmdistribution: component weights must not be all zero");
          elseif (s != 1)
            p = p / s;
          endif
          obj.ComponentProportion = p(:)';
        endif
        if (length (size (sigma)) == 3)
          obj.SharedCovariance = false;
        else
          obj.SharedCovariance = true;
        endif
        if (rows (sigma) == 1 && columns (mu) > 1)
          obj.DiagonalCovariance = true;
          obj.CovarianceType = 'diagonal';
        else
          obj.DiagonalCovariance = false;       ## full
          obj.CovarianceType = 'full';
        endif

        if (! isempty (extra))
          obj.AIC                   = extra.AIC;
          obj.BIC                   = extra.BIC;
          obj.Converged             = extra.Converged;
          obj.NegativeLogLikelihood = extra.NegativeLogLikelihood;
          obj.NlogL                 = extra.NegativeLogLikelihood;
          obj.NumIterations         = extra.NumIterations;
          obj.RegularizationValue   = extra.RegularizationValue;
        endif
      endfunction

      ## -*- texinfo -*-
      ## @deftypefn {gmdistribution} {@var{c} =} cdf (@var{obj}, @var{X})
      ##
      ## Cumulative distribution function of a Gaussian mixture distribution.
      ##
      ## @var{X} is an n-by-d matrix of points at which the distribution is
      ## evaluated, one point per row, where d is the number of variables of
      ## @var{obj}.  @var{c} is an n-by-1 vector holding the value of the
      ## cumulative distribution function at each of them, the components'
      ## cumulative distributions summed with the mixing proportions in
      ## @code{ComponentProportion} as weights.
      ##
      ## @end deftypefn
      function c = cdf (obj, X)
        X = checkX (obj, X, 'cdf');
        p_x_l = zeros (rows (X), obj.NumComponents);
        if (obj.SharedCovariance)
          if (obj.DiagonalCovariance)
            sig = diag (obj.Sigma);
          else
            sig = obj.Sigma;
          endif
        endif
        for i = 1:obj.NumComponents
          if (! obj.SharedCovariance)
            if (obj.DiagonalCovariance)
              sig = diag (obj.Sigma(:,:,i));
            else
              sig = obj.Sigma(:,:,i);
            endif
          endif
          p_x_l(:,i) = mvncdf (X,obj.mu(i,:),sig)*obj.ComponentProportion(i);
        endfor
        c = sum (p_x_l, 2);
      endfunction

      ## -*- texinfo -*-
      ## @deftypefn  {gmdistribution} {@var{idx} =} cluster (@var{obj}, @var{X})
      ## @deftypefnx {gmdistribution} {[@var{idx}, @var{nlogl}] =} cluster (@var{obj}, @var{X})
      ## @deftypefnx {gmdistribution} {[@var{idx}, @var{nlogl}, @var{P}] =} cluster (@var{obj}, @var{X})
      ## @deftypefnx {gmdistribution} {[@var{idx}, @var{nlogl}, @var{P}, @var{logpdf}] =} cluster (@var{obj}, @var{X})
      ## @deftypefnx {gmdistribution} {[@var{idx}, @var{nlogl}, @var{P}, @var{logpdf}, @var{M}] =} cluster (@var{obj}, @var{X})
      ##
      ## Assign each observation to a mixture component.
      ##
      ## @var{X} is an n-by-d matrix of observations, one per row, where d is
      ## the number of variables of @var{obj}.  Each is assigned to the
      ## component under which it is most probable.
      ##
      ## @var{idx} is an n-by-1 vector of component indices.  @var{nlogl} is
      ## the negative log-likelihood of @var{X} under the mixture.  @var{P} is
      ## an n-by-k matrix of posterior probabilities, one column per component,
      ## whose rows sum to one.  @var{logpdf} is an n-by-1 vector holding the
      ## logarithm of the mixture density at each observation.  @var{M} is an
      ## n-by-k matrix of squared Mahalanobis distances from each observation
      ## to each component mean.
      ##
      ## @end deftypefn
      function [idx, nlogl, P, logpdf, M] = cluster (obj,X)
        X = checkX (obj, X, 'cluster');
        [p_x_l, M] = componentProb (obj, X);
        [~, idx] = max (p_x_l, [], 2);
        if (nargout >= 2)
          PDF = sum (p_x_l, 2);
          logpdf = log (PDF);
          nlogl = -sum (logpdf);
          if (nargout >= 3)
            P = bsxfun (@rdivide, p_x_l, PDF);
          endif
        endif
      endfunction


      ## -*- texinfo -*-
      ## @deftypefn {gmdistribution} {@var{D} =} mahal (@var{obj}, @var{X})
      ##
      ## Squared Mahalanobis distance to each mixture component.
      ##
      ## @var{X} is an n-by-d matrix of observations, one per row, where d is
      ## the number of variables of @var{obj}.  @var{D} is an n-by-k matrix
      ## holding the squared Mahalanobis distance from each observation to the
      ## mean of each of the k components, measured in the covariance of the
      ## component it is taken to.
      ##
      ## @end deftypefn
      function D = mahal (obj,X)
        X = checkX (obj, X, 'mahal');
        [~, D] = componentProb (obj,X);
      endfunction

      ## -*- texinfo -*-
      ## @deftypefn {gmdistribution} {@var{c} =} pdf (@var{obj}, @var{X})
      ##
      ## Probability density function of a Gaussian mixture distribution.
      ##
      ## @var{X} is an n-by-d matrix of points at which the density is
      ## evaluated, one point per row, where d is the number of variables of
      ## @var{obj}.  @var{c} is an n-by-1 vector holding the density at each of
      ## them, the components' densities summed with the mixing proportions in
      ## @code{ComponentProportion} as weights.
      ##
      ## @end deftypefn
      function c = pdf (obj,X)
        X = checkX (obj, X, 'pdf');
        p_x_l = componentProb (obj, X);
        c = sum (p_x_l, 2);
      endfunction

      ## -*- texinfo -*-
      ## @deftypefn {gmdistribution} {@var{c} =} posterior (@var{obj}, @var{X})
      ##
      ## Posterior probability of each mixture component.
      ##
      ## @var{X} is an n-by-d matrix of observations, one per row, where d is
      ## the number of variables of @var{obj}.  @var{c} is an n-by-k matrix
      ## whose (i,j) element is the probability that observation i was drawn
      ## from component j, so its rows sum to one.
      ##
      ## @end deftypefn
      function c = posterior (obj,X)
        X = checkX (obj, X, 'posterior');
        p_x_l = componentProb (obj, X);
        c = bsxfun (@rdivide, p_x_l, sum (p_x_l, 2));
      endfunction

      ## -*- texinfo -*-
      ## @deftypefn  {gmdistribution} {@var{c} =} random (@var{obj})
      ## @deftypefnx {gmdistribution} {@var{c} =} random (@var{obj}, @var{n})
      ##
      ## Random numbers from a Gaussian mixture distribution.
      ##
      ## @var{n} is the number of observations to draw and defaults to 1.
      ## @var{c} is an @var{n}-by-d matrix holding one observation per row,
      ## where d is the number of variables of @var{obj}.  Each row is drawn
      ## from a component chosen with probability @code{ComponentProportion}.
      ##
      ## @end deftypefn
      function c = random (obj,n)
        if nargin == 1
          n = 1;
        endif
        c = zeros (n, obj.NumVariables);
        classes = randsample (obj.NumComponents, n, true, ...
                              obj.ComponentProportion);
        if (obj.SharedCovariance)
          if (obj.DiagonalCovariance)
            sig = diag (obj.Sigma);
          else
            sig = obj.Sigma;
          endif
        endif
        for i = 1:obj.NumComponents
          idx = (classes == i);
          k = sum (idx);
          if (k > 0)
            if (! obj.SharedCovariance)
              if (obj.DiagonalCovariance)
                sig = diag (obj.Sigma(:,:,i));
              else
                sig = obj.Sigma(:,:,i);
              endif
            endif
            # [sig] forces [sig] not to have class "diagonal",
            # since mvnrnd uses automatic broadcast,
            # which fails on structured matrices
            c(idx,:) = mvnrnd ([obj.mu(i,:)], [sig], k);
          endif
        endfor
      endfunction
    endmethods

    ########################################
    methods (Hidden)
      function c = disp (obj)
        msg = ['Gaussian mixture distribution with %d ', ...
               "components in %d dimension(s)\n"];
        fprintf (msg, obj.NumComponents, columns (obj.mu));
        for i = 1:obj.NumComponents
          fprintf ("Clust %d: weight %d\n\tMean: ", ...
                          i, obj.ComponentProportion(i));
          fprintf ("%g ", obj.mu(i,:));
          fprintf ("\n");
          if (! obj.SharedCovariance)
            fprintf ("\tVariance:");
            if (! obj.DiagonalCovariance)
              if (columns (obj.mu) > 1)
                fprintf ("\n");
              endif
              disp (squeeze (obj.Sigma(:,:,i)))
            else
              fprintf (" diag(");
              fprintf ("%g ", obj.Sigma(:,:,i));
              fprintf (")\n");
            endif
          endif
        endfor
        if (obj.SharedCovariance)
          fprintf ("Shared variance\n");
          if (! obj.DiagonalCovariance)
            obj.Sigma
          else
            fprintf (" diag(");
            fprintf ("%g ", obj.Sigma);
            fprintf (")\n");
          endif
        endif
        if (! isempty (obj.AIC))
          fprintf ("AIC=%g BIC=%g NLogL=%g Iter=%d Cged=%d Reg=%g\n", ...
                    obj.AIC, obj.BIC, obj.NegativeLogLikelihood, ...
                    obj.NumIterations, obj.Converged, obj.RegularizationValue);
        endif
      endfunction

      function c = display (obj)
        disp (obj);
      endfunction
    endmethods

    ########################################
    methods(Static)
      ## -*- texinfo -*-
      ## @deftypefn  {gmdistribution} {@var{obj} =} fit (@var{X}, @var{k})
      ## @deftypefnx {gmdistribution} {@var{obj} =} fit (@var{X}, @var{k}, @var{Name}, @var{Value})
      ##
      ## Fit a Gaussian mixture distribution to data.
      ##
      ## @var{X} is an n-by-d matrix of observations, one per row, and @var{k}
      ## is the number of components to fit.  Any @var{Name}-@var{Value} pair
      ## accepted by @code{fitgmdist} may follow, which is the function this
      ## method calls and where the options are documented.  @var{obj} is the
      ## fitted @code{gmdistribution} object.
      ##
      ## @end deftypefn
      function c = fit (X, k, varargin)
        c = fitgmdist (X, k, varargin{:});
      endfunction
    endmethods

    ########################################
    methods(Access = private)
      ## Probability density of (row of) X *and* component l
      ## Second argument is an array of the Mahalanobis distances
      function [p_x_l, M] = componentProb (obj, X)
        M     = zeros (rows (X), obj.NumComponents);
        dets  = zeros (1, obj.NumComponents);   % sqrt(determinant)
        if (obj.SharedCovariance)
          if (obj.DiagonalCovariance)
            r = diag (sqrt (obj.Sigma));
          else
            r = chol (obj.Sigma);
          endif
        endif
        for i = 1:obj.NumComponents
          dev = bsxfun (@minus, X, obj.mu(i,:));
          if (! obj.SharedCovariance)
            if (obj.DiagonalCovariance)
              r = diag (sqrt (obj.Sigma(:,:,i)));
            else
              r = chol (obj.Sigma(:,:,i));
            endif
          endif
          M(:,i) = sumsq (dev / r, 2);
          dets(i) = prod (diag (r));
        endfor
        p_x_l = exp (-M/2);
        coeff = obj.ComponentProportion ./ ...
                    ((2 * pi) ^ (obj.NumVariables / 2) .* dets);
        p_x_l = bsxfun (@times, p_x_l, coeff);
      endfunction

      ########################################
      ## Check format of argument X
      function X = checkX (obj, X, name)
        if (columns (X) != obj.NumVariables)
          if (columns (X) == 1 && rows (X) == obj.NumVariables)
            X = X';
          else
            error ("gmdistribution.%s: X has %d columns instead of %d\n", ...
                   name, columns (X), obj.NumVariables);
          endif
        endif
      endfunction
    endmethods
endclassdef

%!test
%! mu = eye (2);
%! Sigma = eye (2);
%! GM = gmdistribution (mu, Sigma);
%! density = GM.pdf ([0 0; 1 1]);
%! assert_equal (density(1) - density(2), 0, 1e-6);
%!
%! [idx, nlogl, P, logpdf,M] = cluster (GM, eye (2));
%! assert_equal (idx, [1; 2]);
%! [idx2,nlogl2,P2,logpdf2] = GM.cluster (eye (2));
%! assert_equal (nlogl - nlogl2, 0, 1e-6);
%! [idx3,nlogl3,P3] = cluster (GM, eye (2));
%! assert_equal (P - P3, zeros (2), 1e-6);
%! [idx4,nlogl4] = cluster (GM, eye (2));
%! assert_equal (size (nlogl4), [1 1]);
%! idx5 = cluster (GM, eye (2));
%! assert_equal (idx - idx5, zeros (2,1));
%!
%! D = GM.mahal ([1;0]);
%! assert_equal (D - M(1,:), zeros (1,2), 1e-6);
%!
%! P = GM.posterior ([0 1]);
%! assert_equal (P - P2(2,:), zeros (1,2), 1e-6);
%!
%! R = GM.random(20);
%! assert_equal (size (R), [20, 2]);
%!
%! R = GM.random();
%! assert_equal (size (R), [1, 2]);

