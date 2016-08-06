## Copyright (C) 2015 Lachlan Andrew <lachlanbis@gmail.com>
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
##
## -*- texinfo -*-
## @deftypefn {Function File} {@var{GMdist} =} gmdistribution (@var{mu}, @var{Sigma})
## @deftypefnx {Function File} {@var{GMdist} =} gmdistribution (@var{mu}, @var{Sigma}, @var{p})
## @deftypefnx {Function File} {@var{GMdist} =} gmdistribution (@var{mu}, @var{Sigma}, @var{p}, @var{extra})
## Create an object of the  gmdistribution  class which represents a Gaussian
## mixture model with k components of n-dimensional Gaussians.
##
## Input @var{mu} is a k-by-n matrix specifying the n-dimensional mean of each
## of the k components of the distribution.
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
## If @var{p} is specified, it is a vector of length k specifying the proportion
## of each component.  If it is omitted or empty, each component has an equal
## proportion.
##
## Input @var{extra} is used by fitgmdist to indicate the parameters of the
## fitting process.
## @seealso{fitgmdist}
## @end deftypefn
classdef gmdistribution
   properties
      mu                        ## means
      Sigma                     ## covariances
      ComponentProportion       ## mixing proportions
      DistributionName          ## "gaussian mixture distribution"
      NumComponents             ## Number of mixture components
      NumVariables              ## Dimension d of each Gaussian component

      CovarianceType            ## 'diagonal' if DiagonalCovariance, 'full' othw
      SharedCovariance          ## true if all components have equal covariance

      ## Set by a call to gmdistribution.fit  or  fitgmdist
      AIC                       ## Akaike Information Criterion
      BIC                       ## Bayes Information Criterion
      Converged                 ## true  if algorithm converged by MaxIter
      NegativeLogLikelihood     ## Negative of log-likelihood
      NlogL                     ## Negative of log-likelihood
      NumIterations             ## Number of iterations
      RegularizationValue       ## const added to diag of cov to make +ve def
   endproperties

   properties (Access = private)
      DiagonalCovariance        ## bool summary of "CovarianceType"
   endproperties

   methods
      ########################################
      ## Constructor
      function obj = gmdistribution (mu,sigma,p = [],extra = [])
        obj.DistributionName = "gaussian mixture distribution";
        obj.mu = mu;
        obj.Sigma = sigma;
        obj.NumComponents = rows (mu);
        obj.NumVariables = columns (mu);
        if (isempty (p))
          obj.ComponentProportion = ones (1,obj.NumComponents) / obj.NumComponents;
        else
          if any (p < 0)
            error ("gmmdistribution: component weights must be non-negative");
          endif
          s = sum(p);
          if (s == 0)
            error ("gmmdistribution: component weights must not be all zero");
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

        if (!isempty (extra))
          obj.AIC                   = extra.AIC;
          obj.BIC                   = extra.BIC;
          obj.Converged             = extra.Converged;
          obj.NegativeLogLikelihood = extra.NegativeLogLikelihood;
          obj.NlogL                 = extra.NegativeLogLikelihood;
          obj.NumIterations         = extra.NumIterations;
          obj.RegularizationValue   = extra.RegularizationValue;
        endif
      endfunction

        ########################################
        ## Cumulative distribution function for Gaussian mixture distribution
      function c = cdf (obj, X)
        X = checkX (obj, X, "cdf");
        p_x_l = zeros (rows (X), obj.NumComponents);
        if (obj.SharedCovariance)
          if (obj.DiagonalCovariance)
            sig = diag (obj.Sigma);
          else
            sig = obj.Sigma;
          endif
        endif
        for i = 1:obj.NumComponents
          if (!obj.SharedCovariance)
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

        ########################################
        ## Construct clusters from Gaussian mixture distribution
        ##
      function [idx,nlogl,P,logpdf,M] = cluster (obj,X)
        X = checkX (obj, X, "cluster");
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

        ########################################
        ## Display Gaussian mixture distribution object
      function c = disp (obj)
          fprintf("Gaussian mixture distribution with %d components in %d dimension(s)\n", obj.NumComponents, columns (obj.mu));
          for i = 1:obj.NumComponents
              fprintf("Clust %d: weight %d\n\tMean: ", i, obj.ComponentProportion(i));
              fprintf("%g ", obj.mu(i,:));
              fprintf("\n");
              if (!obj.SharedCovariance)
                  fprintf("\tVariance:");
                  if (!obj.DiagonalCovariance)
                      if columns (obj.mu) > 1
                        fprintf("\n");
                      endif
                      disp(squeeze(obj.Sigma(:,:,i)))
                  else
                      fprintf(" diag(");
                      fprintf("%g ", obj.Sigma(:,:,i));
                      fprintf(")\n");
                  endif
              endif
          endfor
          if (obj.SharedCovariance)
              fprintf("Shared variance\n");
              if (!obj.DiagonalCovariance)
                  obj.Sigma
              else
                  fprintf(" diag(");
                  fprintf("%g ", obj.Sigma);
                  fprintf(")\n");
              endif
          endif
          if (!isempty (obj.AIC))
              fprintf("AIC=%g BIC=%g NLogL=%g Iter=%d Cged=%d Reg=%g\n", ...
                  obj.AIC, obj.BIC, obj.NegativeLogLikelihood, ...
                  obj.NumIterations, obj.Converged, obj.RegularizationValue);
          endif
      endfunction

        ########################################
        ## Display Gaussian mixture distribution object
      function c = display (obj)
          disp(obj);
      endfunction

        ########################################
        ## Mahalanobis distance to component means
      function D = mahal (obj,X)
        X = checkX (obj, X, "mahal");
        [~, D] = componentProb (obj,X);
      endfunction

        ########################################
        ## Probability density function for Gaussian mixture distribution
      function c = pdf (obj,X)
        X = checkX (obj, X, "pdf");
        p_x_l = componentProb (obj, X);
        c = sum (p_x_l, 2);
      endfunction

        ########################################
        ## Posterior probabilities of components
      function c = posterior (obj,X)
        X = checkX (obj, X, "posterior");
        p_x_l = componentProb (obj, X);
        c = bsxfun(@rdivide, p_x_l, sum (p_x_l, 2));
      endfunction

        ########################################
        ## Random numbers from Gaussian mixture distribution
      function c = random (obj,n)
          c = zeros (n, obj.NumVariables);
          classes = randsample (obj.NumVariables, n, true, obj.ComponentProportion);
          if (obj.SharedCovariance)
              if (obj.DiagonalCovariance)
                  sig = diag (obj.Sigma);
              else
                  sig = obj.Sigma;
              endif
          endif
          for i = 1:obj.NumComponents
              idx = (classes == i);
              k = sum(idx);
              if k > 0
                if (!obj.SharedCovariance)
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
    methods (Static)
      ## Gaussian mixture parameter estimates
      function c = fit  (X,k,varargin)
          c = fitgmdist (X,k,varargin);
      endfunction
    endmethods

    ########################################
    methods (Access = private)
      ## Probability density of (row of) X *and* component l
      ## Second argument is an array of the Mahalonis distances
      function [p_x_l, M] = componentProb (obj, X)
        M     = zeros (rows (X), obj.NumComponents);
        dets  = zeros (1, obj.NumComponents);   % sqrt(determinant)
        if (obj.SharedCovariance)
          if (obj.DiagonalCovariance)
            r = diag (sqrt(obj.Sigma));
          else
            r = chol (obj.Sigma);
          endif
        endif
        for i = 1:obj.NumComponents
          dev = bsxfun (@minus, X, obj.mu(i,:));
          if (!obj.SharedCovariance)
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
        coeff = obj.ComponentProportion ./ ((2*pi)^(obj.NumVariables/2).*dets);
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
          end
        endif
      endfunction
    endmethods
endclassdef

%!test
%! mu = eye(2);
%! Sigma = eye(2);
%! GM = gmdistribution (mu, Sigma);
%! density = GM.pdf ([0 0; 1 1]);
%! assert (density(1) - density(2), 0, 1e-6);
%!
%! [idx, nlogl, P, logpdf,M] = cluster (GM, eye(2));
%! assert (idx, [1; 2]);
%! [idx2,nlogl2,P2,logpdf2] = GM.cluster (eye(2));
%! assert (nlogl - nlogl2, 0, 1e-6);
%! [idx3,nlogl3,P3] = cluster (GM, eye(2));
%! assert (P - P3, zeros (2), 1e-6);
%! [idx4,nlogl4] = cluster (GM, eye(2));
%! assert (size (nlogl4), [1 1]);
%! idx5 = cluster (GM, eye(2));
%! assert (idx - idx5, zeros (2,1));
%!
%! D = GM.mahal ([1;0]);
%! assert (D - M(1,:), zeros (1,2), 1e-6);
%!
%! P = GM.posterior ([0 1]);
%! assert (P - P2(2,:), zeros (1,2), 1e-6);
%!
%! R = GM.random(20);
%! assert (size(R), [20, 2]);

