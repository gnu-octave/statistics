## Copyright (C) 2013 Fernando Damian Nieuwveldt <fdnieuwveldt@gmail.com>
##
## This program is free software; you can redistribute it and/or
## modify it under the terms of the GNU General Public License
## as published by the Free Software Foundation; either version 3
## of the License, or (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn {Function File} {[@var{COEFF}]} = princomp(@var{X})
## @deftypefnx {Function File} {[@var{COEFF},@var{SCORE}]} = princomp(@var{X})
## @deftypefnx {Function File} {[@var{COEFF},@var{SCORE},@var{latent}]} = princomp(@var{X})
## @deftypefnx {Function File} {[@var{COEFF},@var{SCORE},@var{latent},@var{tsquare}]} = princomp(@var{X})
## @deftypefnx {Function File} {[...]} = princomp(@var{X},'econ')
## @itemize @bullet
## @item
## princomp performs principal component analysis on a NxP data matrix X
## @item
## @var{COEFF} : returns the principal component coefficients
## @item
## @var{SCORE} : returns the principal component scores, the representation of X 
## in the principal component space
## @item
## @var{LATENT} : returns the principal component variances, i.e., the 
## eigenvalues of the covariance matrix X.
## @item
## @var{TSQUARE} : returns Hotelling's T-squared Statistic for each observation in X 
## @item
## [...] = princomp(X,'econ') returns only the elements of latent that are not 
## necessarily zero, and the corresponding columns of COEFF and SCORE, that is, 
## when n <= p, only the first n-1. This can be significantly faster when p is 
## much larger than n. In this case the svd will be applied on the transpose of 
## the data matrix X
##
## @end itemize
##
## @subheading References
##
## @enumerate
## @item
## Jolliffe, I. T., Principal Component Analysis, 2nd Edition, Springer, 2002
## 
## @end enumerate
## @end deftypefn

function [COEFF,SCORE,latent,tsquare] = princomp(X,varargin)

   if (nargin < 1 || nargin > 2)
      print_usage ();
   endif
	
   if (nargin == 2 && ! strcmpi (varargin{:}, "econ"))
      error ("princomp: if a second input argument is present, it must be the string  'econ'");
   endif

   [nobs nvars] = size(X);
	 
   # Center the columns to mean zero
   Xcentered = bsxfun(@minus,X,mean(X));

   # Check if there are more variables then observations
   if nvars <= nobs
		
      [U,S,COEFF] = svd(Xcentered);

   else

      # Calculate the svd on the transpose matrix, much faster
      if (nargin == 2 && strcmpi ( varargin{:} , "econ"))
	     [COEFF,S,V] = svd(Xcentered' , 'econ');
      else
	     [COEFF,S,V] = svd(Xcentered');
      endif	

   endif

   if nargout > 1

      # Get the Scores
      SCORE = Xcentered*COEFF;
	
      # Get the rank of the SCORE matrix	
      r = rank(SCORE); 

      # Only use the first r columns, pad rest with zeros if economy != 'econ'
      SCORE = SCORE(:,1:r) ; 
	 
      if !(nargin == 2 && strcmpi ( varargin{:} , "econ"))
	    SCORE = [SCORE, zeros(nobs , nvars-r)];
      else
	    COEFF   = COEFF(: , 1:r);   
      endif

    endif

    # This is the same as the eigenvalues of the covariance matrix of X
    latent  = (diag(S'*S)/(size(Xcentered,1)-1))(1:r);

    if nargout > 2
      if !(nargin == 2 && strcmpi ( varargin{:} , "econ"))
	  latent= [latent;zeros(nvars-r,1)];
      endif
    endif

    if nargout > 3
 	# Calculate the Hotelling T-Square statistic for the observations
	tsquare = sumsq(zscore(SCORE(:,1:r)),2);
    endif

endfunction

%!shared COEFF,SCORE,latent,tsquare,m,x

%!test
%! x=[1,2,3;2,1,3]';
%! [COEFF,SCORE,latent,tsquare] = princomp(x);
%! m=[sqrt(2),sqrt(2);sqrt(2),-sqrt(2);-2*sqrt(2),0]/2;
%! m(:,1) = m(:,1)*sign(COEFF(1,1));
%! m(:,2) = m(:,2)*sign(COEFF(1,2));

%!assert(COEFF,m(1:2,:),10*eps);
%!assert(SCORE,-m,10*eps);
%!assert(latent,[1.5;.5],10*eps);
%!assert(tsquare,[4;4;4]/3,10*eps);

%!test
%! x=x';
%! [COEFF,SCORE,latent,tsquare] = princomp(x);
%! m=[sqrt(2),sqrt(2),0;-sqrt(2),sqrt(2),0;0,0,2]/2;
%! m(:,1) = m(:,1)*sign(COEFF(1,1));
%! m(:,2) = m(:,2)*sign(COEFF(1,2));
%! m(:,3) = m(:,3)*sign(COEFF(3,3));

%!assert(COEFF,m,10*eps);
%!assert(SCORE(:,1),-m(1:2,1),10*eps);
%!assert(SCORE(:,2:3),zeros(2),10*eps);
%!assert(latent,[1;0;0],10*eps);
%!xtest
%! assert(tsquare,[0.5;0.5],10*eps)
