## Copyright (C) 2007 Roman Stanchak
##
## This file is part of Octave.
##
## Octave is free software; you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2, or (at your option)
## any later version.
##
## Octave is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with Octave; see the file COPYING.  If not, write to the Free
## Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
## 02110-1301, USA.

## -*- texinfo -*-
## @deftypefn {Function File} hist3(@var{X})
## @deftypefnx {Function File} hist3(@var{X}, @var{nbins})
## @deftypefnx {Function File} hist3(@var{X}, 'Nbins', @var{nbins})
## @deftypefnx {Function File} hist3(@var{X}, @var{centers})
## @deftypefnx {Function File} hist3(@var{X}, 'Centers', @var{centers})
## @deftypefnx {Function File} hist3(@var{X}, 'Edges', @var{edges})
## @deftypefnx {Function File} {@var{N} =} hist3(@var{X}, ...)
## @deftypefnx {Function File} {[@var{N}, @var{C}] =} hist3(@var{X}, ...)
## Plots a 2D histogram of the N x 2 matrix @var{X} with 10 equally spaced
## bins in both the x and y direction using the @code{mesh} function
##
## The number of equally spaced bins to compute histogram can be specified with
## @var{nbins}. If @var{nbins} is a 2 element vector, use the two values as the
## number of bins in the x and y axis, respectively, otherwise, use the same
## value for each.
##
## The centers of the histogram bins can be specified with @var{centers}.
## @var{centers} should be a cell array containing two arrays of the bin
## centers on the x and y axis, respectively.
##
## The edges of the histogram bins can be specified with @var{edges}.
## @var{edges} should be a cell array containing two arrays of the bin edges
## on the x and y axis, respectively.
##
## @var{N} returns the 2D array of bin counts, and does not plot the
## histogram
##
## @var{N} and @var{C} returns the 2D array of bin counts in @var{N} and the
## bin centers in the 2 element cell array @var{C}, and does not plot the
## histogram
##
## @seealso{hist, mesh}
## @end deftypefn

## Authors: Paul Kienzle (segments borrowed from hist2d),
##          Roman Stanchak (addition of matlab compatible syntax, bin edge arg)

function varargout = hist3(varargin)
  methods={'Nbins', 'Centers', 'Edges'};
  method=1;
  xbins=10;
  ybins=10;
  edges={};
  M=varargin{1};
  if nargin>=2,
    % is a binning method is specified?
    if isstr(varargin{2}),
      method = find(strcmp(methods,varargin{2}));
      if isempty(method),
        error('Unknown property string');
      elseif nargin < 3,
        error('Expected an additional argument');
      elseif method==2,
        xbins = varargin{3}{1};
        ybins = varargin{3}{2};
      elseif method==3
        edges = varargin{3};
      end
    elseif iscell(varargin{2}) % second argument contains centers
      method = 2;
      xbins = varargin{2}{1};
      ybins = varargin{2}{2};
    elseif isscalar(varargin{2}),
      xbins = ybins = varargin{2};
    elseif isvector(varargin{2}), % second argument contain number of bins
      xbins = varargin{2}(1);
      ybins = varargin{2}(2);
    else
      error('Unsupported type for 2nd argument');
    end
  end

  % If n bins, find centers based on n+1 bin edges
  if method==1,
    lo = min(M);
    hi = max(M);
    if isscalar(xbins)
      xbins = linspace(lo(1),hi(1),xbins+1);
        xbins = (xbins(1:end-1)+xbins(2:end))/2;
    end
    if isscalar(ybins)
      ybins = linspace(lo(2),hi(2),ybins+1);
      ybins = (ybins(1:end-1)+ybins(2:end))/2;
    end
    method=2;
  end

  if method==2,
    % centers specified, compute edges
    xcut = (xbins(1:end-1)+xbins(2:end))/2;
    ycut = (ybins(1:end-1)+ybins(2:end))/2;
    xidx = lookup(xcut,M(:,1))+1;
    yidx = lookup(ycut,M(:,2))+1;
  else,
    % edges specified.  Filter points outside edge range
    xidx = lookup(edges{1},M(:,1));
    yidx = lookup(edges{2},M(:,2));
    idx = find(xidx>0);
    idx = intersect(idx, find(xidx<length(edges{1})));
    idx = intersect(idx, find(yidx>0));
    idx = intersect(idx, find(yidx<length(edges{2})));
    xidx=xidx(idx);
    yidx=yidx(idx);

    % compute centers for plotting
    xbins = (edges{1}(1:end-1)+edges{1}(2:end))/2;
    ybins = (edges{2}(1:end-1)+edges{2}(2:end))/2;
  end
  counts = sparse(yidx,xidx,1,length(ybins), length(xbins), 'sum');

  if nargout
    varargout{1} = full(counts);
    if nargout>1
      varargout{2} = {xbins,ybins};
    end
  else
    mesh(xbins,ybins,full(counts));
  end
end
