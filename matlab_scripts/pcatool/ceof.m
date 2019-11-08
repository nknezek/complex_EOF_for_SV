% Complex (Hilbert) EOF:
%
%  [lamda, loadings, pcs, per] = ceof(data, nkp);
%
%   Outputs:
%  lamda    = eigenvalues (should be pretty close to real)
%  loadings = First 10 Complex Loadings (row dimension)
%  pcs      = First 10 Complex Principal Components (column dim)
%  per      = percent variance explained (real)
%
%  Inputs:
%  data     = data, so that size(data) = [ntime nspace]
%  nkp      = number of modes to output (default = 10);
%
%  Note:  pcs can be found by performing the following:
%    pcs = data * loadings(:,1:10);
%
%  Normalization is such that:
%    loadings' * loadings = diag(ones(1,npt));
%    pcs' * pcs = diag(lamda);
%
%  For display purposes, the following patterns and time
%    series go together:
%         real(loadings) goes with real(pcs);
%         imag(loadings) goes with imag(pcs) = hilbert(real(pcs))
%
%  Also, one can divide the pcs by sqrt(lamda) and multiply the
%    loadings by sqrt(lamda) to get actual amplitudes.  Recall,
%    std(real(pcs)) should equal std(imag(pcs)).



 function [lamda, loadings, pcs, per] = ceof(data, nkp);

  a = 0.1:.1:10;
  for i = 0:199
    data((i+1),:) = sin(pi*(0.5*a + 0.1*i)) + 5*(rand(1,100)-0.5);
  end
  data = (data - ones(200,1)*mean(data));
  
  if nargin < 2; nkp = 10; end;
  if nargin < 1; error('Need to input data'); end;

  [ntim, npt] = size(data);

  disp('Calculating hilbert matrix')
  %data = data + j * hilbert(data);
  data = hilbert(data);
  disp('Done with hilbert matrix, calculating covariance matrix')
  c = data' * data / ntim;

  disp('Covariance matrix computed, starting eig(c)')

  [loadings, lamda] = eig(c);
  l = diag(lamda);

  [lamda,k] = sort(l'); loadings = loadings(:,k);
  lamda      = fliplr(lamda);
  loadings   = fliplr(loadings);
  loadings = loadings(:,1:nkp);
  per = real(lamda*100/sum(lamda));
  pcs = data * loadings;

%  loadings = loadings(:,1:100);
%  pcs = data * loadings(:,1:10);


