function s=asym(X)


error(nargchk(1,2,nargin))
sz = size(X);
if nargin<2||isempty(dim),
  % Use 1'st non-singleton dimension or dimension 1
  dim = find(sz~=1, 1 );
  if isempty(dim), dim = 1; end
end


hil=imag(hilbert(X-repmat(mean(X),size(X,1),1)));
% hil=sqrt(hil.*conj(hil));

rsz = ones(size(sz)); rsz(dim)=sz(dim);
mu  = mean(X,dim);
if isscalar(mu)
    Xmu  = (X-mu);
else
    Xmu  = (X-repmat(mu,rsz));
end

s   = mean(hil.^3,dim)./mean(Xmu.^2,dim).^(1.5);
end