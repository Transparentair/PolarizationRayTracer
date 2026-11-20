function [etaLocal,etaCoherent,valid] = prt_overlap(Es,El,weights)
% Fields must already share a physical plane, Cartesian frame and quadrature.
assert(isequal(size(Es),size(El)) && size(Es,1)==3);
if nargin<3, weights=ones(1,size(Es,2)); end
weights=reshape(weights,1,[]);
assert(numel(weights)==size(Es,2) && isreal(weights) && all(isfinite(weights)) && all(weights>=0));
ps=sum(abs(Es).^2,1); pl=sum(abs(El).^2,1);
valid=all(isfinite(Es),1) & all(isfinite(El),1) & ps>0 & pl>0 & weights>0;
inner=sum(conj(El).*Es,1); etaLocal=NaN(size(ps));
etaLocal(valid)=min(1,abs(inner(valid)).^2./(ps(valid).*pl(valid)));
w=weights(valid); denom=sum(w.*ps(valid))*sum(w.*pl(valid));
if denom>0, etaCoherent=min(1,abs(sum(w.*inner(valid)))^2/denom);
else, etaCoherent=NaN; end
end
