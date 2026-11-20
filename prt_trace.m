function [out,beforeLast] = prt_trace(surfaces,P0,k,E,wavelength,includeOPLPhase,initialOPL)
% Numeric 3-by-N batch tracer; surface n1/n2 are explicit, including mirrors.
% wavelength OPTIONAL, same unit as positions; e.g. 1550 nm = 0.00155 mm.
% Optional surfaces(h).layers uses [n,kappa,thickness] in wavelength units.
% includeOPLPhase defaults to true when wavelength is supplied.
if nargin<5, wavelength=[]; end
if nargin<6, includeOPLPhase=~isempty(wavelength); end
assert(size(P0,1)==3 && isequal(size(P0),size(k),size(E)));
assert(isreal(P0) && isreal(k) && all(isfinite(P0(:))) && all(isfinite(k(:))) && all(isfinite(E(:))));
assert(all(vecnorm(k)>0)); k=k./vecnorm(k);
E=E-k.*sum(k.*E,1); % Explicit transverse projection of supplied global input.
count=size(k,2); valid=true(1,count); opl=zeros(1,count);
if nargin>=7
    validateattributes(initialOPL,{'numeric'},{'real','finite','size',[1,count]});
    opl=initialOPL;
end
history=cell(1,numel(surfaces));
beforeLast=[];
for h=1:numel(surfaces)
    if nargout>1 && h==numel(surfaces)
        beforeLast=struct('E',E,'position',P0,'direction',k,'valid',valid,'OPL',opl);
    end
    ids=find(valid); if isempty(ids), break; end
    s=surfaces(h);
    [P,N,t,hit,iterations]=prt_intersect(P0(:,ids),k(:,ids),s);
    valid(ids(~hit))=false; good=ids(hit);
    coating=[];
    if ~isempty(good)
        layers=[];
        if isfield(s,'layers'), layers=s.layers; end
        [E(:,good),k(:,good),coeff]=prt_interface(E(:,good),k(:,good),N(:,hit),s.n1,s.n2,s.mode,layers,wavelength);
        if ~isempty(layers), coating=coeff; end
        P0(:,good)=P(:,hit); opl(good)=opl(good)+s.n1*t(hit);
    end
    history{h}=struct('rayIndices',ids,'position',P,'valid',hit,'iterations',iterations,'coating',coating);
end
E(:,~valid)=NaN; P0(:,~valid)=NaN; k(:,~valid)=NaN; opl(~valid)=NaN;
% Geometric OPL phase only. Common phase removed to reduce argument size.
if includeOPLPhase
    validateattributes(wavelength,{'numeric'},{'scalar','real','positive','finite'});
    if any(valid)
        reference=opl(find(valid,1));
        E(:,valid)=E(:,valid).*exp(-2i*pi*(opl(valid)-reference)/wavelength);
    end
end
out=struct('E',E,'position',P0,'direction',k,'valid',valid,'OPL',opl,'history',{history});
end
