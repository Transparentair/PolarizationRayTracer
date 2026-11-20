function [rs,rp,ts,tp] = prt_coating(n1,n2,ci,layers,wavelength)
% Coherent film stack, rows [Re(n), kappa, thickness] in wavelength units.
% Layer order is incidence side to substrate. Convention n-i*kappa.
validateattributes(wavelength,{'numeric'},{'scalar','real','positive','finite'});
assert(isnumeric(layers) && isreal(layers) && size(layers,2)==3 && ...
    all(isfinite(layers(:))) && all(layers(:,1)>0) && ...
    all(layers(:,2)>=0) && all(layers(:,3)>=0),'prt:layer','Invalid film parameters.');
layers=layers(layers(:,3)>0,:);
n=[n1; complex(layers(:,1),-layers(:,2)); n2];
u=n1^2*(1-ci.^2);
qb=forwardQ(n(end),u);
rs=zeros(size(ci)); rp=rs; ts=ones(size(ci)); tp=ts;
prop=ones(size(ci));
for j=numel(n)-1:-1:1
    na=n(j); nb=n(j+1);
    qa=forwardQ(na,u);
    if na==nb
        aS=zeros(size(ci)); aP=aS; bS=ones(size(ci)); bP=bS;
    else
        denS=qa+qb;
        denP=nb^2*qa+na^2*qb;
        aS=(qa-qb)./denS;
        aP=(nb^2*qa-na^2*qb)./denP;
        bS=2*qa./denS;
        bP=2*na*nb*qa./denP;
    end
    % Cascade from the substrate. Only decaying propagation factors occur.
    returnS=rs.*prop.^2;
    returnP=rp.*prop.^2;
    ds=1+aS.*returnS;
    dp=1+aP.*returnP;
    rs=(aS+returnS)./ds;
    rp=(aP+returnP)./dp;
    ts=bS.*prop.*ts./ds;
    tp=bP.*prop.*tp./dp;
    if j>1
        prop=exp(-2i*pi*qa*layers(j-1,3)/wavelength);
    end
    qb=qa;
end
assert(all(isfinite([rs rp ts tp])),'prt:filmSingularity', ...
    '膜系遇到精确临界角或奇异点。请检查折射率和入射角。');
end

function q=forwardQ(n,u)
q=sqrt(complex(n^2-u));
flip=imag(q)>0 | (imag(q)==0 & real(q)<0);
q(flip)=-q(flip);
end
