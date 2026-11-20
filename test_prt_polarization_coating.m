function test_prt_polarization_coating
% Analytic checks and an independent characteristic-matrix comparison.
ec=[1;1i]/sqrt(2);
for state={{'jones',[real(ec(1)),imag(ec(1)),real(ec(2)),imag(ec(2))]}, ...
        {'amplitude',[1/sqrt(2),1/sqrt(2),90]}, ...
        {'ellipse',[23,45,1]}, {'stokes',[1,0,0,1]}}
    data=state{1};
    modes=prt_polarization(data{1},data{2});
    assert(norm(modes*modes'-ec*ec','fro')<1e-12);
end
modes=prt_polarization('stokes',[2,0.4,0.6,-0.8]);
assert(norm(modes*modes'-[1.2,0.3+0.4i;0.3-0.4i,0.8],'fro')<1e-12);
assert(size(prt_polarization('stokes',[1 0 0 0]),2)==2);
wave=550;
ci=cosd(0:2:78);
k=[sqrt(1-ci.^2);zeros(size(ci));ci];
E=repmat([0;1;0],1,numel(ci));
[~,~,bare]=prt_interface(E,k,[0;0;-1],1,1.52,'transmit');
[rs,rp,ts,tp]=prt_coating(1,1.52,ci,[1.38 0 0],wave);
assert(max(abs([rs-bare.rs,rp-bare.rp,ts-bare.ts,tp-bare.tp]))<1e-12);
nf=sqrt(1.52);
[rs,rp,ts,tp]=prt_coating(1,1.52,1,[nf 0 wave/(4*nf)],wave);
assert(max(abs([rs rp]))<1e-12 && max(abs(1.52*abs([ts tp]).^2-1))<1e-12);
films=[2.15 0 92;1.38 0 105;1.82 0 74];
[rs,rp,ts,tp]=prt_coating(1,1.52,ci,films,wave);
flux=sqrt(1.52^2-(1-ci.^2))./ci;
assert(max(abs([abs(rs).^2+flux.*abs(ts).^2-1,abs(rp).^2+flux.*abs(tp).^2-1]))<1e-12);
for j=1:numel(ci)
    n=[1;films(:,1);1.52];
    q=sqrt(n.^2-(1-ci(j)^2));
    for p=0:1
        adm=q;
        if p, adm=n.^2./q; end
        M=eye(2);
        for m=1:size(films,1)
            d=2*pi*q(m+1)*films(m,3)/wave;
            Y=adm(m+1);
            M=M*[cos(d),1i*sin(d)/Y;1i*Y*sin(d),cos(d)];
        end
        v=M*[1;adm(end)];
        rTang=(adm(1)*v(1)-v(2))/(adm(1)*v(1)+v(2));
        tTang=2*adm(1)/(adm(1)*v(1)+v(2));
        if p
            assert(abs(rp(j)+rTang)<1e-12);
            assert(abs(tp(j)-tTang*ci(j)/(q(end)/n(end)))<1e-12);
        else
            assert(abs(rs(j)-rTang)<1e-12 && abs(ts(j)-tTang)<1e-12);
        end
    end
end
[rs,rp,ts,tp]=prt_coating(1,1.52,ci,[1.7 0.25 350],wave);
assert(all([abs(rs).^2+flux.*abs(ts).^2,abs(rp).^2+flux.*abs(tp).^2]<1));
[rs,rp,ts,tp]=prt_coating(1,1.52,ci,[1.7 10 1e8],wave);
assert(all(isfinite([rs rp ts tp])) && all(ts==0) && all(tp==0));
[rs,rp]=prt_coating(1.5,1,cosd(65),[1.2 0 120],wave);
assert(max(abs(abs([rs rp])-1))<1e-12);
% Use the same path as the GUI: unpolarized input, oblique incidence across a sphere.
row={'透射',0,10,0,8,1,1.52,0,0,0,0,0,1};
opts=struct('diameter',12,'startZ',-10,'samples',21,'wavelength',wave, ...
    'polarization',6,'polarizationName','Stokes','phase',false, ...
    'customForm','stokes','customValues',[1 0 0 0],'coatings',{{[nf 0 wave/(4*nf)]}});
[r,model,pupil]=prt_gui_trace(row,opts);
center=find(all(pupil==0,1),1);
assert(all(r.valid) && size(r.modeFields,3)==2 && all(isnan(r.E(:))));
assert(norm(r.stokesXY(:,center)-[1/1.52;0;0;0])<1e-12);
assert(max(abs(r.history{1}.coating.rp))>1e-4 && ...
    abs(model.surfaces.layers(3)-wave/(4*nf)*1e-6)<1e-15);
fprintf('PASS: polarization forms, mixed states, film coefficients and GUI trace integration\n');
end
