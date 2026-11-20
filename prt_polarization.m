function modes = prt_polarization(form, v)
% Columns are mutually incoherent Jones modes in XY, including amplitude.
% Ellipse uses [azimuth_deg, ellipticity_deg, S0].
% Stokes uses S3=2*imag(conj(Ex)*Ey), matching prt_ellipse.
v=v(:);
validateattributes(v,{'numeric'},{'real','finite','nonempty'});
switch form
    case 'jones'
        assert(numel(v)==4,'Jones 需要 Ex、Ey 的实部与虚部。');
        modes=[complex(v(1),v(2));complex(v(3),v(4))];
    case 'amplitude'
        assert(numel(v)==3 && all(v(1:2)>=0),'振幅必须非负。');
        modes=[v(1);v(2)*exp(1i*deg2rad(v(3)))];
    case 'ellipse'
        assert(numel(v)==3 && abs(v(2))<=45 && v(3)>0, ...
            '椭圆率角须在 ±45° 内，S0 须大于零。');
        a=deg2rad(v(1)); b=deg2rad(v(2));
        modes=sqrt(v(3))*[cos(a)*cos(b)-1i*sin(a)*sin(b); ...
                         sin(a)*cos(b)+1i*cos(a)*sin(b)];
    case 'stokes'
        assert(numel(v)==4 && v(1)>0,'Stokes 需要四个分量，且 S0 须大于零。');
        p=v(2:4)/v(1);
        assert(norm(p)<=1+1e-12,'须满足 S1² + S2² + S3² ≤ S0²。');
        C=0.5*[1+p(1),complex(p(2),-p(3));complex(p(2),p(3)),1-p(1)];
        [U,d]=eig(C,'vector');
        keep=d>1e-14;
        modes=U(:,keep).*reshape(sqrt(d(keep))*sqrt(v(1)),1,[]);
    otherwise
        error('prt:polarizationForm','未知偏振描述形式。');
end
assert(all(isfinite(modes(:))) && any(abs(modes(:))>0),'输入电场不能为零或溢出。');
end
