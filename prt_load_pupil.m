function data = prt_load_pupil(file)
% Each CSV row specifies one launch sample in entrance-pupil local mm.
T=readtable(file,'VariableNamingRule','preserve');
names=T.Properties.VariableNames;
[xy,jones,stokes]=prt_pupil_csv_columns();
assert(height(T)>0 && all(ismember(xy,names)),'CSV 需要 X_mm、Y_mm 和至少一行数据。');
assert(xor(all(ismember(jones,names)),all(ismember(stokes,names))), ...
    '请选择一套完整列：ExRe,ExIm,EyRe,EyIm 或 S0,S1,S2,S3。');
columns=jones;
form='Jones';
if all(ismember(stokes,names)), columns=stokes; form='Stokes'; end
assert(all(varfun(@isnumeric,T(:,[xy columns]),'OutputFormat','uniform')), ...
    '坐标和偏振分量必须为数值。');
a=table2array(T(:,[xy columns])).';
assert(isreal(a) && all(isfinite(a),'all'),'CSV 坐标和偏振分量须为有限实数。');
assert(size(unique(a(1:2,:).','rows'),1)==size(a,2),'CSV 存在重复的入瞳坐标。');
pupil=a(1:2,:);
v=a(3:6,:);
if strcmp(form,'Jones')
    modes=[complex(v(1,:),v(2,:));complex(v(3,:),v(4,:))];
else
    S0=v(1,:); P=vecnorm(v(2:4,:));
    assert(all(S0>=0 & P<=S0.*(1+1e-12)), ...
        'Stokes 每行须满足 S0 ≥ 0 且 sqrt(S1²+S2²+S3²) ≤ S0。');
    unit=zeros(3,numel(P));
    polarized=P>0;
    unit(:,polarized)=v(2:4,polarized)./P(polarized);
    u=complex(zeros(2,numel(P)));
    aDominant=unit(1,:)>=0;
    u(1,aDominant)=sqrt(0.5*(1+unit(1,aDominant)));
    u(2,aDominant)=complex(unit(2,aDominant),unit(3,aDominant))./(2*u(1,aDominant));
    u(2,~aDominant)=sqrt(0.5*(1-unit(1,~aDominant)));
    u(1,~aDominant)=complex(unit(2,~aDominant),-unit(3,~aDominant))./(2*u(2,~aDominant));
    u(1,~polarized)=1;
    P=min(P,S0);
    modes=u.*sqrt(0.5*S0+0.5*P);
    if any(S0-P>1e-14*S0)
        modes(:,:,2)=[-conj(u(2,:));conj(u(1,:))].*sqrt(0.5*(S0-P));
    end
end
data=struct('xy',pupil,'modes',modes,'format',form,'file',file);
end
