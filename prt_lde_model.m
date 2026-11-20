function [surfaces,startZ,stopIndex,objectAtInfinity] = prt_lde_model(rows,coatings)
% Compile a coaxial sequential prescription into the ray tracer's surfaces.
% OBJ is a finite launch plane for collimated rays. Surface 1 is at Z=0.
assert(iscell(rows) && size(rows,2)==13 && size(rows,1)>=2, ...
    'LDE 需要 OBJ、IMA 行以及 13 列参数。');
count=size(rows,1);
assert(strcmp(rows{1,1},'OBJ') && strcmp(rows{end,1},'IMA'), ...
    '第一行必须为 OBJ，最后一行必须为 IMA。');
if nargin<2, coatings=cell(1,count); end
assert(iscell(coatings) && numel(coatings)==count,'膜系数量与 LDE 行数不一致。');
assert(isempty(coatings{1}) && isempty(coatings{end}),'OBJ 和 IMA 不能镀膜。');
distance=rows{1,4};
validateattributes(distance,{'numeric'},{'real','nonnan','scalar','nonnegative'});
objectAtInfinity=isinf(distance);
if objectAtInfinity
    % Use a finite numerical launch plane upstream of the first optical face.
    distance=10;
    if count>2 && isfinite(rows{2,7}), distance=max(distance,2*rows{2,7}); end
end
startZ=-distance;
n=readIndex(rows{1,5},0);
surfaces=repmat(struct('radius',Inf,'conic',0,'A',zeros(1,4),'z',0, ...
    'aperture',Inf,'n1',1,'n2',1,'mode','transmit','layers',[], ...
    'number',0,'isStop',false,'comment',''),1,count-1);
stopIndex=[];
z=0;
for row=2:count
    h=row-1;
    s=surfaces(h);
    s.number=h;
    s.z=z;
    s.n1=n;
    s.n2=n;
    s.comment=char(string(rows{row,2}));
    s.aperture=rows{row,7};
    validateattributes(s.aperture,{'numeric'},{'real','scalar','positive','nonnan'}, ...
        'LDE','Clear Semi-Dia.');
    if row==count
        surfaces(h)=s;
        break
    end
    kind=rows{row,1};
    assert(any(strcmp(kind,{'Standard','Even Asphere'})), ...
        '第 %d 面请选择 Standard 或 Even Asphere。',h);
    s.radius=rows{row,3};
    validateattributes(s.radius,{'numeric'},{'real','scalar','nonnan'},'LDE','Radius');
    if s.radius==0, s.radius=Inf; end
    thickness=rows{row,4};
    validateattributes(thickness,{'numeric'},{'real','finite','scalar'},'LDE','Thickness');
    s.conic=rows{row,8};
    validateattributes(s.conic,{'numeric'},{'real','finite','scalar'},'LDE','Conic');
    if strcmp(kind,'Even Asphere')
        A=rows(row,9:12);
        assert(all(cellfun(@(a) isnumeric(a) && isreal(a) && isscalar(a) && isfinite(a),A)), ...
            '第 %d 面的 A1–A4 必须是有限实数。',h);
        s.A=cell2mat(A);
    end
    flag=rows{row,13};
    assert((islogical(flag) || isnumeric(flag)) && isscalar(flag) && any(flag==[0 1]), ...
        '第 %d 面 STOP 须为逻辑值。',h);
    s.isStop=logical(flag);
    if s.isStop
        assert(isempty(stopIndex) && isfinite(s.aperture),'只能标记一个 STOP，且其半口径须有限。');
        stopIndex=h;
    end
    layer=coatings{row};
    if ~isempty(layer)
        assert(isnumeric(layer) && isreal(layer) && size(layer,2)==3 && all(isfinite(layer(:))) && ...
            all(layer(:,1)>0 & layer(:,2)>=0 & layer(:,3)>0),'第 %d 面膜层参数无效。',h);
        layer(:,3)=layer(:,3)*1e-6;
    end
    s.layers=layer;
    material=strtrim(char(string(rows{row,5})));
    if strcmpi(material,'MIRROR')
        assert(isempty(layer),'第 %d 面有膜层，请指定反射基底：MIRROR:n,k。',h);
        s.mode='ideal-mirror';
    elseif startsWith(material,'MIRROR:','IgnoreCase',true)
        values=str2double(strsplit(extractAfter(string(material),':'),','));
        assert(numel(values)==2 && isreal(values) && all(isfinite(values)) && ...
            values(1)>0 && values(2)>=0,'第 %d 面反射材料格式为 MIRROR:n,k。',h);
        s.mode='reflect';
        s.n2=complex(values(1),-values(2));
    else
        n=readIndex(material,h);
        s.n2=n;
    end
    surfaces(h)=s;
    z=z+thickness;
end
end

function n=readIndex(material,h)
material=strtrim(char(string(material)));
if isempty(material) || strcmpi(material,'AIR')
    n=1;
else
    n=str2double(material);
end
assert(isscalar(n) && isreal(n) && isfinite(n) && n>0, ...
    '第 %d 面 Material 请输入 AIR 或当前波长的实折射率。',h);
end
