function [result, model, pupil, initial, elapsed] = prt_gui_trace(rows, options)
%PRT_GUI_TRACE Validate the surface table and trace its circular pupil.
isLde=~isempty(rows) && strcmp(rows{1,1},'OBJ');
objectAtInfinity=false;
if isLde
    coatings=cell(1,size(rows,1));
    if isfield(options,'coatings'), coatings=options.coatings; end
    [surfaces,startZ,stopIndex,objectAtInfinity]=prt_lde_model(rows,coatings);
else
    surfaces=parseSurfaces(rows);
    startZ=options.startZ;
end
overrides=[];
if isfield(options,'globalPoses'), overrides=options.globalPoses; end
[surfaces,source,poses]=prt_scene_model(surfaces,startZ,overrides);
source.atInfinity=objectAtInfinity;
entrance=source; exit=[]; pupilStatus='';
if isLde && isfield(options,'pupilConfig')
    [entrance,exit,pupilStatus]=prt_pupil_frames(surfaces,source,stopIndex,options.pupilConfig);
end
model=struct('surfaces',surfaces, ...
    'pupilDiameter_mm',options.diameter,'startZ_mm',source.origin(3), ...
    'samplesPerAxis',options.samples,'wavelength_nm',options.wavelength, ...
    'polarization',options.polarizationName,'includeOPLPhase',options.phase, ...
    'sourceFrame',source,'globalPoses',poses,'entrancePupilFrame',entrance, ...
    'exitPupilFrame',exit,'pupilStatus',pupilStatus);
if isLde
    model.ldeRows=rows;
    model.ldeCoatings_nm=coatings;
    model.stopSurface=stopIndex;
end
mapped=isfield(options,'pupilData') && ~isempty(options.pupilData);
if mapped
    pupil=options.pupilData.xy;
else
    radius=options.diameter/2;
    [x,y]=ndgrid(linspace(-radius,radius,options.samples));
    mask=hypot(x,y)<=radius;
    pupil=[reshape(x(mask),1,[]);reshape(y(mask),1,[])];
end
count=size(pupil,2);
entrancePoints=entrance.origin+entrance.rotation*[pupil;zeros(1,count)];
launchOffset=dot(source.origin-entrance.origin,source.rotation(:,3));
initial=entrancePoints+launchOffset*source.rotation(:,3);
direction=repmat(source.rotation(:,3),1,count);
% Presets 4/5: right/left circular, with exp(+i*omega*t), looking toward the source.
states=[1 0 1/sqrt(2) 1/sqrt(2) 1/sqrt(2); ...
        0 1 1/sqrt(2) 1i/sqrt(2) -1i/sqrt(2);zeros(1,5)];
if mapped
    inputModes=options.pupilData.modes;
    model.inputMap=rmfield(options.pupilData,'modes');
    model.pupilDiameter_mm=max(options.diameter,2*max(vecnorm(pupil)));
elseif options.polarization==6
    modes=prt_polarization(options.customForm,options.customValues);
    model.customInput=struct('form',options.customForm,'values',options.customValues);
elseif options.polarization==7
    modes=[-1;1]/sqrt(2);
else
    modes=states(1:2,options.polarization);
end
if ~mapped
    model.inputModesXY=modes;
    model.coherencyXY=modes*modes';
    inputModes=repmat(reshape(modes,2,1,[]),1,count,1);
end
modeCount=size(inputModes,3);
hasJonesField=modeCount==1 && (~mapped || strcmp(options.pupilData.format,'Jones'));
if ~isLde && isfield(options,'coatings')
    assert(numel(options.coatings)==numel(model.surfaces),'膜系数量与光面数不一致。');
    for h=1:numel(model.surfaces)
        layers=options.coatings{h};
        if ~isempty(layers), layers(:,3)=layers(:,3)*1e-6; end
        model.surfaces(h).layers=layers;
    end
end
wavelength=options.wavelength*1e-6;
timer=tic;
exitResult=[];
imaR=eye(3); imaOrigin=[0;0;surfaces(end).z];
if isfield(surfaces,'frame') && ~isempty(surfaces(end).frame)
    imaR=surfaces(end).frame.rotation; imaOrigin=surfaces(end).frame.origin;
end
for j=1:modeCount
    field=source.rotation*[inputModes(:,:,j);zeros(1,count)];
    if isempty(exit)
        ray=prt_trace(model.surfaces,initial,direction,field,wavelength,options.phase);
    else
        [ray,pre]=prt_trace(model.surfaces,initial,direction,field,wavelength,options.phase);
        if isempty(pre)
            pre=ray;
        end
        pre.OPL=pre.OPL+surfaces(1).n1*launchOffset;
        pupilRay=prt_exit_pupil(pre,exit,surfaces(end).n1,wavelength,options.phase);
        if j==1
            exitResult=pupilRay;
            exitModes=complex(NaN(3,count,modeCount));
        end
        exitModes(:,:,j)=pupilRay.E;
    end
    if j==1
        result=ray;
        propagated=complex(NaN(3,count,modeCount));
    end
    propagated(:,:,j)=ray.E;
end
result.OPL=result.OPL+surfaces(1).n1*launchOffset;
result.modeFields=propagated;
result.totalFieldSquared=sum(sum(abs(propagated).^2,1),3);
result.stokesXY=fieldStokes(propagated);
result.input=struct('xy',pupil,'stokesXY',fieldStokes(inputModes), ...
    'valid',true(1,count),'OPL',zeros(1,count), ...
    'totalFieldSquared',sum(sum(abs(inputModes).^2,1),3),'modeFields',inputModes);
if isLde
    local=imaR.'*(result.position-imaOrigin);
    localModes=reshape(imaR.'*reshape(propagated,3,[]),size(propagated));
    result.image=struct('xy',local(1:2,:),'stokesXY',fieldStokes(localModes), ...
        'valid',result.valid,'OPL',result.OPL,'totalFieldSquared',result.totalFieldSquared, ...
        'jonesXY',[]);
    if hasJonesField, result.image.jonesXY=localModes(1:2,:,1); end
end
if ~isempty(exitResult)
    exitResult.modeFields=exitModes;
    exitResult.totalFieldSquared=sum(sum(abs(exitModes).^2,1),3);
    exitResult.stokesXY=fieldStokes(exitModes);
    exitResult.jonesXY=[];
    if hasJonesField, exitResult.jonesXY=exitResult.E(1:2,:); end
    if ~hasJonesField, exitResult.E(:)=complex(NaN,NaN); end
end
result.exitPupil=exitResult;
if ~hasJonesField
    result.E(:)=complex(NaN,NaN); % Mixed/phase-undefined input has no coherent Jones map.
end
elapsed=toc(timer);
end

function S=fieldStokes(modes)
Ex=modes(1,:,:); Ey=modes(2,:,:);
px=sum(abs(Ex).^2,3); py=sum(abs(Ey).^2,3);
xy=sum(conj(Ex).*Ey,3);
S=[px+py;px-py;2*real(xy);2*imag(xy)];
end

function surfaces=parseSurfaces(rows)
assert(~isempty(rows),'请至少设置一个光面。');
surfaces=repmat(struct('radius',Inf,'conic',0,'A',zeros(1,4),'z',0, ...
    'aperture',Inf,'n1',1,'n2',1,'mode','transmit'),1,size(rows,1));
expectedN=[];
for h=1:size(rows,1)
    v=rows(h,2:12);
    assert(all(cellfun(@(x) isnumeric(x) && isscalar(x) && isreal(x) && ~isnan(x),v)), ...
        '第 %d 面存在空白或无效数值。',h);
    v=cell2mat(v);
    assert(all(isfinite(v([1,3,5:11]))) && v(2)~=0 && all(v(4:6)>0) && v(7)>=0, ...
        '第 %d 面参数无效。R 不能为零，半口径和实折射率须为正，κ 须非负。',h);
    assert(isempty(expectedN) || abs(v(5)-expectedN)<=1e-10*max(1,expectedN), ...
        '第 %d 面的 n1 与上一面的出射介质折射率不一致。',h);
    reflect=strcmp(rows{h,1},'反射');
    assert(reflect || v(7)==0,'第 %d 面透射模式要求 κ=0。',h);
    mode='transmit';
    expectedN=v(6);
    if reflect
        mode='reflect';
        expectedN=v(5);
    end
    surfaces(h)=struct('radius',v(2),'conic',v(3),'A',v(8:11),'z',v(1), ...
        'aperture',v(4),'n1',v(5),'n2',complex(v(6),-v(7)),'mode',mode);
end
end
