function [surfaces, source, poses] = prt_scene_model(surfaces,startZ,overrides)
% Rows: [useGlobal, X Y Z (mm), Rx Ry Rz (deg), apertureU apertureV (mm)].
% OBJ is row 1. Remaining rows correspond to surfaces, including IMA.
count=numel(surfaces)+1;
poses=[zeros(count,3),[startZ;reshape([surfaces.z],[],1)],zeros(count,5)];
if nargin>=3 && ~isempty(overrides)
    validateattributes(overrides,{'numeric'},{'real','finite','size',[count 9]});
    assert(all(ismember(overrides(:,1),[0 1])),'全局位姿开关应为 0 或 1。');
    use=logical(overrides(:,1));
    poses(use,:)=overrides(use,:);
end
source=makeFrame(poses(1,:));
for h=1:numel(surfaces)
    if poses(h+1,1)
        surfaces(h).frame=makeFrame(poses(h+1,:));
        surfaces(h).apertureCenter=poses(h+1,8:9).';
    end
end
end

function pose=makeFrame(row)
a=row(5:7);
c=cosd(a); s=sind(a);
Rx=[1 0 0;0 c(1) -s(1);0 s(1) c(1)];
Ry=[c(2) 0 s(2);0 1 0;-s(2) 0 c(2)];
Rz=[c(3) -s(3) 0;s(3) c(3) 0;0 0 1];
pose=struct('origin',row(2:4).','rotation',Rz*Ry*Rx);
end
