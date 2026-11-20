function count = prt_save_pupil(file,data,format)
% Export valid ray samples in the selected plane's local XY coordinates.
[xy,jones,stokes,vector]=prt_pupil_csv_columns();
valid=reshape(data.valid,1,[]);
positions=data.xy(:,valid).';
count=size(positions,1);
assert(count>0,'当前平面没有有效追迹点可导出。');
assert(isreal(positions) && all(isfinite(positions),'all'),'有效追迹点含非有限坐标。');
assert(size(unique(positions,'rows'),1)==count, ...
    '当前平面有重合的光线坐标，无法写成逐点偏振图。请改选出瞳，或用诊断 CSV 保留所有光线。');
switch lower(format)
    case 'jones'
        assert(isfield(data,'jonesXY') && ~isempty(data.jonesXY), ...
            '当前结果没有可导出的单一 Jones 复电场。请选择 Stokes。');
        E=data.jonesXY(:,valid).';
        values=[real(E(:,1)),imag(E(:,1)),real(E(:,2)),imag(E(:,2))];
        columns=[xy,jones];
    case 'stokes'
        values=data.stokesXY(:,valid).';
        columns=[xy,stokes];
    case 'jones3d'
        assert(isfield(data,'jonesGlobal') && ~isempty(data.jonesGlobal), ...
            '当前结果没有可导出的单一 Jones 三维复电场。请选择 Stokes。');
        E=data.jonesGlobal(:,valid);
        k=prt_require_transverse(E,data.directionGlobal(:,valid));
        values=[real(E(1,:));imag(E(1,:));real(E(2,:));imag(E(2,:)); ...
                real(E(3,:));imag(E(3,:));k].';
        columns=[xy,vector];
    otherwise
        error('prt:csvFormat','请选择 Jones 或 Stokes。');
end
assert(isreal(values) && all(isfinite(values),'all'),'有效追迹点含非有限偏振数据。');
[fid,message]=fopen(file,'w');
assert(fid>=0,'无法写入文件：%s',message);
cleanup=onCleanup(@() fclose(fid));
fprintf(fid,'%s\n',strjoin(columns,','));
rowFormat=[repmat('%.17g,',1,numel(columns)-1),'%.17g\n'];
fprintf(fid,rowFormat,[positions,values].');
[message,code]=ferror(fid);
assert(code==0,'写入 CSV 失败：%s',message);
clear cleanup
end
