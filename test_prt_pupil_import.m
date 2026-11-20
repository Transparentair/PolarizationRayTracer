function test_prt_pupil_import
file=[tempname '.csv'];
cleanup=onCleanup(@() delete(file)); %#ok<NASGU>
xy=[0 0;0.31 -0.22;-0.47 0.64;0.25 0.43];
J=[1 0 0 0;0 0 0 1;0 0 0 0;1 0 1 0];
writeMap(J,{'ExRe','ExIm','EyRe','EyIm'});
rows={'OBJ','Launch',0,5,'AIR',0,Inf,0,0,0,0,0,false; ...
      'Standard','Plane',0,7,'AIR',0,5,0,0,0,0,0,true; ...
      'IMA','Image',0,0,'',0,Inf,0,0,0,0,0,false};
opts=struct('diameter',2,'samples',81,'wavelength',1550,'polarization',3, ...
    'polarizationName','Map','phase',false,'pupilConfig',struct('mode','auto'));
opts.pupilData=prt_load_pupil(file);
[r,~,p]=prt_gui_trace(rows,opts);
assert(isequal(p,xy.') && all(r.valid));
assert(max(abs(r.E(1:2,:)-opts.pupilData.modes),[],'all')<1e-12);
assert(max(abs(r.input.stokesXY-r.exitPupil.stokesXY),[],'all')<1e-12);
assert(max(abs(r.input.stokesXY-r.image.stokesXY),[],'all')<1e-12);
[exported,kept]=roundTrip(r.image,'Jones');
assert(kept==size(xy,1) && isequal(exported.xy,xy.'));
assert(max(abs(exported.modes-opts.pupilData.modes),[],'all')<1e-12);
masked=r.exitPupil;
masked.valid(2)=false;
[exported,kept]=roundTrip(masked,'Jones');
assert(kept==3 && isequal(exported.xy,xy([1,3,4],:).'));
[x,~,~,arrows]=prt_pupil_ellipses(r.input.xy,r.input.stokesXY,r.input.valid,p);
assert(size(x,2)==3 && size(arrows,2)==3);
S=[1 0.2 0.3 0.4;2 0 0 0;1 -1 0 0;1 0 0 -1];
writeMap(S,{'S0','S1','S2','S3'});
opts.pupilData=prt_load_pupil(file);
r=prt_gui_trace(rows,opts);
assert(all(isnan(r.E),'all') && size(r.modeFields,3)==2);
assert(max(abs(r.input.stokesXY-S.'),[],'all')<1e-12);
assert(max(abs(r.exitPupil.stokesXY-S.'),[],'all')<1e-12);
[exported,kept]=roundTrip(r.exitPupil,'Stokes');
assert(kept==4 && isequal(exported.xy,xy.'));
opts.pupilData=exported;
check=prt_gui_trace(rows,opts);
assert(max(abs(check.input.stokesXY-S.'),[],'all')<1e-12);
assert(isempty(r.exitPupil.jonesXY));
bad=r.image;
bad.xy(:,2)=bad.xy(:,1);
rejected=false;
try
    prt_save_pupil(file,bad,'Stokes');
catch problem
    rejected=contains(problem.message,'重合');
end
assert(rejected);

% Both input and exported 3D fields retain Ez, direction and relative phase.
K=[0 0 0.6 0;0 0 0 0;1 1 0.8 1];
E=[1 0 0.8 1;1i 1 1i 0;0 0 -0.6 0];
components=[real(E(1,:));imag(E(1,:));real(E(2,:));imag(E(2,:)); ...
    real(E(3,:));imag(E(3,:));K].';
[~,~,~,columns]=prt_pupil_csv_columns();
writeMap(components,columns);
opts.pupilData=prt_load_pupil(file);
opts.pupilConfig=struct('mode','manual','values',[0 0 -2 0 0 0;0 0 2 0 0 0]);
opts.phase=true;
r=prt_gui_trace(rows,opts);
assert(max(abs(r.input.jonesGlobal-E),[],'all')<1e-12);
assert(max(abs(r.OPL-9./K(3,:)))<1e-12);
expected=E.*exp(-2i*pi*(9./K(3,:)-9)/(opts.wavelength*1e-6));
assert(max(abs(r.image.jonesGlobal-expected),[],'all')<2e-10);
assert(all(abs(sum(r.image.jonesGlobal.*r.image.directionGlobal,1))<1e-12));
[exported,kept]=roundTrip(r.image,'Jones3D');
assert(kept==4 && max(abs(exported.modes-expected),[],'all')<2e-10);
assert(max(abs(exported.direction-K),[],'all')<1e-12);
assert(abs(r.image.globalComponentsSquared(3,3)-0.36)<1e-12);
rejected=false;
try
    prt_require_transverse([0;0;1],[0;0;1]);
catch problem
    rejected=contains(problem.message,'不垂直');
end
assert(rejected);
fprintf('PASS: spatial Jones/Stokes CSV, irregular sampling and paired plane results\n');

    function writeMap(values,columns)
        T=array2table([xy values],'VariableNames',[{'X_mm','Y_mm'},columns]);
        writetable(T,file);
    end
    function [data,count]=roundTrip(plane,format)
        count=prt_save_pupil(file,plane,format);
        data=prt_load_pupil(file);
    end
end
