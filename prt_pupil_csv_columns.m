function [xy,jones,stokes,vector] = prt_pupil_csv_columns()
% Shared schema for polarization maps.
xy={'X_mm','Y_mm'};
jones={'ExRe','ExIm','EyRe','EyIm'};
stokes={'S0','S1','S2','S3'};
vector=[jones,{'EzRe','EzIm','Kx','Ky','Kz'}];
end
