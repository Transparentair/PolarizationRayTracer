function [x,y,signS3,arrowX,arrowY] = prt_pupil_ellipses(pupil,S,valid,sampling)
% Draw at physical plane coordinates, with ray selection on the launch grid.
% Columns are individual contours. Major semiaxis is a fixed display length.
assert(size(pupil,1)==2 && size(S,1)==4 && ...
    size(pupil,2)==size(S,2) && numel(valid)==size(S,2));
if nargin<4, sampling=pupil; end
linear=hypot(S(2,:),S(3,:));
polarized=hypot(linear,S(4,:));
usable=reshape(valid,1,[]) & all(isfinite(pupil),1);
eligible=usable & all(isfinite(S),1) & S(1,:)>0 & polarized>1e-12*S(1,:);
extent=max(sampling,[],2)-min(sampling,[],2);
bins=round(12*(sampling-min(sampling,[],2))./max(extent,eps));
ids=find(eligible);
[~,pick]=unique(bins(:,ids).','rows','stable');
keep=ids(pick);
psi=0.5*atan2(S(3,keep),S(2,keep));
chi=0.5*atan2(S(4,keep),linear(keep));
signS3=sign(S(4,keep));
signS3(abs(S(4,keep))<=1e-12*polarized(keep))=0;

limits=pupil(:,usable);
if isempty(limits)
    scale=0;
else
    scale=max(max(limits,[],2)-min(limits,[],2));
end
if scale<=1e-12
    x=[]; y=[]; signS3=[]; arrowX=[]; arrowY=[]; return
end
a=0.38*scale/12;
t=linspace(0,2*pi,65).';
u=a*cos(t);
v=a*(sin(t)*tan(chi));
x=pupil(1,keep)+u*cos(psi)-v.*sin(psi);
y=pupil(2,keep)+u*sin(psi)+v.*cos(psi);
if nargout<4, return; end

arrows=nan(8,numel(psi));
major=exp(1i*psi);
elliptical=signS3~=0;
% At t=pi/2 the contour tangent for increasing t is -major.
% The signed minor axis already encodes S3. Time runs toward decreasing
% contour parameter for exp(+i*omega*t), for either sign of S3.
if any(elliptical)
    tip=complex(x(17,elliptical),y(17,elliptical));
    arrows(1:4,elliptical)=arrowHead(tip,major(elliptical),a);
end
linear=~elliptical;
if any(linear)
    tipPositive=complex(x(1,linear),y(1,linear));
    tipNegative=complex(x(33,linear),y(33,linear));
    arrows(1:4,linear)=arrowHead(tipPositive,major(linear),a);
    arrows(5:8,linear)=arrowHead(tipNegative,-major(linear),a);
end
arrowX=real(arrows);
arrowY=imag(arrows);
end

function points=arrowHead(tip,direction,a)
halfWidth=0.075*a;
base=tip-0.21*a*direction;
wing=1i*halfWidth*direction;
points=[base+wing;tip;base-wing;complex(nan(size(tip)),nan(size(tip)))];
end
