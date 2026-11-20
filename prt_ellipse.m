function [psi,chi,S] = prt_ellipse(E)
% XY projection only; use a common transverse frame for tilted rays.
% S3=2 Im(conj(Ex) Ey); no observer-dependent handedness label implied.
Ex=E(1,:); Ey=E(2,:);
S=[abs(Ex).^2+abs(Ey).^2; abs(Ex).^2-abs(Ey).^2; ...
    2*real(conj(Ex).*Ey); 2*imag(conj(Ex).*Ey)];
psi=0.5*atan2(S(3,:),S(2,:));
chi=0.5*asin(max(-1,min(1,S(4,:)./S(1,:))));
zero=S(1,:)==0; psi(zero)=NaN; chi(zero)=NaN;
end
