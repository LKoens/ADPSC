function [Fq]  = ForceQ(h,x,y,N,fq)
%Force quadrapole Determines the folow at (x,y) from the two dimensional force
% quadrapole of strength fq in a doubly periodic array.
%  =======INPUTS=========
%   h       -   aspect ratio of the cell
%   x       -   x position for flow in the cell
%   y       -   y position for flow in the cell
%   N       -   Series trunkation number
%   fq      -   Singularity strength
%  =======OUTPUTS=========
%   FQ      -   Flow velocity represented as a complex number.

%Set initial values
z0 = (1+1i*h)/2;
z=(x+1i*y);
zeta=exp(2*pi*1i*(z-z0));
rho= exp(-2*pi*h); %defines rho

% s= sum((-1).^(-N:N).*sqrt(rho.^((-N:N).*((-N-1):(N-1)))).*(zeta.^(-N:N)));
% zs10= sum((-1).^(-N:N).*sqrt(rho.^((-N:N).*((-N-1):(N-1)))).*(-N:N).*(zeta.^((-N:N))));
% z2s20= sum((-1).^(-N:N).*sqrt(rho.^((-N:N).*((-N-1):(N-1)))).*(-N:N).*((-N:N)-1).*(zeta.^((-N:N))));
% z3s30= sum((-1).^(-N:N).*sqrt(rho.^((-N:N).*((-N-1):(N-1)))).*(-N:N).*((-N:N)-1).*((-N:N)-2).*(zeta.^((-N:N))));
% 
% rs01= sum((-1).^(-N:N).*sqrt(rho.^((-N:N).*((-N-1):(N-1)))).*((-N:N).*((-N-1):(N-1))).*(zeta.^((-N:N)))/2);
% rzs11= sum((-1).^(-N:N).*sqrt(rho.^((-N:N).*((-N-1):(N-1)))).*((-N:N).*((-N-1):(N-1))).*(-N:N).*(zeta.^((-N:N)))/2);
% rz2s21= sum((-1).^(-N:N).*sqrt(rho.^((-N:N).*((-N-1):(N-1)))).*((-N:N).*((-N-1):(N-1))).*(-N:N).*((-N:N)-1).*(zeta.^((-N:N)))/2);

s= sum((-1).^(-N:N).*(sqrt(rho.^(((-N-1):(N-1))))*zeta).^(-N:N));
zs10= sum((-1).^(-N:N).*((sqrt(rho.^(((-N-1):(N-1))))*zeta).^(-N:N)).*(-N:N));
z2s20= sum((-1).^(-N:N).*((sqrt(rho.^(((-N-1):(N-1))))*zeta).^(-N:N)).*(-N:N).*((-N:N)-1));
z3s30= sum((-1).^(-N:N).*((sqrt(rho.^(((-N-1):(N-1))))*zeta).^(-N:N)).*(-N:N).*((-N:N)-1).*((-N:N)-2));

rs01= sum((-1).^(-N:N).*((sqrt(rho.^(((-N-1):(N-1))))*zeta).^(-N:N)).*((-N:N).*((-N-1):(N-1)))/2);
rzs11= sum((-1).^(-N:N).*((sqrt(rho.^(((-N-1):(N-1))))*zeta).^(-N:N)).*((-N:N).*((-N-1):(N-1))).*(-N:N)/2);
rz2s21= sum((-1).^(-N:N).*((sqrt(rho.^(((-N-1):(N-1))))*zeta).^(-N:N)).*((-N:N).*((-N-1):(N-1))).*(-N:N).*((-N:N)-1)/2);


L=zs10/s + z2s20/(s) - (zs10/s)^2;
M=(s-zs10)*(s*(zs10+3*z2s20)-2*zs10^2)/s^3 + z3s30/s;
rLr=(rzs11+rz2s21)/s - (2*zs10*rzs11+rs01*(zs10+z2s20) )/s^2 + 2*rs01*zs10^2/s^3;

modz=log(abs(zeta)^2);

Fq = -conj(fq)*conj(L) - 2*fq*L-fq*modz*M - fq*log(rho^2)*rLr;

%keyboard

end

