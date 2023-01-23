function [Sd]  = SD(h,x,y,N,sd)
%SD Determines the folow at (x,y) from the two dimensional source dipole of
% strength sd in a doubly periodic array.
%  =======INPUTS=========
%   h       -   aspect ratio of the cell
%   x       -   x position for flow in the cell
%   y       -   y position for flow in the cell
%   N       -   Series trunkation number
%   sd      -   Singularity strength
%  =======OUTPUTS=========
%   Sd      -   Flow velocity represented as a complex number.


%Set initial values
z0 = (1+1i*h)/2;
z=(x+1i*y);
zeta=exp(2*pi*1i*(z-z0));
rho= exp(-2*pi*h); %defines rho

% s= sum((-1).^(-N:N).*sqrt(rho.^((-N:N).*((-N-1):(N-1)))).*(zeta.^(-N:N)));
% zs10= sum((-1).^(-N:N).*sqrt(rho.^((-N:N).*((-N-1):(N-1)))).*(-N:N).*(zeta.^((-N:N))));
% z2s20= sum((-1).^(-N:N).*sqrt(rho.^((-N:N).*((-N-1):(N-1)))).*(-N:N).*((-N:N)-1).*(zeta.^((-N:N))));

s= sum((-1).^(-N:N).*((sqrt(rho.^(((-N-1):(N-1))))*zeta).^(-N:N)));
zs10= sum((-1).^(-N:N).*((sqrt(rho.^(((-N-1):(N-1))))*zeta).^(-N:N)).*(-N:N));
z2s20= sum((-1).^(-N:N).*((sqrt(rho.^(((-N-1):(N-1))))*zeta).^(-N:N)).*(-N:N).*((-N:N)-1));

L=zs10/s + z2s20/(s) - (zs10/s)^2;

Sd = -sd*L;

%keyboard

end

