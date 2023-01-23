function [So]  = SO(h,x,y,N,so)
%SO Determines the folow at (x,y) from the two dimensional source octopole
% of strength so in a doubly periodic array.
%  =======INPUTS=========
%   h       -   aspect ratio of the cell
%   x       -   x position for flow in the cell
%   y       -   y position for flow in the cell
%   N       -   Series trunkation number
%   so      -   Singularity strength
%  =======OUTPUTS=========
%   So      -   Flow velocity represented as a complex number.


%Set initial values
z0 = (1+1i*h)/2;
z=(x+1i*y);
zeta=exp(2*pi*1i*(z-z0));
rho= exp(-2*pi*h); %defines rho

% s= sum((-1).^(-N:N).*sqrt(rho.^((-N:N).*((-N-1):(N-1)))).*(zeta.^(-N:N)));
% zs10= sum((-1).^(-N:N).*sqrt(rho.^((-N:N).*((-N-1):(N-1)))).*(-N:N).*(zeta.^((-N:N))));
% z2s20= sum((-1).^(-N:N).*sqrt(rho.^((-N:N).*((-N-1):(N-1)))).*(-N:N).*((-N:N)-1).*(zeta.^((-N:N))));
% z3s30= sum((-1).^(-N:N).*sqrt(rho.^((-N:N).*((-N-1):(N-1)))).*(-N:N).*((-N:N)-1).*((-N:N)-2).*(zeta.^((-N:N))));
% z4s40= sum((-1).^(-N:N).*sqrt(rho.^((-N:N).*((-N-1):(N-1)))).*(-N:N).*((-N:N)-1).*((-N:N)-2).*((-N:N)-3).*(zeta.^((-N:N))));


s= sum((-1).^(-N:N).*((sqrt(rho.^(((-N-1):(N-1))))*zeta).^(-N:N)));
zs10= sum((-1).^(-N:N).*((sqrt(rho.^(((-N-1):(N-1))))*zeta).^(-N:N)).*(-N:N));
z2s20= sum((-1).^(-N:N).*((sqrt(rho.^(((-N-1):(N-1))))*zeta).^(-N:N)).*(-N:N).*((-N:N)-1));
z3s30= sum((-1).^(-N:N).*((sqrt(rho.^(((-N-1):(N-1))))*zeta).^(-N:N)).*(-N:N).*((-N:N)-1).*((-N:N)-2));
z4s40= sum((-1).^(-N:N).*((sqrt(rho.^(((-N-1):(N-1))))*zeta).^(-N:N)).*(-N:N).*((-N:N)-1).*((-N:N)-2).*((-N:N)-3));

n= zs10/s -7*(zs10/s)^2+12*(zs10/s)^3-6*(zs10/s)^4 +7*z2s20/s-18*zs10*z2s20/s^2+12*(z2s20/s)*(zs10/s)^2-3*(z2s20/s)^2 +6*z3s30/s -4*zs10*z3s30/s^2+z4s40/s;

So = so*n;

%keyboard

end

