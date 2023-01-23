function [S]  = Stokeslet(h,x,y,N,F)
%Stokeslet Determines the folow at (x,y) from the two dimensional stokeslet
% of strength F in a doubly periodic array.
%  =======INPUTS=========
%   h       -   aspect ratio of the cell
%   x       -   x position for flow in the cell
%   y       -   y position for flow in the cell
%   N       -   Series trunkation number
%   F      -   Singularity strength
%  =======OUTPUTS=========
%   S      -   Flow velocity represented as a complex number.

%Set initial values
z0 = (1+1i*h)/2;
z=(x+1i*y);
zeta=exp(2*pi*1i*(z-z0));
rho= exp(-2*pi*h); %defines rho
A=prod((1+rho.^(1:N)).^2)/sum(sqrt(rho.^((1:N).*(0:(N-1))))); %Define A

%Derivative of A with respect to rho
Ap =-prod((1+rho.^(1:N)).^2)*sum(0.5*(1:N).*(0:(N-1)).*sqrt(rho.^((1:N).*(0:(N-1))))/rho)/(sum(sqrt(rho.^((1:N).*(0:(N-1)))))^2); 
for i=1:N
temp= (1+rho.^(1:N)).^2;
temp(i) = 2*i*rho^(i-1)*(1+rho^i);
Ap = Ap + prod(temp)/sum(sqrt(rho.^((1:N).*(0:(N-1))))); 
end 

% s= sum((-1).^(-N:N).*sqrt(rho.^((-N:N).*((-N-1):(N-1)))).*(zeta.^(-N:N)));
% zs10= sum((-1).^(-N:N).*sqrt(rho.^((-N:N).*((-N-1):(N-1)))).*(-N:N).*(zeta.^((-N:N))));
% rs01= sum((-1).^(-N:N).*sqrt(rho.^((-N:N).*((-N-1):(N-1)))).*((-N:N).*((-N-1):(N-1))).*(zeta.^((-N:N)))/2);

s= sum((-1).^(-N:N).*((sqrt(rho.^(((-N-1):(N-1))))*zeta).^(-N:N)));
zs10= sum((-1).^(-N:N).*((sqrt(rho.^(((-N-1):(N-1))))*zeta).^(-N:N)).*(-N:N));
rs01= sum((-1).^(-N:N).*((sqrt(rho.^(((-N-1):(N-1))))*zeta).^(-N:N)).*((-N:N).*((-N-1):(N-1)))/2);

P=A*s; %Schottky-Klien prime function for an anulus.
K=zs10/s;
rPr=rho*Ap/A + rs01/s;

modz=log(abs(zeta)^2);

S = -2*conj(F)*log(abs(P)) + real(F)*modz-F*modz*K-F*log(rho^2)*rPr-real(F)*modz^2/(2*log(rho)); 

%keyboard

end

