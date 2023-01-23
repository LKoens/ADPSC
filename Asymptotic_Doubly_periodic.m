function [U,F2,DeltaP,ubar,k,sd,fq,so] = Asymptotic_Doubly_periodic(R,h)
%ASYMPTOTIC_DOUBLY_PERIODIC provides background velocity, force, pressure
%drop, flux, permeability, and velocity for a post of size R in a cell of
% size 1xh
%  =======INPUTS=========
%   R       -   Radius of the post
%   h       -   aspect ratio of the cell
%  =======OUTPUTS=========
%   U       -    Background flow in the absence of the post
%   F       -    Force per unit length in z
%   DeltaP  -    Pressure drop over the cell
%   ubar    -    mean flow of fluid throught the cell
%   k       -    Permeability
%   sd      -    Source dipole strength
%   fq      -    Force quadrapole strength
%   so      -    source octopole strength
%% Warnings

if R>0.15
disp('Radius too large. Solution may have >5% error')
end

if (R/h)>0.15
disp('Cell confinement h too small. Solution may have >5% error')
end

%% Set initial values
N=100;
rho= exp(-2*pi*h); %defines rho
A=prod((1+rho.^(1:N)).^2)/sum(sqrt(rho.^((1:N).*(0:(N-1))))); %Define A matrix

%Derivative of A with respect to rho
Ap =-prod((1+rho.^(1:N)).^2)*sum(0.5*(1:N).*(0:(N-1)).*sqrt(rho.^((1:N).*(0:(N-1))))/rho)/(sum(sqrt(rho.^((1:N).*(0:(N-1)))))^2); 
for i=1:N
temp= (1+rho.^(1:N)).^2;
temp(i) = 2*i*rho^(i-1)*(1+rho^i);
Ap = Ap + prod(temp)/sum(sqrt(rho.^((1:N).*(0:(N-1))))); 
end 

s10= sum((-1).^(-N:N).*sqrt(rho.^((-N:N).*((-N-1):(N-1)))).*(-N:N));
rs11= sum((-1).^(-N:N).*sqrt(rho.^((-N:N).*((-N-1):(N-1)))).*((-N:N).*((-N-1):(N-1))).*(-N:N)/2);
s30= sum((-1).^(-N:N).*sqrt(rho.^((-N:N).*((-N-1):(N-1)))).*(-N:N).*((-N:N)-1).*((-N:N)-2));
rs31= sum((-1).^(-N:N).*sqrt(rho.^((-N:N).*((-N-1):(N-1)))).*((-N:N).*((-N-1):(N-1))).*(-N:N).*((-N:N)-1).*((-N:N)-2)/2);

B = 1+4*s30/s10;

S00 = 1+ log(4*pi^2*R^2*A^2*s10^2) + 2*log(rho)*(rho*Ap/A + rs11/s10);
S20 = pi^2*(B + 12/log(rho))/3;
S2m2= S20/2;
S22 = pi^2*(3*B + 12/log(rho) +4*log(rho^2)*(rs31/s10 -rs11*s30/(s10)^2))/6;

CFU  = 3/(3*S00+(R^2)*(pi^2*B+3*S20));
CSDF = 4*pi^2*(1+R^2*S2m2);
CFQF = 4*pi^2*S22;
CSOFQ= -4*pi^2/3;

%keyboard

%% Singularity strengths

U=1/(1+Stokeslet(h,1/2,h,N,CFU) + SD(h,1/2,h,N,R^2*CSDF*CFU) + ForceQ(h,1/2,h,N,R^4*CFQF*CFU)+SO(h,1/2,h,N,R^6*CSOFQ*CFQF*CFU)); %background flow without the post
F=CFU*U; %force on the post from the fluid
sd= R^2*CSDF*F;
fq=R^4*CFQF*F;
so=R^6*CSOFQ*CFQF*F;
F2= -8*pi*F;
DeltaP = -4*F/h; %pressure drop over cell

%% mean flow 
zeta=-sqrt(rho);
sv2= sum((-1).^(-N:N).*((sqrt(rho.^(((-N-1):(N-1))))*zeta).^(-N:N)));
zs10v2=  sum((-1).^(-N:N).*((sqrt(rho.^(((-N-1):(N-1))))*zeta).^(-N:N)).*(-N:N));
P=A*sv2; %Schottky-Klien prime function for an anulus.
K=zs10v2/sv2;
Sflux= log(rho)^2/(12*pi) + 2*log(rho)*log(P/(rho^(1/4)))/(2*pi) - log(rho^2)*sum((1:N).*log((1+rho.^((1:N)-1/2))./(1+rho.^((1:N) +1/2))))/pi;
SDflux= (2*K-1)/(2*pi);

ubar = U + Sflux*F/h + SDflux*(sd+fq)/h;

%% Permiability

k =U/DeltaP - Sflux/4 - SDflux*(R^2*CSDF+R^4*CFQF)/4;


end

