function [u,ux, uy] = velocity_in_cell(R,h,x,y)
%ASYMPTOTIC_DOUBLY_PERIODIC provides background velocity, force, pressure
%drop, flux, permeability, and velocity for a post of size R in a cell of
% size 1xh
%  =======INPUTS=========
%   R       -   Radius of the post
%   h       -   aspect ratio of the cell
%   x       -   x matrix of positions for flow in the cell
%   y       -   y matrix of positions for flow in the cell
%  =======OUTPUTS=========
%   u        -    complex velocity in the form ux - i uy
%   ux       -    flow velocity in x at x,y in the cell
%   uy       -    flow velocity in y at x,y in the cell


%% Determine cell properties
N=100;
[U,F2,~,~,~,sd,fq,so]= Asymptotic_Doubly_periodic(R,h);
F=-F2/(8*pi);

dim = size(x);
u= zeros(dim);

%% Flow velocity

for i=1:dim(1)
    for j=1:dim(2)
        if (x(i,j)-1/2)^2+(y(i,j)-h/2)^2<R^2
            disp('Point selected within cylinder. Velocity set to 0')
            u(i,j)=0;
        else
            u(i,j)=U+Stokeslet(h,x(i,j),y(i,j),N,F) + SD(h,x(i,j),y(i,j),N,sd) +ForceQ(h,x(i,j),y(i,j),N,fq)+SO(h,x(i,j),y(i,j),N,so);
        end
    end
end
ux=real(u);
uy=-imag(u);


end

