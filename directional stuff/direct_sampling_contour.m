function [xcont, ycont]=direct_sampling_contour(x,y,P,dtdeg)

% calculates direct sampling contour
% Inputs:
%   x,y - sample of data
%   P - non-exceedance probability of contour
%   dtdeg - directional step in degrees

% initialisation
dt=dtdeg*pi/180;
t=(dt:dt:2*pi)';
Nx=length(x);
Nt=length(t);
Np=length(P);
r=zeros(Nt,Np);
xcont=zeros(Nt,Np);
ycont=zeros(Nt,Np);

% find radius for each angle
for i=1:length(t)
    if Nx>=1e6
        fprintf('   angle %u/%u\n',i,Nt)
    end
    z=x*cos(t(i))+y*sin(t(i));
    r(i,:)=quantile(z,P);
end

% find intersection of lines
t=[t;dt];
r=[r;r(1,:)];
denominator = sin(t(2:end)).*cos(t(1:end-1)) - sin(t(1:end-1)).*cos(t(2:end));
for j=1:Np
    xcont(:,j)=(sin(t(2:end)).*r(1:end-1,j)-sin(t(1:end-1)).*r(2:end,j))./denominator;
    ycont(:,j)=(-cos(t(2:end)).*r(1:end-1,j)+cos(t(1:end-1)).*r(2:end,j))./denominator;
end