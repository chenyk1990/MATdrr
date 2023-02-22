% Demonstration script for 
% 3D dealiased reconstruction via the damped rank-reduction method
%
%  Copyright (C) 2023 The University of Texas at Austin
%  Copyright (C) 2023 Yangkang Chen
%
%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published
%  by the Free Software Foundation, either version 3 of the License, or
%  any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details: http://www.gnu.org/licenses/
%
%  References:   
%
%  [0] Huang, W., D. Feng, and Y. Chen, 2020, De‐aliased and de‐noise Cadzow filtering for seismic data reconstruction, Geophysical Prospecting, 68, 443-571.
%  [1] Chen, Y., W. Huang, D. Zhang, W. Chen, 2016, An open-source matlab code package for improved rank-reduction 3D seismic data denoising and reconstruction, Computers & Geosciences, 95, 59-66.
%  [2] Chen, Y., D. Zhang, Z. Jin, X. Chen, S. Zu, W. Huang, and S. Gan, 2016, Simultaneous denoising and reconstruction of 5D seismic data via damped rank-reduction method, Geophysical Journal International, 206, 1695-1717.
%  [3] Huang, W., R. Wang, Y. Chen, H. Li, and S. Gan, 2016, Damped multichannel singular spectrum analysis for 3D random noise attenuation, Geophysics, 81, V261-V270.
%  [4] Chen et al., 2017, Preserving the discontinuities in least-squares reverse time migration of simultaneous-source data, Geophysics, 82, S185-S196.
%  [5] Chen et al., 2019, Obtaining free USArray data by multi-dimensional seismic reconstruction, Nature Communications, 10:4434.

clc;clear;close all;
addpath(genpath('../matdrr'));


%% generate 3D synthetic data
a1=zeros(300,20);
[n,m]=size(a1);
a3=a1;
a4=a1;

k=0;
a=0.1;
b=1;
for t=-0.055:0.002:0.055
    k=k+1;
    b1(k)=(1-2*(pi*30*t).^2).*exp(-(pi*30*t).^2);
    b2(k)=(1-2*(pi*40*t).^2).*exp(-(pi*40*t).^2);
    b3(k)=(1-2*(pi*40*t).^2).*exp(-(pi*40*t).^2);
    b4(k)=(1-2*(pi*30*t).^2).*exp(-(pi*30*t).^2);
end
for i=1:m
  t1(i)=round(140);
  t3(i)=round(-6*i+180);
  t4(i)=round(6*i+10);
  a1(t1(i):t1(i)+k-1,i)=b1; 
  a3(t3(i):t3(i)+k-1,i)=b1; 
  a4(t4(i):t4(i)+k-1,i)=b1;
end

temp=a1(1:300,:)+a3(1:300,:)+a4(1:300,:);
for j=1:20
    a4=zeros(300,20);
    for i=1:m
  t4(i)=round(6*i+10+3*j); 
  a4(t4(i):t4(i)+k-1,i)=b1;
  
  t1(i)=round(140-2*j);
  a1(t1(i):t1(i)+k-1,i)=b1;
    end
    shot(:,:,j)=a1(1:300,:)+a3(1:300,:)+a4(1:300,:);
end
plane3d=shot;
d=plane3d/max(max(max(plane3d)));


%% add band-limitted noise and decimate data
randn('state',201315);
[nt,nx,ny]=size(d);
noise=randn(nt,nx,ny);
noise=drr_bandpass(noise,0.004,0,60);
dn=0.2*noise+d;

%% Doing the aliased reconstruction (densification)
d3=drr3drecon_dealiase(dn,1,50,0.004,3,4,10,2,2,1);%takes about 20 minutes

dy=1;dx=1;dt=0.004;
[n1,n22,n33]=size(dn);
yy=[1:n33]*dy;
xx=[1:n22]*dx;
zz=[1:n1]*dt;

[n1,n2,n3]=size(d3);
dy=1;dx=1;dt=0.004;
y=[1:n3]*dy;
x=[1:n2]*dx;
z=[1:n1]*dt;
figure('units','normalized','Position',[0.2 0.4 0.8, 0.6],'color','w');
subplot(1,2,1);drr_plot3d(dn,[100,5,5],zz,xx,yy);caxis([-0.5,0.5]);xlabel('X (sample)','Fontsize',15);ylabel('Y (sample)','Fontsize',15);zlabel('Time (s)','Fontsize',15);title('Raw','Fontsize',15,'fontweight','normal');set(gca,'Linewidth',2,'Fontsize',15);text(-2.5,-2.5, -0.3,'a)','color','k','Fontsize',40,'fontweight','bold','HorizontalAlignment','left');
subplot(1,2,2);drr_plot3d(d3,[100,10,10],z,x,y);caxis([-0.5,0.5]);xlabel('X (sample)','Fontsize',15);ylabel('Y (sample)','Fontsize',15);zlabel('Time (s)','Fontsize',15);title('Densified','Fontsize',15,'fontweight','normal');set(gca,'Linewidth',2,'Fontsize',15);text(-5,-5, -0.3,'b)','color','k','Fontsize',40,'fontweight','bold','HorizontalAlignment','left');

% subplot(3,2,5);drr_plot3d(diffr1,[100,10,10],z,x,y);caxis([-0.5,0.5]);xlabel('X (sample)','Fontsize',15);ylabel('Y (sample)','Fontsize',15);zlabel('Time (s)','Fontsize',15);title('LDRR diffraction','Fontsize',15,'fontweight','normal');set(gca,'Linewidth',2,'Fontsize',15);
% subplot(3,2,6);drr_plot3d(data-diffr,[100,10,10],z,x,y);caxis([-0.5,0.5]);xlabel('X (sample)','Fontsize',15);ylabel('Y (sample)','Fontsize',15);zlabel('Time (s)','Fontsize',15);title('LDRR reflection','Fontsize',15,'fontweight','normal');set(gca,'Linewidth',2,'Fontsize',15);
print(gcf,'-dpng','-r300','test_matdrr_drr3drecon_dealiase.png');
print(gcf,'-depsc','-r200','test_matdrr_drr3drecon_dealiase.eps');



