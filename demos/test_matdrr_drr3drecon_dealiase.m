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


%% create data
wt=300;
wx=80;
p=2.2;
a1=zeros(wt,wx);
[n,m]=size(a1);
a2=a1;
a3=a1;
a4=a1;
a5=a1;
a6=a1;
a7=a1;
a8=a1;
a9=a1;
a10=a1;
a11=a1;
a12=a1;
a13=a1;
a14=a1;
k=0;
a=0.1;
b=1;
for t=-0.054:0.001:0.056
    k=k+1;
    b1(k)=(1-2*(pi*35*t).^2).*exp(-(pi*35*t).^2);
    b2(k)=(1-2*(pi*20*t).^2).*exp(-(pi*20*t).^2);
    b3(k)=(1-2*(pi*60*t).^2).*exp(-(pi*60*t).^2);
    b4(k)=(1-2*(pi*30*t).^2).*exp(-(pi*30*t).^2);
end
randn('state',1);
c=0.2*round(randn(1,200));
c=0.2*ones(1,200);
c([1:10:199])=1;
c([2:20:200])=-1;
%bnew=conv(c,b1);
bnew=sin(0:0.3:20*pi);
kk=size(bnew,2);
for j=1:20
for i=1:m
%     t1(i)=round(a*i^2+b);
%     t2(i)=round(0.5*i^2+100);
%     t3(i)=round(0.02*i^2+200);
%     t4(i)=round(5*sin(i/3)+350);
%     a1(t1(i):t1(i)+k-1,i)=b1;
%     a2(20+i:20+k-1+i,i)=b2;
%     a3(40+2*i:40+k-1+2*i,i)=b3;
%     a4(40+2*i:40+k-1+2*i,i)=b4;
%     a5(300-i:300+k-1-i,i)=b2;
%     a6(470-3*i:470+k-1-3*i,i)=b3;
%     a7(t2(i):t2(i)+k-1,i)=b3;
%     a8(t3(i):t3(i)+k-1,i)=b4;
%     a9(30:30+k-1,i)=b4;
%  a14(t4(i):t4(i)+k-1,i)=b2;

%   t1(i)=round(0.005*(i-1)^2+110);
%   t11(i)=1;
%   t2(i)=round(0.009*(i-1)^2+50);
%   t3(i)=round(0.003*(i-1)^2+170);
  t4(i)=round(-1*i+185-1*j);
  t5(i)=round(2*i+15);
  t6(i)=round(50+j*3);
%   a1(t1(i):t1(i)+k-1,i)=b1;
%   % a1(t11(i):t11(i)+k-1,i)=b1;
%   a2(t2(i):t2(i)+k-1,i)=b1; 
%   a3(t3(i):t3(i)+k-1,i)=b1; 
  a4(t4(i):t4(i)+k-1,i)=b2;
  a5(t5(i):t5(i)+k-1,i)=-b2;
  a6(t6(i):t6(i)+k-1,i)=b1;
end

%% add band-limitted noise and decimate data
randn('state',201315+j);
%shot=0*randn(wt,wx)+a2(1:wt,:)+a1(1:wt,:)+a3(1:wt,:)+a4(1:wt,:);%+a5(1:500,:)+a6(1:500,:)+a7(1:500,:)+a8(1:500,:)+a9(1:500,:)+a10(1:500,:)+a11(1:500,:)+a12(1:500,:)+a13(1:500,:)+a14(1:500,:);%+randn(100,20,10);
noise=1*randn(wt,wx);
% f_n=fft(noise);
% f_n(1:3,:)=0;
% f_n(15:100,:)=0;
% noise=real(ifft(f_n));
[o] =  yc_bp(noise,0.001,2,10,90,100);
shotn(:,:,j)=0.6*o+a5(1:wt,:)+a4(1:wt,:)+a6(1:wt,:);
shotc(:,:,j)=a5(1:wt,:)+a4(1:wt,:)+a6(1:wt,:);
end
shotn2=shotn;
jj=1;
for j=2:2:20
    shotn(:,:,j)=0;
    ii=1;
for i=2:2:wx
    shotn(:,i,j-1)=0;
    shotn_in(:,ii,jj)=shotn(:,i-1,j-1);
    ii=ii+1;
end
jj=jj+1;
end

%% Doing the reconstruction

dn=shotn2;
% dn=shotn2(51:150,1:40,1:20);
d3=drr3drecon_dealiase(dn,1,50,0.004,3,4,10,2,2,1);%takes about 20 minutes

figure;imagesc(dn(:,:,1));
figure;imagesc(d3(:,:,1));

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
subplot(1,2,1);drr_plot3d(dn,[100,10,5],zz,xx,yy);caxis([-1,1]);xlabel('X (sample)','Fontsize',15);ylabel('Y (sample)','Fontsize',15);zlabel('Time (s)','Fontsize',15);title('Raw','Fontsize',15,'fontweight','normal');set(gca,'Linewidth',2,'Fontsize',15);
subplot(1,2,2);drr_plot3d(d3,[100,20,10],z,x,y);caxis([-1,1]);xlabel('X (sample)','Fontsize',15);ylabel('Y (sample)','Fontsize',15);zlabel('Time (s)','Fontsize',15);title('Densified','Fontsize',15,'fontweight','normal');set(gca,'Linewidth',2,'Fontsize',15);

% subplot(3,2,5);drr_plot3d(diffr1,[100,10,10],z,x,y);caxis([-1,1]);xlabel('X (sample)','Fontsize',15);ylabel('Y (sample)','Fontsize',15);zlabel('Time (s)','Fontsize',15);title('LDRR diffraction','Fontsize',15,'fontweight','normal');set(gca,'Linewidth',2,'Fontsize',15);
% subplot(3,2,6);drr_plot3d(data-diffr,[100,10,10],z,x,y);caxis([-1,1]);xlabel('X (sample)','Fontsize',15);ylabel('Y (sample)','Fontsize',15);zlabel('Time (s)','Fontsize',15);title('LDRR reflection','Fontsize',15,'fontweight','normal');set(gca,'Linewidth',2,'Fontsize',15);
print(gcf,'-dpng','-r300','test_matdrr_drr3drecon_dealiase.png');
print(gcf,'-deps','-r300','test_matdrr_drr3drecon_dealiase.eps');



