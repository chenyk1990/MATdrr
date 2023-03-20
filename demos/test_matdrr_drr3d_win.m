% Demonstration script for 
% 3D seismic denoising via the localized damped rank-reduction method
%
%  Copyright (C) 2015 The University of Texas at Austin
%  Copyright (C) 2015 Yangkang Chen
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
%  Wang et al., 2020, Separation and imaging of seismic diffractions using a localized rank-reduction method with adaptively selected ranks, 85, V497–V506.
%  Chen et al., 2022, 3D seismic diffraction separation and imaging using the local rank-reduction method, IEEE Transactions on Geoscience and Remote Sensing, 60, 4507110.
%  Chen et al., 2019, Obtaining free USArray data by multi-dimensional seismic reconstruction, Nature Communications, 10:4434.
%  Bai et al., 2020, Seismic signal enhancement based on the lowrank methods, Geophysical Prospecting, 68, 2783-2807.
%  Chen et al., 2020, Five-dimensional seismic data reconstruction using the optimally damped rank-reduction method, Geophysical Journal International, 222, 1824-1845.
%  Chen, Y., W. Huang, D. Zhang, W. Chen, 2016, An open-source matlab code package for improved rank-reduction 3D seismic data denoising and reconstruction, Computers & Geosciences, 95, 59-66.
%  Chen, Y., D. Zhang, Z. Jin, X. Chen, S. Zu, W. Huang, and S. Gan, 2016, Simultaneous denoising and reconstruction of 5D seismic data via damped rank-reduction method, Geophysical Journal International, 206, 1695-1717.
%  Huang, W., R. Wang, Y. Chen, H. Li, and S. Gan, 2016, Damped multichannel singular spectrum analysis for 3D random noise attenuation, Geophysics, 81, V261-V270.
%  Chen et al., 2017, Preserving the discontinuities in least-squares reverse time migration of simultaneous-source data, Geophysics, 82, S185-S196.
% 

clc;clear;close all;
addpath(genpath('../matdrr'));

%% generate synthetic data
cmp=zeros(50,32,32);
[nz,nx,ny]=size(cmp);

x=[1:32];
y=[1:32];
[x,y]=meshgrid(x,y);
z=10*(sqrt(1+((x-16)./10).^2+((y-16)./10).^2)-1)+20;
% figure;surf(y);set(gca,'ydir','reverse');

for ix=1:nx
    for iy=1:ny
        cmp(round(z(ix,iy))+1,ix,iy)=1;
    end  
end

cmp=reshape(cmp,nz,nx*ny);
wav=drr_ricker(30,0.004,0.2);
for ix=1:nx*ny
cmp(:,ix)=conv(cmp(:,ix),wav,'same');
end
cmp=reshape(cmp,nz,nx,ny);
dc=drr_scale(cmp,3);

randn('state',201920);
dn=dc+0.1*randn(size(dc));

figure;imagesc(dc(:,:,16));colormap(gray);
figure;imagesc(squeeze(dc(:,16,:)));colormap(gray);
figure;imagesc(squeeze(dc(:,:)));colormap(gray);
figure;imagesc(squeeze(dn(:,:)));colormap(gray);

%% parameters
[n1,n2,n3]=size(dn);
param.dt=0.004;
param.flow=0;
param.fhigh=245;
param.N=2;
param.K=4;
param.verb=1;

%% Global DRR
param.N=6;
d1=drr3d(dn,param.flow,param.fhigh,param.dt,param.N,param.K,param.verb);

%% Local DRR (Option I, integrated version)
n1win=50;n2win=16;n3win=16;
r1=0.5;r2=0.5;r3=0.5;
param.N=2;
d2=drr3d_win(dn,param.flow,param.fhigh,param.dt,param.N,param.K,param.verb,n1win,n2win,n3win,r1,r2,r3);

%% Local DRR (Option II, separated version)
d3=drr_win3d(@localdrr3d, param, dn, n1win, n2win, n3win, r1, r2, r3);

%% SNR comparison
s0=drr_snr(dc,dn,2);
s1=drr_snr(dc,d1,2);
s2=drr_snr(dc,d2,2)
fprintf('SNR of Noisy is %g\n',drr_snr(dc,dn,2));
fprintf('SNR of DRR is %g\n',drr_snr(dc,d1,2));
fprintf('SNR of LDRR is %g\n',drr_snr(dc,d2,2));

%% verify correctness of two options (output is zero)
fprintf('Difference between two versions is %g\n',norm(d2(:)-d3(:)));


%% plot the results
[n1,n2,n3]=size(dn);
dy=1;dx=1;dt=0.004;
y=[1:n3]*dy;
x=[1:n2]*dx;
z=[1:n1]*dt;
figure('units','normalized','Position',[0.2 0.4 0.8, 0.6],'color','w');
subplot(2,3,1);drr_plot3d(dn,[25,16,16],z,x,y);caxis([-0.4,0.4]);xlabel('X (sample)','Fontsize',15);ylabel('Y (sample)','Fontsize',15);zlabel('Time (s)','Fontsize',15);title(strcat('Noisy (SNR=',num2str(s0),' dB )'),'Fontsize',15,'fontweight','normal');set(gca,'Linewidth',2,'Fontsize',15);text(-1.2,-0.2, -0.06,'a)','color','k','Fontsize',30,'fontweight','bold','HorizontalAlignment','left');
subplot(2,3,2);drr_plot3d(d1,[25,16,16],z,x,y);caxis([-0.4,0.4]);xlabel('X (sample)','Fontsize',15);ylabel('Y (sample)','Fontsize',15);zlabel('Time (s)','Fontsize',15);title(strcat('DRR (SNR=',num2str(s1),' dB)'),'Fontsize',15,'fontweight','normal');set(gca,'Linewidth',2,'Fontsize',15);text(-1.2,-0.2, -0.06,'b)','color','k','Fontsize',30,'fontweight','bold','HorizontalAlignment','left');
subplot(2,3,3);drr_plot3d(d2,[25,16,16],z,x,y);caxis([-0.4,0.4]);xlabel('X (sample)','Fontsize',15);ylabel('Y (sample)','Fontsize',15);zlabel('Time (s)','Fontsize',15);title(strcat('LDRR (SNR=',num2str(s2),' dB)'),'Fontsize',15,'fontweight','normal');set(gca,'Linewidth',2,'Fontsize',15);text(-1.2,-0.2, -0.06,'c)','color','k','Fontsize',30,'fontweight','bold','HorizontalAlignment','left');
subplot(2,3,5);drr_plot3d(dn-d1,[25,16,16],z,x,y);caxis([-0.4,0.4]);xlabel('X (sample)','Fontsize',15);ylabel('Y (sample)','Fontsize',15);zlabel('Time (s)','Fontsize',15);title(strcat('Noise of DRR'),'Fontsize',15,'fontweight','normal');set(gca,'Linewidth',2,'Fontsize',15);text(-1.2,-0.2, -0.06,'d)','color','k','Fontsize',30,'fontweight','bold','HorizontalAlignment','left');
subplot(2,3,6);drr_plot3d(dn-d2,[25,16,16],z,x,y);caxis([-0.4,0.4]);xlabel('X (sample)','Fontsize',15);ylabel('Y (sample)','Fontsize',15);zlabel('Time (s)','Fontsize',15);title(strcat('Noise of LDRR'),'Fontsize',15,'fontweight','normal');set(gca,'Linewidth',2,'Fontsize',15);text(-1.2,-0.2, -0.06,'e)','color','k','Fontsize',30,'fontweight','bold','HorizontalAlignment','left');
print(gcf,'-dpng','-r300','test_matdrr_drr3d_win.png');
print(gcf,'-depsc','-r200','test_matdrr_drr3d_win.eps');








