% Demonstration script for 
% 3D diffraction separation via the localized damped rank-reduction method
% (LDRR)
% 
%  Copyright (C) 2020 Yangkang Chen and Wei Chen
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
%  Chen et al., 2022, 3D seismic diffraction separation and imaging using the local rank-reduction method, IEEE Transactions on Geoscience and Remote Sensing, 60, 4507110.
%  Wang et al., 2020, Separation and imaging of seismic diffractions using a localized rank-reduction method with adaptively selected ranks, 85, V497–V506.
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


%% Please download data from https://github.com/aaspip/data/blob/main/diffr_syn_3d.mat
% data is observed data
% diffr is the ground-truth diffraction data
load diffr_syn_3d.mat

[n1,n2,n3]=size(data);

%% perform diffraction separation
%This example is introduced in Chen et al., 2022.
lf=0;hf=100;dt=0.004;verb=1;N=12;N2=50;K=4; % N is a scalar (Nmax) or a vector (Nmin,Nmax);
n1win=200;n2win=40;n3win=40;r1=0.5;r2=0.5;r3=0.5;mode=2;%mode=2 means using the singular value ratio criterion
d1=drr3d_win_auto(data,lf,hf,dt,[N,N2],K,verb,n1win,n2win,n3win,r1,r2,r3,mode);%reflection from LDRR 
diffr1=data-d1; %diffraction from LDRR 

% Quick visualization
figure;drr_imagesc([data(:,:,20),diffr(:,:,20),diffr1(:,:,20),data(:,:,20)-diffr(:,:,20),d1(:,:,20)]);

%% SNR of the separated diffraction
drr_snr(diffr,diffr1,2)

%% plot the results
y=[0:n3-1]*dy;
x=[0:n2-1]*dx;
z=[0:n1-1]*dt;
figure('units','normalized','Position',[0.2 0.4 0.5, 0.8],'color','w');
subplot(3,2,1);drr_plot3d(data,[100,20,20],z,x,y);caxis([-0.5,0.5]);xlabel('X (km)','Fontsize',15);ylabel('Y (km)','Fontsize',15);zlabel('Time (s)','Fontsize',15);title('Data','Fontsize',15,'fontweight','normal');set(gca,'Linewidth',2,'Fontsize',15);text(-0.5,-0.2, -0.3,'a)','color','k','Fontsize',30,'fontweight','bold','HorizontalAlignment','left');
subplot(3,2,3);drr_plot3d(diffr,[100,20,20],z,x,y);caxis([-0.5,0.5]);xlabel('X (km)','Fontsize',15);ylabel('Y (km)','Fontsize',15);zlabel('Time (s)','Fontsize',15);title('Ground-truth diffraction','Fontsize',15,'fontweight','normal');set(gca,'Linewidth',2,'Fontsize',15);text(-0.5,-0.2, -0.3,'b)','color','k','Fontsize',30,'fontweight','bold','HorizontalAlignment','left');
subplot(3,2,4);drr_plot3d(data-diffr,[100,20,20],z,x,y);caxis([-0.5,0.5]);xlabel('X (km)','Fontsize',15);ylabel('Y (km)','Fontsize',15);zlabel('Time (s)','Fontsize',15);title('Ground-truth reflection','Fontsize',15,'fontweight','normal');set(gca,'Linewidth',2,'Fontsize',15);text(-0.5,-0.2, -0.3,'c)','color','k','Fontsize',30,'fontweight','bold','HorizontalAlignment','left');
subplot(3,2,5);drr_plot3d(diffr1,[100,20,20],z,x,y);caxis([-0.5,0.5]);xlabel('X (km)','Fontsize',15);ylabel('Y (km)','Fontsize',15);zlabel('Time (s)','Fontsize',15);title('LDRR diffraction','Fontsize',15,'fontweight','normal');set(gca,'Linewidth',2,'Fontsize',15);text(-0.5,-0.2, -0.3,'d)','color','k','Fontsize',30,'fontweight','bold','HorizontalAlignment','left');
subplot(3,2,6);drr_plot3d(data-diffr,[100,20,20],z,x,y);caxis([-0.5,0.5]);xlabel('X (km)','Fontsize',15);ylabel('Y (km)','Fontsize',15);zlabel('Time (s)','Fontsize',15);title('LDRR reflection','Fontsize',15,'fontweight','normal');set(gca,'Linewidth',2,'Fontsize',15);text(-0.5,-0.2, -0.3,'e)','color','k','Fontsize',30,'fontweight','bold','HorizontalAlignment','left');
print(gcf,'-dpng','-r300','test_matdrr_drr3d_diffraction.png');
print(gcf,'-depsc','-r300','test_matdrr_drr3d_diffraction.eps');

