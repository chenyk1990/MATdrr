%  5D seismic data denoising via the optimally damped rank-reduction method
%  
%  Copyright (C) 2022 The University of Texas at Austin
%  Copyright (C) 2022 Yangkang Chen
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
%  [1] Bai et al., 2020, Seismic signal enhancement based on the lowrank methods, Geophysical Prospecting, 68, 2783-2807.
%  [2] Chen et al., 2020, Five-dimensional seismic data reconstruction using the optimally damped rank-reduction method, Geophysical Journal International, 222, 1824-1845.
%  [3] Chen, Y., W. Huang, D. Zhang, W. Chen, 2016, An open-source matlab code package for improved rank-reduction 3D seismic data denoising and reconstruction, Computers & Geosciences, 95, 59-66.
%  [4] Chen, Y., D. Zhang, Z. Jin, X. Chen, S. Zu, W. Huang, and S. Gan, 2016, Simultaneous denoising and reconstruction of 5D seismic data via damped rank-reduction method, Geophysical Journal International, 206, 1695-1717.
%  [5] Huang, W., R. Wang, Y. Chen, H. Li, and S. Gan, 2016, Damped multichannel singular spectrum analysis for 3D random noise attenuation, Geophysics, 81, V261-V270.
%  [6] Chen et al., 2017, Preserving the discontinuities in least-squares reverse time migration of simultaneous-source data, Geophysics, 82, S185-S196.
%  [7] Chen et al., 2019, Obtaining free USArray data by multi-dimensional seismic reconstruction, Nature Communications, 10:4434.
%  [8] Oboue et al., 2021, Robust damped rank-reduction method for simultaneous denoising and reconstruction of 5-D seismic data, Geophysics, 86, V71–V89.

clc;clear;close all;
addpath(genpath('../matdrr'));

%% download and load data
% https://github.com/aaspip/data/blob/main/yc_synth5d.mat
load yc_synth5d.mat

d=data5d;d=d/max(max(max(max(max(d)))));
[nt,nhx,nhy,nx,ny]=size(d);
dt=0.004;

%% exploring the data
%1) ploting CMP gather
figure;drr_imagesc(reshape(d(:,:,:,1,1),100,10*10));

%2) ploting common offset gather
figure;drr_imagesc(reshape(d(:,5,5,:,:),100,10*10));

%% simultaneous denoising and reconstruction
randn('state',201314);
var=0.25;
dn=d+var*randn(size(d));

%% decimate
[nt,nhx,nhy,nx,ny]=size(d);
ratio=0.3;
mask=drr_genmask(reshape(d,nt,nhx*nhy*nx*ny),ratio,'c',201415);
mask=reshape(mask,nt,nhx,nhy,nx,ny);
d0=dn.*mask;

%% RR5D
flow=1;fhigh=100;dt=0.004;N=6;Niter=10;mode=1;verb=1;
a=(Niter-(1:Niter))/(Niter-1); %linearly decreasing
d1=drr5drecon(d0,mask,flow,fhigh,dt,N,50,Niter,eps,verb,mode,a);

figure;
subplot(3,1,1);drr_imagesc(reshape(d(:,:,:,1,1),nt,10*10));
subplot(3,1,2);drr_imagesc(reshape(d0(:,:,:,1,1),nt,10*10));
subplot(3,1,3);drr_imagesc(reshape(d1(:,:,:,1,1),nt,10*10));

figure;
subplot(3,1,1);drr_imagesc(reshape(d(:,5,5,:,:),nt,10*10));
subplot(3,1,2);drr_imagesc(reshape(d0(:,5,5,:,:),nt,10*10));
subplot(3,1,3);drr_imagesc(reshape(d1(:,5,5,:,:),nt,10*10));

%% DRR5D
flow=1;fhigh=100;dt=0.004;N=6;K=4;Niter=10;mode=1;verb=1;
d2=drr5drecon(d0,mask,flow,fhigh,dt,N,K,Niter,eps,verb,mode,a);

figure;
subplot(3,1,1);drr_imagesc(reshape(d(:,:,:,1,1),nt,10*10));
subplot(3,1,2);drr_imagesc(reshape(d0(:,:,:,1,1),nt,10*10));
subplot(3,1,3);drr_imagesc(reshape(d2(:,:,:,1,1),nt,10*10));

figure;
subplot(3,1,1);drr_imagesc(reshape(d(:,5,5,:,:),nt,10*10));caxis([-0.3,0.3]);
subplot(3,1,2);drr_imagesc(reshape(d0(:,5,5,:,:),nt,10*10));caxis([-0.3,0.3]);
subplot(3,1,3);drr_imagesc(reshape(d2(:,5,5,:,:),nt,10*10));caxis([-0.3,0.3]);

%% ODRR5D
flow=1;fhigh=100;dt=0.004;N=6;K=4;Niter=10;mode=1;verb=1;O=1;
d3=odrr5drecon(d0,mask,flow,fhigh,dt,N,K,O,Niter,eps,verb,mode,a);

figure;
subplot(3,1,1);drr_imagesc(reshape(d(:,:,:,1,1),nt,10*10));
subplot(3,1,2);drr_imagesc(reshape(d0(:,:,:,1,1),nt,10*10));
subplot(3,1,3);drr_imagesc(reshape(d3(:,:,:,1,1),nt,10*10));

figure;
subplot(3,1,1);drr_imagesc(reshape(d(:,5,5,:,:),nt,10*10));caxis([-0.3,0.3]);
subplot(3,1,2);drr_imagesc(reshape(d0(:,5,5,:,:),nt,10*10));caxis([-0.3,0.3]);
subplot(3,1,3);drr_imagesc(reshape(d3(:,5,5,:,:),nt,10*10));caxis([-0.3,0.3]);

%% calculate Signal-to-noise Ratio (SNR)
drr_snr(d(:,:),dn(:,:)) %-8.6178
drr_snr(d(:,:),d0(:,:)) %-4.5929
drr_snr(d(:,:),d1(:,:)) %7.4525
drr_snr(d(:,:),d2(:,:)) %10.5141
drr_snr(d(:,:),d3(:,:)) %11.8747


%% calculate Signal-to-noise Ratio (SNR)
snrnn=drr_snr(d(:,:),dn(:,:)) %
snr00=drr_snr(d(:,:),d0(:,:)) %
snr11=drr_snr(d(:,:),d1(:,:)) %
snr22=drr_snr(d(:,:),d2(:,:)) %
snr33=drr_snr(d(:,:),d3(:,:)) %

t=[0:100-1]*0.004;
%% visualization
figure('units','normalized','Position',[0.2 0.4 0.8, 1],'color','w');
subplot(6,2,1);drr_imagesc(reshape(d(:,:,:,1,1),100,10*10),0.1,2,1:100,t);ylabel('Time (s)','Fontsize',12,'fontweight','normal');set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','normal');title(strcat('Clean'),'Fontsize',15,'fontweight','normal');text(-9,-0.1,'a)','color','k','Fontsize',24,'fontweight','bold','HorizontalAlignment','left');
subplot(6,2,3);drr_imagesc(reshape(dn(:,:,:,1,1),100,10*10),0.1,2,1:100,t);ylabel('Time (s)','Fontsize',12,'fontweight','normal');set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','normal');title(strcat('Noisy (SNR=',num2str(snrnn),' dB )'),'Fontsize',15,'fontweight','normal');text(-9,-0.1,'c)','color','k','Fontsize',24,'fontweight','bold','HorizontalAlignment','left');
subplot(6,2,5);drr_imagesc(reshape(d0(:,:,:,1,1),100,10*10),0.1,2,1:100,t);ylabel('Time (s)','Fontsize',12,'fontweight','normal');set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','normal');title(strcat('Incomplete (SNR=',num2str(snr00),' dB )'),'Fontsize',15,'fontweight','normal');text(-9,-0.1,'e)','color','k','Fontsize',24,'fontweight','bold','HorizontalAlignment','left');
subplot(6,2,7);drr_imagesc(reshape(d1(:,:,:,1,1),100,10*10),0.1,2,1:100,t);ylabel('Time (s)','Fontsize',12,'fontweight','normal');set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','normal');title(strcat('RR (SNR=',num2str(snr11),' dB )'),'Fontsize',15,'fontweight','normal');text(-9,-0.1,'g)','color','k','Fontsize',24,'fontweight','bold','HorizontalAlignment','left');
subplot(6,2,9);drr_imagesc(reshape(d2(:,:,:,1,1),100,10*10),0.1,2,1:100,t);ylabel('Time (s)','Fontsize',12,'fontweight','normal');set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','normal');title(strcat('DRR (SNR=',num2str(snr22),' dB )'),'Fontsize',15,'fontweight','normal');text(-9,-0.1,'i)','color','k','Fontsize',24,'fontweight','bold','HorizontalAlignment','left');
subplot(6,2,11);drr_imagesc(reshape(d3(:,:,:,1,1),100,10*10),0.1,2,1:100,t);xlabel('Offset trace','Fontsize',12,'fontweight','normal');ylabel('Time (s)','Fontsize',12,'fontweight','normal');set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','normal');title(strcat('ODRR (SNR=',num2str(snr33),' dB )'),'Fontsize',15,'fontweight','normal');text(-9,-0.1,'k)','color','k','Fontsize',24,'fontweight','bold','HorizontalAlignment','left');

subplot(6,2,2);drr_imagesc(reshape(d(:,5,5,:,:),100,10*10),0.1,2,1:100,t);ylabel('Time (s)','Fontsize',12,'fontweight','normal');set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','normal');title(strcat('Clean'),'Fontsize',15,'fontweight','normal');text(-9,-0.1,'b)','color','k','Fontsize',24,'fontweight','bold','HorizontalAlignment','left');
subplot(6,2,4);drr_imagesc(reshape(dn(:,5,5,:,:),100,10*10),0.1,2,1:100,t);ylabel('Time (s)','Fontsize',12,'fontweight','normal');set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','normal');title(strcat('Noisy (SNR=',num2str(snrnn),' dB )'),'Fontsize',15,'fontweight','normal');text(-9,-0.1,'d)','color','k','Fontsize',24,'fontweight','bold','HorizontalAlignment','left');
subplot(6,2,6);drr_imagesc(reshape(d0(:,5,5,:,:),100,10*10),0.1,2,1:100,t);ylabel('Time (s)','Fontsize',12,'fontweight','normal');set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','normal');title(strcat('Incomplete (SNR=',num2str(snr00),' dB )'),'Fontsize',15,'fontweight','normal');text(-9,-0.1,'f)','color','k','Fontsize',24,'fontweight','bold','HorizontalAlignment','left');
subplot(6,2,8);drr_imagesc(reshape(d1(:,5,5,:,:),100,10*10),0.1,2,1:100,t);ylabel('Time (s)','Fontsize',12,'fontweight','normal');set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','normal');title(strcat('RR (SNR=',num2str(snr11),' dB )'),'Fontsize',15,'fontweight','normal');text(-9,-0.1,'h)','color','k','Fontsize',24,'fontweight','bold','HorizontalAlignment','left');
subplot(6,2,10);drr_imagesc(reshape(d2(:,5,5,:,:),100,10*10),0.1,2,1:100,t);ylabel('Time (s)','Fontsize',12,'fontweight','normal');set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','normal');title(strcat('DRR (SNR=',num2str(snr22),' dB )'),'Fontsize',15,'fontweight','normal');text(-9,-0.1,'j)','color','k','Fontsize',24,'fontweight','bold','HorizontalAlignment','left');
subplot(6,2,12);drr_imagesc(reshape(d3(:,5,5,:,:),100,10*10),0.1,2,1:100,t);xlabel('Midpoint trace','Fontsize',12,'fontweight','normal');ylabel('Time (s)','Fontsize',12,'fontweight','normal');set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','normal');title(strcat('ODRR (SNR=',num2str(snr33),' dB )'),'Fontsize',15,'fontweight','normal');text(-9,-0.1,'l)','color','k','Fontsize',24,'fontweight','bold','HorizontalAlignment','left');

print(gcf,'-dpng','-r300','test_matdrr_odrr5drecon.png');
print(gcf,'-depsc','-r200','test_matdrr_odrr5drecon.eps');


