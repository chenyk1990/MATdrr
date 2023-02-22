% Demonstration script for
% 3D seismic denoising and reconstruction via the optimally damped rank-reduction method
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

%% without noise
dn=d;

%% decimate
[nt,nx,ny]=size(d);
ratio=0.5;
mask=drr_genmask(reshape(d,nt,nx*ny),ratio,'c',201415);
mask=reshape(mask,nt,nx,ny);
d0=dn.*mask;

%% reconstruct (without denoising)
flow=0;fhigh=125;dt=0.004;N=3;Niter=10;mode=0;verb=1;
d1=drr3drecon(d0,mask,flow,fhigh,dt,50,N,Niter,eps,verb,mode);

% 2D quick comparison (clean,noisy,observed,reconstructed using RR)
figure;drr_imagesc([d(:,:,9),d0(:,:,9),d1(:,:,9)]);caxis([-0.5,0.5]);

%% simultaneous denoising and reconstruction
% adding noise
randn('state',201314);
var=0.2;
dn=d+var*randn(size(d));
d0=dn.*mask;

%% using RR (when K is suffiently large, see derivations in the references)
flow=0;fhigh=250;dt=0.002;N=6;Niter=10;mode=1;verb=1;
a=(Niter-(1:Niter))/(Niter-1); %linearly decreasing
d1=drr3drecon(d0,mask,flow,fhigh,dt,N,100,Niter,eps,verb,mode,a);

% 2D quick comparison
% figure;drr_imagesc([d(:,:,9),d0(:,:,9),d1(:,:,9)]);caxis([-0.5,0.5]);

%% using DRR
flow=0;fhigh=250;dt=0.002;N=6;Niter=10;mode=1;verb=1;K=2;
a=(Niter-(1:Niter))/(Niter-1); %linearly decreasing
d2=drr3drecon(d0,mask,flow,fhigh,dt,N,K,Niter,eps,verb,mode,a);

% 2D quick comparison
% figure;drr_imagesc([d(:,:,9),d0(:,:,9),d2(:,:,9)]);caxis([-0.5,0.5]);

%% using ODRR
flow=0;fhigh=250;dt=0.002;N=6;Niter=10;mode=1;verb=1;K=2;O=1;
a=(Niter-(1:Niter))/(Niter-1); %linearly decreasing
d3=odrr3drecon(d0,mask,flow,fhigh,dt,N,K,O,Niter,eps,verb,mode,a);

% 2D quick comparison
% figure;drr_imagesc([d(:,:,9),d0(:,:,9),d3(:,:,9)]);caxis([-0.5,0.5]);

%% calculate Signal-to-noise Ratio (SNR)
snr0=drr_snr(d,d0,2) %observed data
snrn=drr_snr(d,dn,2) %noisy data
snr1=drr_snr(d,d1,2) %RR method
snr2=drr_snr(d,d2,2) %DRR method
snr3=drr_snr(d,d3,2) %ODRR method

%SNR results when N=3 (might be slightly different for different PC platforms)
%d0: -5.9853
%d1: 1.4503
%d2: 6.4816
%d3: 7.6997 

%SNR results when N=6 (when rank over-estimated) (might be slightly different for different PC platforms)
%d0: -5.9853
%d1: -0.9812
%d2: 4.9273
%d3: 7.4988

t=[0:300-1]*0.004;
figure('units','normalized','Position',[0.2 0.4 0.4, 1],'color','w');
subplot(6,1,1);drr_imagesc(d(:,:),0.1,2,1:400,t);ylabel('Time (s)','Fontsize',12,'fontweight','normal');set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','normal');title(strcat('Clean'),'Fontsize',15,'fontweight','normal');text(-40,-0.3,'a)','color','k','Fontsize',25,'fontweight','bold','HorizontalAlignment','left');
subplot(6,1,2);drr_imagesc(dn(:,:),0.1,2,1:400,t);ylabel('Time (s)','Fontsize',12,'fontweight','normal');set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','normal');title(strcat('Noisy (SNR=',num2str(snrn),' dB )'),'Fontsize',15,'fontweight','normal');text(-40,-0.3,'b)','color','k','Fontsize',25,'fontweight','bold','HorizontalAlignment','left');
subplot(6,1,3);drr_imagesc(d0(:,:),0.1,2,1:400,t);ylabel('Time (s)','Fontsize',12,'fontweight','normal');set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','normal');title(strcat('Incomplete (SNR=',num2str(snr0),' dB )'),'Fontsize',15,'fontweight','normal');text(-40,-0.3,'c)','color','k','Fontsize',25,'fontweight','bold','HorizontalAlignment','left');
subplot(6,1,4);drr_imagesc(d1(:,:),0.1,2,1:400,t);ylabel('Time (s)','Fontsize',12,'fontweight','normal');set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','normal');title(strcat('RR (SNR=',num2str(snr1),' dB )'),'Fontsize',15,'fontweight','normal');text(-40,-0.3,'d)','color','k','Fontsize',25,'fontweight','bold','HorizontalAlignment','left');
subplot(6,1,5);drr_imagesc(d2(:,:),0.1,2,1:400,t);ylabel('Time (s)','Fontsize',12,'fontweight','normal');set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','normal');title(strcat('DRR (SNR=',num2str(snr2),' dB )'),'Fontsize',15,'fontweight','normal');text(-40,-0.3,'e)','color','k','Fontsize',25,'fontweight','bold','HorizontalAlignment','left');
subplot(6,1,6);drr_imagesc(d3(:,:),0.1,2,1:400,t);xlabel('Trace','Fontsize',12,'fontweight','normal');ylabel('Time (s)','Fontsize',12,'fontweight','normal');set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','normal');title(strcat('ODRR (SNR=',num2str(snr3),' dB )'),'Fontsize',15,'fontweight','normal');text(-40,-0.3,'f)','color','k','Fontsize',25,'fontweight','bold','HorizontalAlignment','left');
print(gcf,'-dpng','-r300','test_matdrr_odrr3drecon.png');
print(gcf,'-depsc','-r200','test_matdrr_odrr3drecon.eps');


