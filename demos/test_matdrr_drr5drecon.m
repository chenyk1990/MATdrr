%  5D seismic data denoising via rank reduction (DMSSA)
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
figure;imagesc(reshape(d(:,:,:,1,1),100,10*10));colormap(jet);

%2) ploting common offset gather
figure;imagesc(reshape(d(:,5,5,:,:),100,10*10));colormap(jet);

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
flow=5;fhigh=100;dt=0.004;N=6;Niter=10;mode=1;verb=1;iflb=0;
a=(Niter-(1:Niter))/(Niter-1); %linearly decreasing
d1=drr5drecon(d0,mask,flow,fhigh,dt,N,50,Niter,eps,verb,mode,iflb,a);

figure;
subplot(3,1,1);imagesc(reshape(d(:,:,:,1,1),nt,10*10));
subplot(3,1,2);imagesc(reshape(d0(:,:,:,1,1),nt,10*10));
subplot(3,1,3);imagesc(reshape(d1(:,:,:,1,1),nt,10*10));

figure;
subplot(3,1,1);imagesc(reshape(d(:,5,5,:,:),nt,10*10));
subplot(3,1,2);imagesc(reshape(d0(:,5,5,:,:),nt,10*10));
subplot(3,1,3);imagesc(reshape(d1(:,5,5,:,:),nt,10*10));

%% DRR5D
flow=5;fhigh=100;dt=0.004;N=6;K=2;Niter=10;mode=1;verb=1;iflb=0;
d2=drr5drecon(d0,mask,flow,fhigh,dt,N,K,Niter,eps,verb,mode,iflb,a);

figure;
subplot(3,1,1);imagesc(reshape(d(:,:,:,1,1),nt,10*10));
subplot(3,1,2);imagesc(reshape(d0(:,:,:,1,1),nt,10*10));
subplot(3,1,3);imagesc(reshape(d2(:,:,:,1,1),nt,10*10));

figure;
subplot(3,1,1);imagesc(reshape(d(:,5,5,:,:),nt,10*10));caxis([-0.3,0.3]);colormap(jet);
subplot(3,1,2);imagesc(reshape(d0(:,5,5,:,:),nt,10*10));caxis([-0.3,0.3]);colormap(jet);
subplot(3,1,3);imagesc(reshape(d2(:,5,5,:,:),nt,10*10));caxis([-0.3,0.3]);colormap(jet);


drr_snr(d(:,:),dn(:,:)) %-8.6178
drr_snr(d(:,:),d0(:,:)) %-4.5929
drr_snr(d(:,:),d1(:,:)) %7.5638
drr_snr(d(:,:),d2(:,:)) %11.5626




