%  5D seismic data denoising via the damped rank-reduction method
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
%% add noise
randn('state',201516);
dn=d+.2*randn(nt,nhx,nhy,nx,ny);

figure;
subplot(2,1,1);drr_imagesc(reshape(d(:,:,:,1,1),100,10*10));
subplot(2,1,2);drr_imagesc(reshape(dn(:,:,:,1,1),100,10*10));

%% denoise (Traditonal RR)
flow=1;fhigh=100;dt=0.004;N=20;
d1=drr5d(dn,flow,fhigh,dt,N,100,1);
figure;
subplot(3,1,1);drr_imagesc(reshape(d(:,:,:,1,1),100,10*10));
subplot(3,1,2);drr_imagesc(reshape(dn(:,:,:,1,1),100,10*10));
subplot(3,1,3);drr_imagesc(reshape(d1(:,:,:,1,1),100,10*10));

%% denoise (DRR)
flow=1;fhigh=100;dt=0.004;N=20;K=4;
d2=drr5d(dn,flow,fhigh,dt,N,K,1);
figure;
subplot(3,1,1);drr_imagesc(reshape(d(:,:,:,1,1),100,10*10));
subplot(3,1,2);drr_imagesc(reshape(dn(:,:,:,1,1),100,10*10));
subplot(3,1,3);drr_imagesc(reshape(d2(:,:,:,1,1),100,10*10));

%% denoise (ODRR)
flow=1;fhigh=100;dt=0.004;N=20;K=4;O=1;
d3=odrr5d(dn,flow,fhigh,dt,N,K,O,1);
figure;
subplot(3,1,1);drr_imagesc(reshape(d(:,:,:,1,1),100,10*10));
subplot(3,1,2);drr_imagesc(reshape(dn(:,:,:,1,1),100,10*10));
subplot(3,1,3);drr_imagesc(reshape(d3(:,:,:,1,1),100,10*10));

%% calculate Signal-to-noise Ratio (SNR)
drr_snr(d(:,:),dn(:,:))
drr_snr(d(:,:),d1(:,:))
drr_snr(d(:,:),d2(:,:))
drr_snr(d(:,:),d3(:,:))


% SNR when N=20; (when rank over-estimated)
% sn=-6.6662;
% s1=10.4512;
% s2=18.2299;
% s3=19.0101;





