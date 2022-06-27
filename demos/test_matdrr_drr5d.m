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
%% add noise
randn('state',201516);
dn=d+.2*randn(nt,nhx,nhy,nx,ny);

figure;
subplot(2,1,1);imagesc(reshape(d(:,:,:,1,1),100,10*10));colormap(jet);
subplot(2,1,2);imagesc(reshape(dn(:,:,:,1,1),100,10*10));colormap(jet);

%% denoise (Traditonal MSSA)
flow=5;fhigh=100;dt=0.004;N=4;
d1=drr5d(dn,flow,fhigh,dt,N,100,1);
figure;
subplot(3,1,1);imagesc(reshape(d(:,:,:,1,1),100,10*10));colormap(jet);
subplot(3,1,2);imagesc(reshape(dn(:,:,:,1,1),100,10*10));colormap(jet);
subplot(3,1,3);imagesc(reshape(d1(:,:,:,1,1),100,10*10));colormap(jet);

%% denoise (DMSSA)
flow=5;fhigh=100;dt=0.004;N=4;NN=2;
d2=drr5d(dn,flow,fhigh,dt,N,NN,1);
figure;
subplot(3,1,1);imagesc(reshape(d(:,:,:,1,1),100,10*10));colormap(jet);
subplot(3,1,2);imagesc(reshape(dn(:,:,:,1,1),100,10*10));colormap(jet);
subplot(3,1,3);imagesc(reshape(d2(:,:,:,1,1),100,10*10));colormap(jet);

s0=reshape(d(:,:,:,1,1),100,10*10);
sn=reshape(dn(:,:,:,1,1),100,10*10);
s1=reshape(d1(:,:,:,1,1),100,10*10);
s2=reshape(d2(:,:,:,1,1),100,10*10);

drr_snr(s0,sn) %-6.7209
drr_snr(s0,s1) %16.5590
drr_snr(s0,s2) %18.2507





