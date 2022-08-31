% Demonstration script for 
% 2D denoising via the localized damped rank-reduction method
% (LDRR)
% 
%  Copyright (C) 2018 Yangkang Chen
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

%% download synthetic data
%please download the data from https://github.com/aaspip/data/blob/main/hevents.mat
load hevents.mat

d=d/max(max(d));
randn('state',201314);
dn=d+0.1*randn(size(d));

%% parameters
[n1,n2]=size(dn);
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
n1win=50;n2win=20;n3win=1;
r1=0.5;r2=0.5;r3=0.5;
param.N=2;
d2=drr3d_win(dn,param.flow,param.fhigh,param.dt,param.N,param.K,param.verb,n1win,n2win,n3win,r1,r2,r3);

%% Local DRR (Option II, separated version)
d3=drr_win3d(@localdrr3d, param, dn, n1win, n2win, n3win, r1, r2, r3);

%% SNR comparison
fprintf('SNR of DRR is %g\n',drr_snr(d,d1));
fprintf('SNR of LDRR is %g\n',drr_snr(d,d2));

%% verify correctness of two options (output is zero)
fprintf('Difference between two versions is %g\n',norm(d2-d3));

%% plot results
figure('units','normalized','Position',[0.2 0.4 0.55, 0.45]);
subplot(1,6,1);drr_imagesc(d);caxis([-0.5,0.5]);axis off;title('Clean');
subplot(1,6,2);drr_imagesc(dn);caxis([-0.5,0.5]);axis off;title('Noisy');
subplot(1,6,3);drr_imagesc(d1);caxis([-0.5,0.5]);axis off;title('DRR');
subplot(1,6,4);drr_imagesc(d2);caxis([-0.5,0.5]);axis off;title('LDRR');
subplot(1,6,5);drr_imagesc(dn-d1);caxis([-0.5,0.5]);axis off;title('DRR');
subplot(1,6,6);drr_imagesc(dn-d2);caxis([-0.5,0.5]);axis off;title('LDRR');
print(gcf,'-dpng','-r300','test_matdrr_drr2d_win.png');