%  Demo script for OTG teleseismic denoising and reconstruction
%  as introduced in Chen et al., 2019, Obtaining free USArray data by multi-dimensional seismic reconstruction, Nature Communications, 10:4434.
%  Another version of this example is at 
%  https://github.com/chenyk1990/reproducible_research/tree/master/nonMada/usarray
% 
%  This takes about 10 minutes
% 
%  Written by Yangkang Chen
%  Feb, 2018
%  Modified on Dec, 2020
%  Further polished on July, 2022
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
% REFERENCES
%  Chen et al., 2019, Obtaining free USArray data by multi-dimensional seismic reconstruction, Nature Communications, 10:4434.
%  Wang et al., 2020, Separation and imaging of seismic diffractions using a localized rank-reduction method with adaptively selected ranks, 85, V497?V506.
%  Chen et al., 2017, Preserving the discontinuities in least-squares reverse time migration of simultaneous-source data, Geophysics, 82, S185-S196.
%  Bai et al., 2020, Seismic signal enhancement based on the lowrank methods, Geophysical Prospecting, 68, 2783-2807.
%  Chen et al., 2020, Five-dimensional seismic data reconstruction using the optimally damped rank-reduction method, Geophysical Journal International, 222, 1824-1845.
%  Chen, Y., W. Huang, D. Zhang, W. Chen, 2016, An open-source matlab code package for improved rank-reduction 3D seismic data denoising and reconstruction, Computers & Geosciences, 95, 59-66.
%  Chen, Y., D. Zhang, Z. Jin, X. Chen, S. Zu, W. Huang, and S. Gan, 2016, Simultaneous denoising and reconstruction of 5D seismic data via damped rank-reduction method, Geophysical Journal International, 206, 1695-1717.
%  Huang, W., R. Wang, Y. Chen, H. Li, and S. Gan, 2016, Damped multichannel singular spectrum analysis for 3D random noise attenuation, Geophysics, 81, V261-V270.

clc;clear;close all;
addpath(genpath('~/MATdrr'));

%% Please download data from https://github.com/aaspip/data/blob/main/usarray_200901181411_wfm.mat
load('usarray_200901181411_wfm.mat');
% d, dists(shot receiver distance/offset in degree), stla, stlo, t 

% figure;yc_imagesc([d(:,:)]);

%% rm bad trace
inds=[18,41,70];
d(:,inds)=[];d=drr_scale(d);
stlo(inds)=[];
stla(inds)=[];
d0=d(:,105:433);
stlo0=stlo(105:433);
stla0=stla(105:433);



d0=d0(1001:3000,:);

%% 3D processing/reconstruction
mla=[33,49];
mlo=[-116,-102];
%binning
[d3d,x1,y1,mask]=drr_bin3d(d0,stlo0,stla0,16,28,mlo(1),mla(1),mlo(2),mla(2));
[stlo1,stla1]=meshgrid(x1,y1);
% figure;plot(stlo0,stla0,'bv');hold on;
% plot(stlo1(:),stla1(:),'r*');
% print(gcf,'-depsc','-r200','grids.eps');

% figure;imagesc(squeeze(mask(1,:,:))');colorbar;set(gca,'YDir','normal');
% 
% figure;
% subplot(1,2,1);imagesc(squeeze(d3d(:,13,:)));caxis([-0.05,0.05]);colormap(gray);
% subplot(1,2,2);imagesc(squeeze(d3d(:,:,10)));caxis([-0.05,0.05]);colormap(gray);
% 

%% test if mask is correct
% tt=(d3d.*mask-d3d);
% norm(tt(:)) %0->correct
% 
% ratio=size(find(mask==1))/size(mask(:));
% fprintf('Sampling ratio is %g\n',ratio);
% 
% figure;
% subplot(1,2,1);imagesc(squeeze(mask(:,13,:)));caxis([0,1]);colormap(jet);colorbar;
% subplot(1,2,2);imagesc(squeeze(mask(:,:,10)));caxis([0,1]);colormap(jet);colorbar;

%% global processing
flow=0;fhigh=0.5;dt=1;N=8;Niter=10;mode=1;verb=1;eps=0.00001;K=4;
a=(Niter-(1:Niter))/(Niter-1); %linearly decreasing
d1=drr3drecon(d3d,mask,flow,fhigh,dt,N,K,Niter,eps,verb,mode,a);


figure;yc_imagesc([d3d(:,:),d1(:,:)]);

d2=drr3drecon_otg(d0,stlo0,stla0,16,28,mlo(1),mla(1),mlo(2),mla(2),flow,fhigh,dt,N,K,Niter,eps,verb,mode);

figure;yc_imagesc([d3d(:,:),d1(:,:),d2(:,:)]);

% %% local processing (first way, directly utilize drr_win3dmask.m; second way, use drr3drecon_win.m, which is an integrated function, instead)
% % [n1,n2,n3]=size(d3d);
% % param.dt=dt;
% % param.flow=flow;
% % param.fhigh=fhigh;
% % param.N=N;
% % param.K=4;
% % param.niter=Niter;
% % param.a=(Niter-(1:Niter))/(Niter-1); %linearly decreasing;
% % param.eps=0.00001;
% % param.mode=1;
% % param.verb=1;
% % n1win=2000;n2win=n2;n3win=n3;
% % r1=0.5;r2=0.5;r3=0.5;
% % d2=drr_win3dmask(@localdrr3drecon, mask, param, d3d, n1win, n2win, n3win, r1, r2, r3);
% % 
% % param.amode=2;
% % d3=drr_win3dmask(@localdrr3drecon_auto, mask, param, d3d, n1win, n2win, n3win, r1, r2, r3);
% 
% %% fixed-rank
% [n1,n2,n3]=size(d3d);
% n1win=2000;n2win=n2;n3win=n3;r1=0.5;r2=0.5;r3=0.5;
% d2=drr3drecon_win(d3d,mask,flow,fhigh,dt,N,K,Niter,eps,verb,mode,a,n1win,n2win,n3win,r1,r2,r3);
% 
% %% automatically chosen ranks (cleaner)
% amode=2;eps2=nan;%if amode=2, then eps3 is useless
% %N=[4,N]; %N can also have a lower limit
% d3=drr3drecon_win_auto(d3d,mask,flow,fhigh,dt,N,K,Niter,eps,verb,mode,amode,a,eps2,n1win,n2win,n3win,r1,r2,r3);
% 
% 

%% plot the results
x1=750;
x2=1200;
y1=45;
y2=49.5;

ilon=12;
dtest=squeeze(d3d(:,ilon,:));
figure('units','normalized','Position',[0.2 0.4 1.0, 1.2],'color','w');
subplot(2,3,1);
drr_wigbh(dtest,stla1,t(1001:3000),20);
ylim([31.5,50.2]);
ylabel('Latitude (^o)','Fontsize',20);
xlabel('Time (s)','Fontsize',20);
title('Raw','Fontsize',20);
yc_framebox(x1,x2,y1,y2,'r',2);
text(0,52,'a)','color','k','Fontsize',30,'fontweight','bold','HorizontalAlignment','left');

subplot(2,3,2);
dtest=squeeze(d1(:,ilon,:));
drr_wigbh(dtest,stla1,t(1001:3000),20);
ylim([31.5,50.2]);
ylabel('Latitude (^o)','Fontsize',20);
xlabel('Time (s)','Fontsize',20);
title('DRR','Fontsize',20);
yc_framebox(x1,x2,y1,y2,'r',2);
text(0,52,'b)','color','k','Fontsize',30,'fontweight','bold','HorizontalAlignment','left');

subplot(2,3,3);
dtest=squeeze(d2(:,ilon,:));
drr_wigbh(dtest,stla1,t(1001:3000),20);
ylim([31.5,50.2]);
ylabel('Latitude (^o)','Fontsize',20);
xlabel('Time (s)','Fontsize',20);
title('OTG','Fontsize',20);
yc_framebox(x1,x2,y1,y2,'r',2);
text(0,52,'c)','color','k','Fontsize',30,'fontweight','bold','HorizontalAlignment','left');

%% zoomed comparison
subplot(2,3,4);
dtest=squeeze(d3d(:,ilon,:));
drr_wigbh(dtest,stla1,t(1001:3000),20);
ylim([y1,y2]);xlim([x1,x2]);
ylabel('Latitude (^o)','Fontsize',20);
xlabel('Time (s)','Fontsize',20);
title('Raw','Fontsize',20);
text(660,50,'d)','color','k','Fontsize',30,'fontweight','bold','HorizontalAlignment','left');

% print(gcf,'-depsc','-r400','us_lon1_z.eps');  

subplot(2,3,5);
dtest=squeeze(d1(:,ilon,:));
drr_wigbh(dtest,stla1,t(1001:3000),20);
ylim([y1,y2]);xlim([x1,x2]);
ylabel('Latitude (^o)','Fontsize',20);
xlabel('Time (s)','Fontsize',20);
title('DRR','Fontsize',20);
text(660,50,'e)','color','k','Fontsize',30,'fontweight','bold','HorizontalAlignment','left');

subplot(2,3,6);
dtest=squeeze(d2(:,ilon,:));
drr_wigbh(dtest,stla1,t(1001:3000),20);
ylim([y1,y2]);xlim([x1,x2]);
ylabel('Latitude (^o)','Fontsize',20);
xlabel('Time (s)','Fontsize',20);
title('OTG','Fontsize',20);
text(660,50,'f)','color','k','Fontsize',30,'fontweight','bold','HorizontalAlignment','left');


%% following are all annotations
% Create textbox
print(gcf,'-dpng','-r300','test_matdrr_drr3drecon_otg_usarray.png');
print(gcf,'-depsc','-r300','test_matdrr_drr3drecon_otg_usarray.eps');


