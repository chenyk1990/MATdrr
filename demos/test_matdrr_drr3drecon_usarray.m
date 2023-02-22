%  Demo script for teleseismic denoising and reconstruction
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
addpath(genpath('../matdrr'));

%% Please download data from https://github.com/aaspip/data/blob/main/usarray_200901181411_wfm.mat
load('usarray_200901181411_wfm.mat');
% d, dists(shot receiver distance/offset in degree), stla, stlo, t 

%% rm bad trace
inds=[18,41,70];
d(:,inds)=[];d=drr_scale(d);
stlo(inds)=[];
stla(inds)=[];
d0=d(:,105:433);
stlo0=stlo(105:433);
stla0=stla(105:433);

%% 3D processing/reconstruction
mla=[33,49];
mlo=[-116,-102];
%binning
[d3d,x1,y1,mask]=drr_bin3d(d0,stlo0,stla0,16,28,mlo(1),mla(1),mlo(2),mla(2));
[stlo1,stla1]=meshgrid(x1,y1);
figure;plot(stlo0,stla0,'bv');hold on;
plot(stlo1(:),stla1(:),'r*');
print(gcf,'-depsc','-r200','grids.eps');

figure;imagesc(squeeze(mask(1,:,:))');colorbar;set(gca,'YDir','normal');

figure;
subplot(1,2,1);imagesc(squeeze(d3d(:,13,:)));caxis([-0.05,0.05]);colormap(gray);
subplot(1,2,2);imagesc(squeeze(d3d(:,:,10)));caxis([-0.05,0.05]);colormap(gray);


%% test if mask is correct
tt=(d3d.*mask-d3d);
norm(tt(:)) %0->correct

ratio=size(find(mask==1))/size(mask(:));
fprintf('Sampling ratio is %g\n',ratio);

figure;
subplot(1,2,1);imagesc(squeeze(mask(:,13,:)));caxis([0,1]);colormap(jet);colorbar;
subplot(1,2,2);imagesc(squeeze(mask(:,:,10)));caxis([0,1]);colormap(jet);colorbar;

%% global processing
flow=0;fhigh=0.5;dt=1;N=8;Niter=10;mode=1;verb=1;eps=0.00001;K=4;
a=(Niter-(1:Niter))/(Niter-1); %linearly decreasing
d1=drr3drecon(d3d,mask,flow,fhigh,dt,N,K,Niter,eps,verb,mode,a);

%% local processing (first way, directly utilize drr_win3dmask.m; second way, use drr3drecon_win.m, which is an integrated function, instead)
% [n1,n2,n3]=size(d3d);
% param.dt=dt;
% param.flow=flow;
% param.fhigh=fhigh;
% param.N=N;
% param.K=4;
% param.niter=Niter;
% param.a=(Niter-(1:Niter))/(Niter-1); %linearly decreasing;
% param.eps=0.00001;
% param.mode=1;
% param.verb=1;
% n1win=2000;n2win=n2;n3win=n3;
% r1=0.5;r2=0.5;r3=0.5;
% d2=drr_win3dmask(@localdrr3drecon, mask, param, d3d, n1win, n2win, n3win, r1, r2, r3);
% 
% param.amode=2;
% d3=drr_win3dmask(@localdrr3drecon_auto, mask, param, d3d, n1win, n2win, n3win, r1, r2, r3);

%% fixed-rank
[n1,n2,n3]=size(d3d);
n1win=2000;n2win=n2;n3win=n3;r1=0.5;r2=0.5;r3=0.5;
d2=drr3drecon_win(d3d,mask,flow,fhigh,dt,N,K,Niter,eps,verb,mode,a,n1win,n2win,n3win,r1,r2,r3);

%% automatically chosen ranks (cleaner)
amode=2;eps2=nan;%if amode=2, then eps3 is useless
%N=[4,N]; %N can also have a lower limit
d3=drr3drecon_win_auto(d3d,mask,flow,fhigh,dt,N,K,Niter,eps,verb,mode,amode,a,eps2,n1win,n2win,n3win,r1,r2,r3);

%% plot the results
ilon=12;
dtest=squeeze(d3d(:,ilon,:));
figure('units','normalized','Position',[0.2 0.4 1.0, 1.2],'color','w');
subplot(2,3,1);
drr_wigbh(dtest,stla1,t(1:5400),20);
ylim([31.5,50.2]);
ylabel('Latitude (^o)','Fontsize',20);
xlabel('Time (s)','Fontsize',20);
title('Raw','Fontsize',20);
text(-1400,53,'a)','color','k','Fontsize',30,'fontweight','bold','HorizontalAlignment','left');

subplot(2,3,2);
dtest=squeeze(d1(:,ilon,:));
drr_wigbh(dtest,stla1,t(1:5400),20);
ylim([31.5,50.2]);
ylabel('Latitude (^o)','Fontsize',20);
xlabel('Time (s)','Fontsize',20);
title('Global','Fontsize',20);
text(-1400,53,'b)','color','k','Fontsize',30,'fontweight','bold','HorizontalAlignment','left');

subplot(2,3,3);
dtest=squeeze(d2(:,ilon,:));
drr_wigbh(dtest,stla1,t(1:5400),20);
ylim([31.5,50.2]);
ylabel('Latitude (^o)','Fontsize',20);
xlabel('Time (s)','Fontsize',20);
title('Local','Fontsize',20);
text(-1400,53,'c)','color','k','Fontsize',30,'fontweight','bold','HorizontalAlignment','left');


%% zoomed comparison
subplot(2,3,4);
dtest=squeeze(d3d(1088:2587,ilon,8:12));
drr_wigbh(dtest,stla1(8:12,1),t(1088:2587),8);
ylim([37,39]);xlim([500,1400]);
ylabel('Latitude (^o)','Fontsize',20);
xlabel('Time (s)','Fontsize',20);
title('Raw','Fontsize',20);
text(350,39.2,'d)','color','k','Fontsize',30,'fontweight','bold','HorizontalAlignment','left');
% print(gcf,'-depsc','-r400','us_lon1_z.eps');  

subplot(2,3,5);
dtest1=squeeze(d1(1088:2587,ilon,8:12));
drr_wigbh(dtest1,stla1(8:12,1),t(1088:2587),8);
ylim([37,39]);xlim([500,1400]);
ylabel('Latitude (^o)','Fontsize',20);
xlabel('Time (s)','Fontsize',20);
title('Global','Fontsize',20);
text(350,39.2,'e)','color','k','Fontsize',30,'fontweight','bold','HorizontalAlignment','left');
% print(gcf,'-depsc','-r400','us_lon2_z.eps');  

subplot(2,3,6);
dtest2=squeeze(d2(1088:2587,ilon,8:12));
drr_wigbh(dtest2,stla1(8:12,1),t(1088:2587),8);
ylim([37,39]);xlim([500,1400]);
ylabel('Latitude (^o)','Fontsize',20);
xlabel('Time (s)','Fontsize',20);
text(350,39.2,'f)','color','k','Fontsize',30,'fontweight','bold','HorizontalAlignment','left');
title('Local','Fontsize',20);



%% following are all annotations
% Create textbox
annotation(gcf,'textbox',...
    [0.268055555555556 0.897013245033114 0.0491319444444445 0.0331125827814569],...
    'Color',[1 0 0],...
    'String',{'Rayleigh'},...
    'FontWeight','bold',...
    'FontSize',12,...
    'EdgeColor',[1 1 1]);

% Create textbox
annotation(gcf,'textbox',...
    [0.232638888888889 0.898337748344372 0.0258680555555556 0.0331125827814569],...
    'Color',[1 0 0],...
    'String',{'SS'},...
    'FontWeight','bold',...
    'FontSize',12,...
    'EdgeColor',[1 1 1]);

% Create textbox
annotation(gcf,'textbox',...
    [0.209722222222222 0.897013245033114 0.0203125 0.0331125827814569],...
    'Color',[1 0 0],...
    'String',{'S'},...
    'FontWeight','bold',...
    'FontSize',12,...
    'EdgeColor',[1 1 1]);

% Create textbox
annotation(gcf,'textbox',...
    [0.179166666666667 0.898337748344372 0.0206597222222222 0.0331125827814569],...
    'Color',[1 0 0],...
    'String',{'P'},...
    'FontWeight','bold',...
    'FontSize',12,...
    'EdgeColor',[1 1 1]);

% Create rectangle
annotation(gcf,'rectangle',...
    [0.169444444444444 0.682119205298013 0.0388888888888889 0.0450331125827814],...
    'Color',[0 1 0],...
    'LineWidth',3);

% Create textbox
annotation(gcf,'textbox',...
    [0.459722222222223 0.899662251655633 0.0206597222222222 0.0331125827814569],...
    'Color',[1 0 0],...
    'String',{'P'},...
    'FontWeight','bold',...
    'FontSize',12,...
    'EdgeColor',[1 1 1]);

% Create rectangle
annotation(gcf,'rectangle',...
    [0.450000000000001 0.683443708609272 0.0388888888888889 0.0450331125827814],...
    'Color',[0 1 0],...
    'LineWidth',3);

% Create textbox
annotation(gcf,'textbox',...
    [0.490277777777778 0.898337748344375 0.0203125 0.0331125827814569],...
    'Color',[1 0 0],...
    'String',{'S'},...
    'FontWeight','bold',...
    'FontSize',12,...
    'EdgeColor',[1 1 1]);

% Create textbox
annotation(gcf,'textbox',...
    [0.513194444444445 0.899662251655634 0.0258680555555556 0.0331125827814569],...
    'Color',[1 0 0],...
    'String',{'SS'},...
    'FontWeight','bold',...
    'FontSize',12,...
    'EdgeColor',[1 1 1]);

% Create textbox
annotation(gcf,'textbox',...
    [0.548611111111112 0.898337748344376 0.0491319444444445 0.0331125827814569],...
    'Color',[1 0 0],...
    'String',{'Rayleigh'},...
    'FontWeight','bold',...
    'FontSize',12,...
    'EdgeColor',[1 1 1]);

% Create textbox
annotation(gcf,'textbox',...
    [0.741666666666667 0.899662251655631 0.0206597222222222 0.0331125827814569],...
    'Color',[1 0 0],...
    'String',{'P'},...
    'FontWeight','bold',...
    'FontSize',12,...
    'EdgeColor',[1 1 1]);

% Create rectangle
annotation(gcf,'rectangle',...
    [0.731944444444445 0.683443708609272 0.0388888888888889 0.0450331125827814],...
    'Color',[0 1 0],...
    'LineWidth',3);

% Create textbox
annotation(gcf,'textbox',...
    [0.772222222222223 0.898337748344374 0.0203125 0.0331125827814569],...
    'Color',[1 0 0],...
    'String',{'S'},...
    'FontWeight','bold',...
    'FontSize',12,...
    'EdgeColor',[1 1 1]);

% Create textbox
annotation(gcf,'textbox',...
    [0.795138888888889 0.899662251655633 0.0258680555555556 0.0331125827814569],...
    'Color',[1 0 0],...
    'String',{'SS'},...
    'FontWeight','bold',...
    'FontSize',12,...
    'EdgeColor',[1 1 1]);

% Create textbox
annotation(gcf,'textbox',...
    [0.830555555555556 0.898337748344374 0.0491319444444445 0.0331125827814569],...
    'Color',[1 0 0],...
    'String',{'Rayleigh'},...
    'FontWeight','bold',...
    'FontSize',12,...
    'EdgeColor',[1 1 1]);

% Create line
annotation(gcf,'line',[0.184027777777778 0.182638888888889],...
    [0.905960264900662 0.582781456953642],'Color',[1 0 0],'LineWidth',3,...
    'LineStyle','--');

% Create line
annotation(gcf,'line',[0.215277777777778 0.2125],...
    [0.901986754966887 0.582781456953643],'Color',[1 0 0],'LineWidth',3,...
    'LineStyle','--');

% Create line
annotation(gcf,'line',[0.2625 0.249305555555556],...
    [0.901986754966887 0.582781456953642],'Color',[1 0 0],'LineWidth',3,...
    'LineStyle','--');

% Create line
annotation(gcf,'line',[0.238888888888889 0.231944444444445],...
    [0.894039735099338 0.584105960264901],'Color',[1 0 0],'LineWidth',3,...
    'LineStyle','--');

% Create line
annotation(gcf,'line',[0.465277777777778 0.463888888888889],...
    [0.905960264900662 0.582781456953642],'Color',[1 0 0],'LineWidth',3,...
    'LineStyle','--');

% Create line
annotation(gcf,'line',[0.496527777777778 0.493750000000001],...
    [0.900662251655629 0.581456953642384],'Color',[1 0 0],'LineWidth',3,...
    'LineStyle','--');

% Create line
annotation(gcf,'line',[0.520138888888889 0.513194444444445],...
    [0.894039735099338 0.584105960264901],'Color',[1 0 0],'LineWidth',3,...
    'LineStyle','--');

% Create line
annotation(gcf,'line',[0.543750000000001 0.530555555555557],...
    [0.901986754966887 0.582781456953642],'Color',[1 0 0],'LineWidth',3,...
    'LineStyle','--');

% Create line
annotation(gcf,'line',[0.746527777777778 0.745138888888889],...
    [0.905960264900662 0.582781456953642],'Color',[1 0 0],'LineWidth',3,...
    'LineStyle','--');

% Create line
annotation(gcf,'line',[0.777777777777778 0.775],...
    [0.903311258278146 0.5841059602649],'Color',[1 0 0],'LineWidth',3,...
    'LineStyle','--');

% Create line
annotation(gcf,'line',[0.80138888888889 0.794444444444445],...
    [0.895364238410596 0.585430463576159],'Color',[1 0 0],'LineWidth',3,...
    'LineStyle','--');

% Create line
annotation(gcf,'line',[0.825 0.811805555555556],...
    [0.903311258278146 0.5841059602649],'Color',[1 0 0],'LineWidth',3,...
    'LineStyle','--');

print(gcf,'-dpng','-r300','test_matdrr_drr3drecon_usarray.png');
print(gcf,'-depsc','-r200','test_matdrr_drr3drecon_usarray.eps');


