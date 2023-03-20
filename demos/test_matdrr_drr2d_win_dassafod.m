% DEMO script to generate the DAS example
% BY Yangkang Chen
% Jan, 2023
% This script takes about 10 minutes
%
% Dependency MATdrr
% svn co https://github.com/aaspip/MATdrr/trunk ./MATdrr
% or git clone https://github.com/aaspip/MATdrr ./
%
% Related pacakge:
% https://github.com/chenyk1990/dasmrrcoh
% https://github.com/chenyk1990/dasmrrcoh-dataonly

clc;clear;close all;
addpath(genpath('../matdrr'));

if ~isdir('fig')
    mkdir('fig');
end

if ~isdir('processed')
    mkdir('processed');
end

%% Download data first
%% The whole dataset in the folder "raw" can be downloaded from https://github.com/chenyk1990/dasmrrcoh-dataonly
%% A complete processing workflow can be found at https://github.com/chenyk1990/dasmrrcoh-dataonly or https://github.com/chenyk1990/dasmrrcoh
%
% https://github.com/chenyk1990/dasmrrcoh-dataonly/tree/main/raw/2017-06-23T22:03:22.380000Z_mag2.17.mat

name='2017-06-23T22:03:22.380000Z_mag2.17.mat'
ieq=3;
load(name);
eq=data;
d_bp=drr_bandpass(eq',1/250,0,20)';
d_bpmf=drr_mf(d_bp,5,1,1);
%% LDRR
n1win=1024;n2win=800;n3win=1;
n1win=512;n2win=200;n3win=1;
r1=0.5;r2=0.5;r3=0.5;
d_bpmfmrr=drr3d_win(d_bpmf',0,50,1/250,2,4,0,n1win,n2win,n3win,r1,r2,r3)';
save(sprintf('processed/eq%d.mat',ieq),'d_bp','d_bpmf','d_bpmfmrr');

[n1,n2]=size(data);
t=[0:n2-1]*(1/250);
x=1:n1;

figure('units','normalized','Position',[0.2 0.4 0.7, 0.6],'color','w');
ax1=subplot(2,2,1);
drr_imagesc(eq,95,1,t,x);colormap(ax1,seis);
name(end-3:end)=[];
title(name,'Interpreter', 'none','Fontsize',14,'fontweight','bold');
% xlabel('Time (s)','Fontsize',14,'fontweight','bold');
ylabel('Channel','Fontsize',14,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',14,'Fontweight','bold');
text(-5,-100,'a)','color','k','Fontsize',18,'fontweight','bold','HorizontalAlignment','center');

ax2=subplot(2,2,2);
drr_imagesc(d_bp,95,1,t,x);colormap(ax2,seis);
title('BP','Fontsize',14,'fontweight','bold');
% xlabel('Time (s)','Fontsize',14,'fontweight','bold');
% ylabel('Channel','Fontsize',14,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',14,'Fontweight','bold');
text(-5,-100,'b)','color','k','Fontsize',18,'fontweight','bold','HorizontalAlignment','center');

ax3=subplot(2,2,3);
drr_imagesc(d_bpmf,95,1,t,x);colormap(ax3,seis);title('BP+MF');
title('BP+MF','Fontsize',14,'fontweight','bold');
% xlabel('Time (s)','Fontsize',14,'fontweight','bold');
ylabel('Channel','Fontsize',14,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',14,'Fontweight','bold');
text(-5,-100,'c)','color','k','Fontsize',18,'fontweight','bold','HorizontalAlignment','center');
xlabel('Time (s)','Fontsize',14,'fontweight','bold');

ax4=subplot(2,2,4);
drr_imagesc(d_bpmfmrr,95,1,t,x);colormap(ax4,seis);title('BP+MF+MRR');
title('BP+MF+MRR','Fontsize',14,'fontweight','bold');
% xlabel('Time (s)','Fontsize',14,'fontweight','bold');
% ylabel('Channel','Fontsize',14,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',14,'Fontweight','bold');
text(-5,-100,'d)','color','k','Fontsize',18,'fontweight','bold','HorizontalAlignment','center');
xlabel('Time (s)','Fontsize',14,'fontweight','bold');
print(gcf,'-depsc','-r300','test_matdrr_drr2d_win_dassafod.eps');
print(gcf,'-dpng','-r300','test_matdrr_drr2d_win_dassafod.png');
