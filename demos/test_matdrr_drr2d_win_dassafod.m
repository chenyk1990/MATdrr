% Script to plot Figure 2
% BY Yangkang Chen
% Jan, 2023
% This script takes about 10 minutes
%
% Dependency MATdrr
% svn co https://github.com/chenyk1990/MATdrr/trunk ./MATdrr
% or git clone https://github.com/chenyk1990/MATdrr ./

clc;clear;close all;
addpath(genpath('~/MATdrr'));
addpath(genpath('~/dasmrrcoh'));

if ~isdir('fig')
    mkdir('fig');
end

if ~isdir('processed')
    mkdir('processed');
end

names=dir('raw/*.mat');

ieq=3;
for ii=ieq:ieq
    load(strcat(names(ii).folder,'/',names(ii).name));
    name=names(ii).name;
    if ii==ieq
        [n1,n2]=size(data);
        t=[0:n2-1]*(1/250);
        x=1:n1;
    end

    if ii==12
        data(find(isnan(data)))=0;
    end
    eq=data;
    d_bp=das_bandpass(eq',1/250,0,20)';
    d_bpmf=das_mf(d_bp,5,1,1);
    %% LDRR
    n1win=1024;n2win=800;n3win=1;
    n1win=512;n2win=200;n3win=1;
    r1=0.5;r2=0.5;r3=0.5;
    d_bpmfmrr=drr3d_win(d_bpmf',0,50,1/250,2,4,0,n1win,n2win,n3win,r1,r2,r3)';

    save(sprintf('processed/eq%d.mat',ii),'d_bp','d_bpmf','d_bpmfmrr');

    fprintf('event %d/%d is done\n',ii,length(names));
    %     close gcf;
end


%% waveform
% figure('units','normalized','Position',[0.2 0.4 0.8, 1],'color','w');
% subplot(4,2,1);
% das_imagesc(eq,95,1,t,x);title(name,'Interpreter', 'none');
% subplot(4,2,2);
% das_imagesc(d_bp,95,1,t,x);title('BP');
% subplot(4,2,3);
% das_imagesc(d_bpmf,95,1,t,x);title('MF (5)');%xlim([7.5,8.5]);
% subplot(4,2,5);
% das_imagesc(d_bpmf,95,1,t,x);title('MF (5)');xlim([0,30]);ylim([0,200]);
% subplot(4,2,4);
% das_imagesc(d_bpmfmrr,95,1,t,x);title('LDRR');
% subplot(4,2,6);
% das_imagesc(d_bpmfmrr,95,1,t,x);title('LDRR');xlim([0,30]);ylim([0,200]);

%% coherency
% d_bp=data';
% d_bpmf=datat';
% d_bpmfmrr=datatt';


[nx,nt]=size(eq);
v=linspace(-0.0013,0.0013,100);
% v=linspace(-0.0002,0.0002,200);
Param.v=v;
Param.nt=nt;
Param.h=[0:nx-1];
Param.dt=1/250.0;
Param.type=1;
Param.oper=-1;
c0=das_coh(eq',Param);
c_bp=das_coh(d_bp',Param);
c_bpmf=das_coh(d_bpmf',Param);
c_bpmfmrr=das_coh(d_bpmfmrr',Param);



figure('units','normalized','Position',[0.2 0.4 0.7, 0.6],'color','w');
ax1=subplot(2,2,1);
das_imagesc(eq,95,1,t,x);colormap(ax1,seis);
title(name,'Interpreter', 'none','Fontsize',14,'fontweight','bold');
% xlabel('Time (s)','Fontsize',14,'fontweight','bold');
ylabel('Channel','Fontsize',14,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',14,'Fontweight','bold');
text(-5,-100,'a)','color','k','Fontsize',18,'fontweight','bold','HorizontalAlignment','center');



ax2=subplot(2,2,2);
das_imagesc(d_bp,95,1,t,x);colormap(ax2,seis);
title('BP','Fontsize',14,'fontweight','bold');
% xlabel('Time (s)','Fontsize',14,'fontweight','bold');
% ylabel('Channel','Fontsize',14,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',14,'Fontweight','bold');
text(-5,-100,'b)','color','k','Fontsize',18,'fontweight','bold','HorizontalAlignment','center');

ax3=subplot(2,2,3);
das_imagesc(d_bpmf,95,1,t,x);colormap(ax3,seis);title('BP+MF');
title('BP+MF','Fontsize',14,'fontweight','bold');
% xlabel('Time (s)','Fontsize',14,'fontweight','bold');
ylabel('Channel','Fontsize',14,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',14,'Fontweight','bold');
text(-5,-100,'c)','color','k','Fontsize',18,'fontweight','bold','HorizontalAlignment','center');
xlabel('Time (s)','Fontsize',14,'fontweight','bold');

ax4=subplot(2,2,4);
das_imagesc(d_bpmfmrr,95,1,t,x);colormap(ax4,seis);title('BP+MF+MRR');
title('BP+MF+MRR','Fontsize',14,'fontweight','bold');
% xlabel('Time (s)','Fontsize',14,'fontweight','bold');
% ylabel('Channel','Fontsize',14,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',14,'Fontweight','bold');
text(-5,-100,'d)','color','k','Fontsize',18,'fontweight','bold','HorizontalAlignment','center');
xlabel('Time (s)','Fontsize',14,'fontweight','bold');
print(gcf,'-depsc','-r300','test_matdrr_drr2d_win_dassafod.eps');
print(gcf,'-dpng','-r300','test_matdrr_drr2d_win_dassafod.png');
