%  test ssa in windows  
%
%  Copyright (C) 2018 Yangkang Chen
%  
%  This program is free software; you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation; either version 2 of the License, or
%  (at your option) any later version.
%  
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%  
%  You should have received a copy of the GNU General Public License
%  along with this program; if not, write to the Free Software
%  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  
%  
%  Reference:   Simultaneous seismic data denoising and reconstruction via multichannel
%               singular spectrum analysis, Geophysics, 2011, 76, V25-V32
%           
%               Chen et al., 2019, Obtaining free USArray data by multi-dimensional
%               seismic reconstruction, Nature Communications, 10:4434.
%

clc;clear;close all;
%%  create synthetic data
%d=hevents;
%please download the data from https://github.com/aaspip/data/blob/main/hevents.mat
load hevents.mat

d=d/max(max(d));
dn=d+0.1*randn(size(d));


[n1,n2]=size(dn);
param.dt=0.002;
param.lf=5;
param.hf=245;
param.N=2;
param.verb=1;


n1win=50;n2win=20;
% twin=n1;xwin=n2;
r1=0.5;r2=0.5;
%% Main program goes here !
dssa_win=win2d(@localssa, param, dn, n1win, n2win, r1, r2);

% dssa_win2=process_win(@localssa, param, dn, n1win, n2win);
% norm(dssa_win-dssa_win2)%for some reasons not the same
% figure;imagesc([d,dn,dssa_win,-(dn-dssa_win)]);colorbar;colormap(seis);

param.N=6;
dssa=fxmssa(dn,param.lf,param.hf,param.dt,param.N,param.verb);
figure;imagesc([d,dn,dssa_win,dssa,dn-dssa_win,dn-dssa]);caxis([-0.5,0.5]);colormap(seis);


%% option II
param.N=2;
dssa_win2=fxmssa_win(dn,param.lf,param.hf,param.dt,param.N,param.verb,n1win,n2win,r1,r2);
figure;imagesc([d,dn,dssa_win2,dssa,dn-dssa_win2,dn-dssa]);caxis([-0.5,0.5]);colormap(seis);


