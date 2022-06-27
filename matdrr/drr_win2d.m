function [ dout ] =win2d(oper, param, din, n1win, n2win, r1, r2)
% Processing in 2D windows
% 2D version similar to process_win
% coding strategy follows exactly usr/cm/Mfxmssa_win.c
% 
% din:          input data
% oper:         operator
% param:        parameters of operator
% n1win:        first window length
% n2win:        second window length
% r1:           first overlapping ratio
% r2:           second overlapping ratio
% 
% dout:     output data
%
% Author      : Yangkang Chen
%
% Date        : Jan, 2018
%
%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as pwublished
%  by the Free Software Foundation, either version 3 of the License, or
%  any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details: http://www.gnu.org/licenses/
%
% see also 
% win3d.m,win3d_mask.m
% 
% Example: 
% test/test_win2d_fxmssa.m
% test/test_win3d_fxymssa.m
% 

[n1,n2]=size(din);

if nargin==3
   n1win=n1;
   n2win=n2;
   r1=0.5;
   r2=0.5;
end

if nargin==5
   r1=0.5;
   r2=0.5;
end

nov1=(1-r1)*n1win;  % non-overlapping size 1
nov2=(1-r2)*n2win;  % non-overlapping size 2
ov1=r1*n1win;       % overlapping size 1
ov2=r2*n2win;       % overlapping size 2

n1pad=n1win;        %padding size 1
nw1=1;
while n1pad<n1
    n1pad=n1pad+nov1;nw1=nw1+1;
end

n2pad=n2win;        %padding size 2
nw2=1;
while n2pad<n2
    n2pad=n2pad+nov2;nw2=nw2+1;
end

D1=zeros(n1pad,n2pad);D1(1:n1,1:n2)=din;%copy din into padded D1
D2=zeros(n1pad,n2pad);

for iw2=0:nw2-1
for iw1=0:nw1-1
    s1=iw1*nov1;s2=iw2*nov2;
    dtmp=D1(s1+1:s1+n1win,s2+1:s2+n2win);

    %uncomment this line for checking the correctness (checked 100% correct)
    dtmp = feval(oper,dtmp,param);
    %only valid for space-independent param
    %for reconstruction, with mask, param needs to be changed

    dtmp=win_weight2d(dtmp,iw1,iw2,nw1,nw2,n1win,n2win,ov1,ov2);
    
    D2(s1+1:s1+n1win,s2+1:s2+n2win)=D2(s1+1:s1+n1win,s2+1:s2+n2win)+dtmp;
end
end

dout=D2(1:n1,1:n2);

return


function [dout ]= win_weight2d(din,iw1,iw2,nw1,nw2,n1win,n2win,ov1,ov2)
%       follow exactly usr/cm/win.c
%       float din /*input data*/,
% 		int iw1 /*starting window 1 in dst*/,
% 		int iw2 /*starting window 2 in dst*/,
% 		int nw1 /*no of windows 1 in src*/,
% 		int nw2 /*no of windows 2 in src*/,
% 		int n1win /*window length 1 in src*/,
% 		int n2win /*window legnth 2 in src*/,
% 		int ov1 /*copy length in axis1*/,
% 		int ov2 /*copy length in axis2*/)

if iw2~=0
    for i1=0:n1win-1
        for i2=0:ov2-1
            din(i1+1,i2+1)=din(i1+1,i2+1)*(i2+1)/(ov2+1);
        end
    end
    
end

if iw2~=nw2-1
    
    for i1=0:n1win-1
        for i2=0:ov2-1
            din(i1+1,n2win-ov2+i2+1) = din(i1+1,n2win-ov2+i2+1)*(ov2-i2)/(ov2+1);
        end
    end
    
end

if iw1~=0
    for i2=0:n2win-1
        for i1=0:ov1-1
            din(i1+1,i2+1)=din(i1+1,i2+1)*(i1+1)/(ov1+1);
        end
    end
    
end

if iw1~=nw1-1
    for i2=0:n2win-1
        for i1=0:ov1-1
            din(n1win-ov1+i1+1,i2+1)=din(n1win-ov1+i1+1,i2+1)*(ov1-i1)/(ov1+1);
        end
    end
end

dout=din;
return












