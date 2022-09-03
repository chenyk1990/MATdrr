function [ D1 ] = odrr5d(D,flow,fhigh,dt,N,K,O,verb)
%  DRR5D: Optimally damped rank-reduction method for 5D denoising
%
%  IN   D:   	intput 5D data
%       flow:   processing frequency range (lower)
%       fhigh:  processing frequency range (higher)
%       dt:     temporal sampling interval
%       N:      number of singular value to be preserved
%       K:      damping factor (default: 4)
%       O:	    ODRR flag; 0 for DRR and 1 for ODRR (default: 1)
%
%       verb:   verbosity flag (default: 0)
%
%  OUT  D1:  	output data
%
%  Copyright (C) 2015 The University of Texas at Austin
%  Copyright (C) 2015 Yangkang Chen
%  Further modified by Yangkang Chen and Min Bai 2018-2020
%  Finalized by Yangkang Chen 2022
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

if nargin==0
    error('Input data must be provided!');
end

if nargin==1
    flow=1;
    fhigh=124;
    dt=0.004;
    N=1;
    K=4;
    O=1;
    verb=0;
end;


[nt,nx,ny,nhx,nhy]=size(D);
D1=zeros(nt,nx,ny,nhx,nhy);

nf=2^nextpow2(nt);

% Transform into F-X domain
DATA_FX=fft(D,nf,1);
DATA_FX0=zeros(nf,nx,ny,nhx,nhy);

% First and last nts of the DFT.
ilow  = floor(flow*dt*nf)+1;

if ilow<1;
    ilow=1;
end;

ihigh = floor(fhigh*dt*nf)+1;

if ihigh > floor(nf/2)+1;
    ihigh=floor(nf/2)+1;
end

lx=floor(nx/2)+1;
lxx=nx-lx+1;
ly=floor(ny/2)+1;
lyy=ny-ly+1;
lhx=floor(nhx/2)+1;
lhxx=nhx-lhx+1;
lhy=floor(nhy/2)+1;
lhyy=nhy-lhy+1;
M=zeros(lx*ly*lhx*lhy,lxx*lyy*lhxx,lhyy);

% main loop
for k=ilow:ihigh
    if(ny==1 & nhx==1 & nhy==1)
        M=P_H(DATA_FX(k,:,:,:,:).',lx,ly,lhx,lhy);
    else
        M=P_H(squeeze(DATA_FX(k,:,:,:,:)),lx,ly,lhx,lhy);
    end
    
    if O==0
        M=P_R(M,N,K);
    else
        M=P_OR(M,N,K);
    end
    
    DATA_FX0(k,:,:,:,:)=P_A(M,nx,ny,nhx,nhy,lx,ly,lhx,lhy);
    
    if(mod(k,5)==0 && verb==1)
        fprintf( 'F %d is done!\n\n',k);
    end
end

% Honor symmetries
for k=nf/2+2:nf
    DATA_FX0(k,:,:,:,:) = conj(DATA_FX0(nf-k+2,:,:,:,:));
end

% Back to TX (the output)
D1=real(ifft(DATA_FX0,[],1));
D1=D1(1:nt,:,:,:,:);

return


function [dout]=P_H(din,lx,ly,lhx,lhy)
% forming block Hankel matrix (5D version)
[nx,ny,nhx,nhy]=size(din);
lxx=nx-lx+1;
lyy=ny-ly+1;
lhxx=nhx-lhx+1;
lhyy=nhy-lhy+1;

for ky=1:nhy
    for kx=1:nhx
        for j=1:ny
            r1=hankel(din(1:lx,j,kx,ky),[din(lx:nx,j,kx,ky)]);
            if j<ly
                for id=1:j
                    r1o(1+(j-1)*lx-(id-1)*lx:j*lx-(id-1)*lx,1+(id-1)*lxx:lxx+(id-1)*lxx) = r1;
                end
            else
                for id=1:(ny-j+1)
                    r1o((ly-1)*lx+1-(id-1)*lx:ly*lx-(id-1)*lx,(j-ly)*lxx+1+(id-1)*lxx:(j-ly+1)*lxx+(id-1)*lxx)=r1;
                end
            end
        end
        r2=r1o;r1o=[];
        if kx<lhx
            for id=1:kx
                r2o(1+(kx-1)*lx*ly-(id-1)*lx*ly:kx*lx*ly-(id-1)*lx*ly,1+(id-1)*lxx*lyy:lxx*lyy+(id-1)*lxx*lyy) = r2;
            end
        else
            for id=1:(nhx-kx+1)
                r2o((lhx-1)*lx*ly+1-(id-1)*lx*ly:lhx*lx*ly-(id-1)*lx*ly,(kx-lhx)*lxx*lyy+1+(id-1)*lxx*lyy:(kx-lhx+1)*lxx*lyy+(id-1)*lxx*lyy)=r2;
            end
        end
    end
    r3=r2o;r2o=[];
    if ky<lhy
        for id=1:ky
            r3o(1+(ky-1)*lx*ly*lhx-(id-1)*lx*ly*lhx:ky*lx*ly*lhx-(id-1)*lx*ly*lhx,1+(id-1)*lxx*lyy*lhxx:lxx*lyy*lhxx+(id-1)*lxx*lyy*lhxx) = r3;
        end
    else
        for id=1:(nhy-ky+1)
            r3o((lhy-1)*lx*ly*lhx+1-(id-1)*lx*ly*lhx:lhy*lx*ly*lhx-(id-1)*lx*ly*lhx,(ky-lhy)*lxx*lyy*lhxx+1+(id-1)*lxx*lyy*lhxx:(ky-lhy+1)*lxx*lyy*lhxx+(id-1)*lxx*lyy*lhxx)=r3;
        end
    end
end
dout=r3o;
return

function [dout]=P_R(din,N,K)
% Rank reduction on the block Hankel matrix


%      [U,D,V]=svds(din,N); % a little bit slower for small matrix
%      dout=U*D*V';
% %
[U,D,V]=svds(din,N+1);
for j=1:N
    D(j,j)=D(j,j)*(1-D(N+1,N+1)^K/(D(j,j)^K+0.000000000000001));
end

dout=U(:,1:N)*D(1:N,1:N)*(V(:,1:N)');

return

function [dout]=P_OR(din,N,K)
% Rank reduction on the block Hankel matrix

dout=drr_optshrink_damp(din,N,K);

return

function [dout]=P_A(din,nx,ny,nhx,nhy,lx,ly,lhx,lhy)
% Averaging the block Hankel matrix to output the result (5D version)
lxx=nx-lx+1;
lyy=ny-ly+1;
lhxx=nhx-lhx+1;
lhyy=nhy-lhy+1;
dout=zeros(nx,ny,nhx,nhy);


for ky=1:nhy
    r3o=zeros(lx*ly*lhx,lxx*lyy*lhxx);
    if ky<lhy
        for id=1:ky
            r3o=r3o+din(1+(ky-1)*lx*ly*lhx-(id-1)*lx*ly*lhx:ky*lx*ly*lhx-(id-1)*lx*ly*lhx,1+(id-1)*lxx*lyy*lhxx:lxx*lyy*lhxx+(id-1)*lxx*lyy*lhxx)/ky;
        end
    else
        for id=1:(nhy-ky+1)
            r3o=r3o+din((lhy-1)*lx*ly*lhx+1-(id-1)*lx*ly*lhx:lhy*lx*ly*lhx-(id-1)*lx*ly*lhx,(ky-lhy)*lxx*lyy*lhxx+1+(id-1)*lxx*lyy*lhxx:(ky-lhy+1)*lxx*lyy*lhxx+(id-1)*lxx*lyy*lhxx)/(nhy-ky+1);
        end
    end
    for kx=1:nhx
        r2o=zeros(lx*ly,lxx*lyy);
        if kx<lhx
            for id=1:kx
                r2o=r2o+ r3o(1+(kx-1)*lx*ly-(id-1)*lx*ly:kx*lx*ly-(id-1)*lx*ly,1+(id-1)*lxx*lyy:lxx*lyy+(id-1)*lxx*lyy)/kx;
            end
        else
            for id=1:(nhx-kx+1)
                r2o=r2o+ r3o((lhx-1)*lx*ly+1-(id-1)*lx*ly:lhx*lx*ly-(id-1)*lx*ly,(kx-lhx)*lxx*lyy+1+(id-1)*lxx*lyy:(kx-lhx+1)*lxx*lyy+(id-1)*lxx*lyy)/(nhx-kx+1);
            end
        end
        for j=1:ny
            if j<ly
                for id=1:j
                    dout(:,j,kx,ky) =dout(:,j,kx,ky)+ ave_antid(r2o(1+(j-1)*lx-(id-1)*lx:j*lx-(id-1)*lx,1+(id-1)*lxx:lxx+(id-1)*lxx))/j;
                end
            else
                for id=1:(ny-j+1)
                    dout(:,j,kx,ky) =dout(:,j,kx,ky)+ ave_antid(r2o((ly-1)*lx+1-(id-1)*lx:ly*lx-(id-1)*lx,(j-ly)*lxx+1+(id-1)*lxx:(j-ly+1)*lxx+(id-1)*lxx))/(ny-j+1);
                end
            end
        end
    end
end
return


function [dout] =ave_antid(din);
% averaging along antidiagonals

[n1,n2]=size(din);
nout=n1+n2-1;
dout=zeros(nout,1);
for i=1:nout
    if i<n1
        for id=1:i
            dout(i)=dout(i) + din(i-(id-1),1+(id-1))/i;
        end
    else
        for id=1:nout+1-i
            dout(i)=dout(i) + din(n1-(id-1),1+(i-n1)+(id-1))/(nout+1-i);
        end
    end
end
return







