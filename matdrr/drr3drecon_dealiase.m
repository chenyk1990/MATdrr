function [ D1 ] = drr3drecon_dealiase(D,flow,fhigh,dt,N,K,Niter,a,lambda,verb)
%  DRR3DRECON: DRR for 3D seismic reconstruction (or F-XY domain damped multichannel singular spectrum analysis (DMSSA) for
%  simultaneous denoising and reconstruction)
%
%  IN   D:   	 intput 3D data
%       MASK:   sampling mask (consistent with the POCS based approaches)
%       flow:   processing frequency range (lower)
%       fhigh:  processing frequency range (higher)
%       dt:     temporal sampling interval
%       N:      number of singular value to be preserved
%       K:     damping factor (default: 4)
%       Niter:  number of maximum iteration
%       a:      interpolation interval
%       lambda: parameter to control the increasing rate of the weighting
%       verb:   verbosity flag (default: 0)
%
%  OUT  D1:  	output data
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
%  [0] Huang, W., D. Feng, and Y. Chen, 2020, De‐aliased and de‐noise Cadzow filtering for seismic data reconstruction, Geophysical Prospecting, 68, 443-571.
%  [1] Chen, Y., W. Huang, D. Zhang, W. Chen, 2016, An open-source matlab code package for improved rank-reduction 3D seismic data denoising and reconstruction, Computers & Geosciences, 95, 59-66.
%  [2] Chen, Y., D. Zhang, Z. Jin, X. Chen, S. Zu, W. Huang, and S. Gan, 2016, Simultaneous denoising and reconstruction of 5D seismic data via damped rank-reduction method, Geophysical Journal International, 206, 1695-1717.
%  [3] Huang, W., R. Wang, Y. Chen, H. Li, and S. Gan, 2016, Damped multichannel singular spectrum analysis for 3D random noise attenuation, Geophysics, 81, V261-V270.
%  [4] Chen et al., 2017, Preserving the discontinuities in least-squares reverse time migration of simultaneous-source data, Geophysics, 82, S185-S196.
%  [5] Chen et al., 2019, Obtaining free USArray data by multi-dimensional seismic reconstruction, Nature Communications, 10:4434.

if nargin==0
    error('Input data must be provided!');
end

if nargin==1
    error('Samping mask should be given');
end

if nargin==2
    flow=1;
    fhigh=124;
    dt=0.004;
    N=1;
    K=4;
    Niter=30;
    a=2;
    lambda=2;
    verb=0;
end;

[nt,nx,ny]=size(D);
D1=zeros(nt,nx,ny);

nf=2^nextpow2(nt);

% Transform into F-X domain
DATA_FX=fft(D,nf,1);
DATA_FX0=zeros(nf,nx*a,ny*a);

% First and last nts of the DFT.
ilow  = floor(flow*dt*nf)+1;

if ilow<1;
    ilow=1;
end;

ihigh = floor(fhigh*dt*nf)+1;

if ihigh > floor(nf/2)+1;
    ihigh=floor(nf/2)+1;
end

%% key sizes
lx=floor(nx/2)+1;
ly=floor(ny/2)+1;
lxx=nx-lx+1;
lyy=ny-ly+1;
lx2=nx*(a-1)+lx;
ly2=ny*(a-1)+ly;
M=zeros(lx*ly,lxx*lyy);

%% define G matrix
G=zeros(nx*a,nx);
GG=zeros(nx*a,nx);
for i=1:nx
    G(1+(i-1)*a,i)=1;
    GG(2+(i-1)*a:i*a,i)=ones(a-1,1);
end
G2=zeros(ny*a,ny);
GG2=zeros(ny*a,ny);
G2=zeros(ny*a,ny);
GG2=zeros(ny*a,ny);
for i=1:ny
    G2(1+(i-1)*a,i)=1;
    GG2(2+(i-1)*a:i*a,i)=ones(a-1,1);
end



% main loop
for k=ilow:ihigh
    
    fa_re=zeros(ny,nx);
    i_low=round(k/a);
    temp_low=DATA_FX(i_low,:,:);
    fa_low=reshape(temp_low,[nx,ny]);
    fa_low=fa_low';
    A_low=P_H(fa_low.',lx,ly).';
    [U_low,S_low,V_low]=svd(A_low);

    temp=DATA_FX(k,:,:);
    fa=reshape(temp,[nx,ny]);
    fa_new1=G*fa;
    fa_new1=fa_new1';
    fa_new=G2*fa_new1;
    fa_iter=fa_new;
    
    for iter=1:Niter
        A=h_hankel(fa_iter,lx2,ly2);
        M=U_low(:,1:K)*U_low(:,1:K)'*A;
        fa=anti_Hankel_2(M,lx2,lxx,ly2,lyy,nx*a,ny*a);
        fa_re=fa;
        for alj=1:a:ny*a
            for ali=1:a:nx*a
                fa_re(alj,ali)=0;
            end
        end

        if iter==1
            ra=0;
        else
            ra=(1/(Niter-iter+1)).^lambda;
        end
        fa_iter=(ra)*(fa-fa_re)+(1-ra)*fa_new+fa_re;
    end

    DATA_FX0(k,:,:)=reshape(fa_iter',[1,nx*a,ny*a]);

    if(mod(k,5)==0 && verb==1)
        fprintf( 'F %d is done!\n\n',k);
    end
    
end

% Honor symmetries
for k=nf/2+2:nf
    DATA_FX0(k,:,:) = conj(DATA_FX0(nf-k+2,:,:));
end

% Back to TX (the output)
D1=real(ifft(DATA_FX0,[],1));
D1=D1(1:nt,:,:);

return

function [dout]=P_H(din,lx,ly)
% forming block Hankel matrix
% size(din)
[nx,ny]=size(din);
lxx=nx-lx+1;
lyy=ny-ly+1;

for j=1:ny
    r=hankel(din(1:lx,j),[din(lx:nx,j)]);
    if j<ly
        for id=1:j
            dout(1+(j-1)*lx-(id-1)*lx:j*lx-(id-1)*lx,1+(id-1)*lxx:lxx+(id-1)*lxx) = r;
        end
    else
        for id=1:(ny-j+1)
            dout((ly-1)*lx+1-(id-1)*lx:ly*lx-(id-1)*lx,(j-ly)*lxx+1+(id-1)*lxx:(j-ly+1)*lxx+(id-1)*lxx)=r;
        end
    end
end
return

function [dout]=P_RD(din,N,K)
% Rank reduction on the block Hankel matrix


%      [U,D,V]=svds(din,N); % a little bit slower for small matrix
%      dout=U*D*V';
% %
[U,D,V]=svd(din);
for j=1:N
    D(j,j)=D(j,j)*(1-D(N+1,N+1)^K/(D(j,j)^K+0.000000000000001));
end

dout=U(:,1:N)*D(1:N,1:N)*(V(:,1:N)');

return

function [dout]=P_A(din,nx,ny,lx,ly)
% Averaging the block Hankel matrix to output the result
lxx=nx-lx+1;
lyy=ny-ly+1;
dout=zeros(nx,ny);

for j=1:ny
    if j<ly
        for id=1:j
            dout(:,j) =dout(:,j)+ ave_antid(din(1+(j-1)*lx-(id-1)*lx:j*lx-(id-1)*lx,1+(id-1)*lxx:lxx+(id-1)*lxx))/j;
        end
    else
        for id=1:(ny-j+1)
            dout(:,j) =dout(:,j)+ ave_antid(din((ly-1)*lx+1-(id-1)*lx:ly*lx-(id-1)*lx,(j-ly)*lxx+1+(id-1)*lxx:(j-ly+1)*lxx+(id-1)*lxx))/(ny-j+1);
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


function A=h_hankel(fa,Nx,Lx)
%Construct the hankel matrix
%A=hankel(fa,50,40), A=hankel(fa,51,38)
if nargin==2
    Lx=1;
end
[n,m]=size(fa);
Ny=m+1-Nx;
Ly=n+1-Lx;
for i=1:n
    for j=1:Ny     
        H(i,(j-1)*Nx+1:j*Nx)=fa(i,j:Nx+j-1);%Hankel matrix
    end
end
for k=0:(Ly-1)
    for i=1:Lx
        for j=1:Ny
            A((j+k*Ny),(i-1)*Nx+1:i*Nx)=H(i+k,(j-1)*Nx+1:j*Nx);
        end
    end
end
return

function fa=anti_Hankel_2(M,Nx,Ny,Lx,Ly,n,m)
a=Nx-Ny;
b=Lx-Ly;
M_new=zeros(Ny*Ly,Nx*Lx+(Ly-1)*Nx);
temp=zeros(Ny,Nx,m);
JS=ones(Ly,Lx);
JS_new=zeros(Ly,Lx+Ly-1);
for i=1:Ly
    M_new(1+(i-1)*Ny:i*Ny,1+(i-1)*Nx:Nx*Lx+(i-1)*Nx)=M(1+(i-1)*Ny:i*Ny,:);
    JS_new(i,i:Lx+i-1)=JS(i,:);
end
JS_1=sum(JS_new);
for i=1:Lx+Ly-1
    sump=zeros(Ny,Nx);
    for j=1:Ly
        sump=sump+M_new(1+(j-1)*Ny:j*Ny,1+(i-1)*Nx:i*Nx);
    end
    temp(:,:,i)=sump/JS_1(1,i);
end

N_new=zeros(Ny,Nx+Ny-1);
JS=ones(Ny,Nx);
JS_new=zeros(Ny,Nx+Ny-1);
for i=1:m
    w=temp(:,:,i);
    for j=1:Ny
        JS_new(j,j:Nx+j-1)=JS(j,:);
        N_new(j,j:Nx+j-1)=w(j,:);
    end
    fa(i,:)=sum(N_new)./sum(JS_new);
end

return





