function [] = drr_plot3d(d3d,frames,z,x,y)
%%drr_plot3d: plot beautiful 3D slices
% 
% By Yangkang Chen
% Jan, 2022
% 
% INPUT
% d3d: 3D input
% frames: vecotr [3x1], frame1,2,3
% z:axis
% x:axis
% y:axis
%
% OUTPUT
% NO
%
% DEMO
% demos/test_matdrr_drr3d_diffraction.m
% 
% REFERENCE
% If you find this plotting function useful, please cite the following
% paper to recognize the authors' credit, where this script was originally created.
% 
% Chen, Y., S. Fomel, and R. Abma, 2022, Joint deblending and source time inversion, 88(1), WA27–WA35.
% 

% load data3ddb.mat
data=d3d;
% figure;yc_imagesc(data(:,:,60));

% load data3ddb0.mat
[nz,nx,ny]=size(d3d);
% d=d3d;

if nargin==2
z=[0:nz-1]*0.004;
x=[1:nx];
y=[1:ny];
end

if nargin==1
frames=[round(nz/2),round(nx/2),round(ny/2)];
z=[0:nz-1]*0.004;
x=[1:nx];
y=[1:ny];
end

f1=frames(1);
f2=frames(2);
f3=frames(3);

data(1,:,:)=d3d(f1,:,:);
data(:,end,:)=d3d(:,f2,:);
data(:,:,1)=d3d(:,:,f3);

dd=data;
%% the following ploting is separate

%% plot isosurface
%% given a input 3D volume of mask (1 for salt, 0 for non-salt): d
%% given a input 3D seismic volume: dd
% randn('state',202020);
% dd=randn(nz,nx,ny);
% %key step
% %shiftdim
% [nz,nx,ny]=size(d);
% z=linspace(-1,1,nz);
% x=linspace(-1,1,nx);
% y=linspace(-1,1,ny);

%shift data
% d2=shiftdim(d,1);%fix d, the creation and plot are separated
% d2=drr_transp(d2,12);

dd2=shiftdim(dd,1);%fix d, the creation and plot are separated
dd2=drr_transp(dd2,12);


[x2,y2,z2]=meshgrid(x,y,z);
% x2=shiftdim(x1,1);
% y2=shiftdim(y1,1);
% z2=shiftdim(z1,1);
[nx,ny,nz]=size(dd2);
% figure('units','normalized','Position',[0.01 0.01 0.45, 0.5],'color','w');
% slice(x2,y2,z2,dd2,[nx/2],[ny/2],[nz/2]); caxis([-0.001,0.001]);shading interp;
x0=[min(x2(1,:,1)),max(x2(1,:,1))];
y0=[min(y2(:,1,1)),max(y2(:,1,1))];
z0=[min(z2(1,1,:)),max(z2(1,1,:))];
h=slice(x2,y2,z2,dd2,[x0(2)],[y0(1)],[z0(1)]); caxis([-0.01,0.01]);shading interp;
set(h,'FaceColor','interp','FaceAlpha','interp');
% set(h,'FaceColor','interp','FaceAlpha',0);
alpha('color');alpha(1);

xlabel('x');
ylabel('y');
zlabel('z');
% set(gca,'YDir','reverse');
set(gca,'ZDir','reverse');
colormap(cseis);
% view(gca,[63 31]);
view(gca,[37 19]);
hold on;

[nx]=length(x(:));
[ny]=length(y(:));
[nz]=length(z(:));
plot3(x,y0(1)*ones(nx,1),z(f1)*ones(nx,1),'b-','linewidth',2);
plot3(x0(2)*ones(ny,1),y,z(f1)*ones(ny,1),'b-','linewidth',2);

plot3(x(f2)*ones(nz,1),y0(1)*ones(nz,1),z,'b-','linewidth',2);
plot3(x(f2)*ones(ny,1),y,z0(1)*ones(ny,1),'b-','linewidth',2);


plot3(x0(2)*ones(nz,1),y(f3)*ones(nz,1),z,'b-','linewidth',2);
plot3(x,y(f3)*ones(nx,1),z0(1)*ones(nx,1),'b-','linewidth',2);


xlim([min(x),max(x)]);
ylim([min(y),max(y)]);
zlim([min(z),max(z)]);

ylabel('Y','Fontsize',16,'fontweight','bold');
xlabel('X','Fontsize',16,'fontweight','bold');
zlabel('Z','Fontsize',16,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',16,'Fontweight','normal');


return

function [map]=cseis()

map = [[0.5*ones(1,40),linspace(0.5,1,88),linspace(1,0,88),zeros(1,40)]',[0.25*ones(1,40),linspace(0.25,1,88),linspace(1,0,88),zeros(1,40)]',[zeros(1,40),linspace(0.,1,88),linspace(1,0,88),zeros(1,40)]'];

return

function [dout]=drr_transp(din,plane)
% drr_transp: Transpose two axes in a dataset
% by Yangkang Chen, Dec 18, 2019
% Modified on Jan, 2020
% 
% INPUT
% din: input dataset
% plane: Two-digit number with axes to transpose. The default is 12
% OUTPUT
% dout: output dataset
%
% DEMO:
% a=magic(3);b=reshape(a,3,1,3);c=drr_transp(b,23);norm(a-c)

if nargin==1
    plane=12;
end

[n1,n2,n3,n4,n5]=size(din);

switch plane
    
    case 12
        dout=zeros(n2,n1,n3,n4,n5);
        
        for i5=1:n5
            for i4=1:n4
                for i3=1:n3
                    dout(:,:,i3,i4,i5)=din(:,:,i3,i4,i5).';
                end
            end
        end
        
    case 23
        dout=zeros(n1,n3,n2,n4,n5);
        for i5=1:n5
            for i4=1:n4
                for i1=1:n1
                    dout(i1,:,:,i4,i5)=squeeze(din(i1,:,:,i4,i5)).';
                end
            end
        end
    case 13
        dout=zeros(n3,n2,n1,n4,n5);
        for i5=1:n5
            for i4=1:n4
                for i2=1:n2
                    dout(:,i2,:,i4,i5)=squeeze(din(:,i2,:,i4,i5)).';
                end
            end
        end
        
     case 14
        dout=zeros(n4,n2,n3,n1,n5);
        for i5=1:n5
            for i3=1:n3
                for i2=1:n2
                    dout(:,i2,i3,:,i5)=squeeze(din(:,i2,i3,:,i5)).';
                end
            end
        end
        
     case 15
        dout=zeros(n5,n2,n3,n4,n1);
        for i4=1:n4
            for i3=1:n3
                for i2=1:n2
                    dout(:,i2,i3,i4,:)=squeeze(din(:,i2,i3,i4,:)).';
                end
            end
        end        
        
    otherwise
        error('Invalid argument value.');
end

return
