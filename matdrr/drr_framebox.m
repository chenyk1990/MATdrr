function [] = drr_framebox(x1,x2,y1,y2,c,lw)
%
% for drawing a frame box
% 
% Jan, 5, 2023, By Yangkang Chen
% 
% x1,x2,y1,y2: intuitive
% 
if nargin==4
c='r';
lw=2;
end

hold on;


hold on;
plot([x1,x2],[y1,y1],'-','color',c,'linewidth',lw);
plot([x1,x2],[y2,y2],'-','color',c,'linewidth',lw);
plot([x1,x1],[y1,y2],'-','color',c,'linewidth',lw);
plot([x2,x2],[y1,y2],'-','color',c,'linewidth',lw);


return