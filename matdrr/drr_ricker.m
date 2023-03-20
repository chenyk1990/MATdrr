function [w,tw] = drr_ricker(f,dt,tlength)
% drr_ricker: Ricker wavelet of central frequency f0.
%
% INPUT:
% f : central freq. in Hz (f <<1/(2dt) )
% dt: sampling interval in sec
% tlength : the duration of wavelet in sec
%
% OUTPUT: 
% w:  the Ricker wavelet
% tw: time axis
%
% Example
%
%   [w,tw] = drr_ricker(10,0.004,0.2);
%    plot(tw,w);
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

if nargin==3
    nw=floor(tlength/dt)+1;
else
    nw=2.2/f/dt;
    nw=2*floor(nw/2)+1;
end
nc=floor(nw/2);
w = zeros(nw,1);

k=[1:1:nw]';

alpha = (nc-k+1).*f*dt*pi;
beta=alpha.^2;
w = (1.-beta.*2).*exp(-beta);

if nargout>1;
    tw = -(nc+1-[1:1:nw])*dt;
end

