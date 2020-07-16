function MyQuiver(xx,yy,uu,vv,scale,varargin)
%MYQUIVER show vector in specific scale 
%
% coded by Xin Liu
% email: liuxin2018@zju.edu.cn
% Apr.15, 2020

switch nargin
    case 5
        Color = 'k';
    case 6
        Color = varargin{1};
end

xs = xx(1:scale:length(xx));
ys = yy(1:scale:length(yy));
us = uu(1:scale:length(yy),1:scale:length(xx));
vs = vv(1:scale:length(yy),1:scale:length(xx));

quiver(xs,ys,us,vs,Color,'MaxHeadSize',1);
end