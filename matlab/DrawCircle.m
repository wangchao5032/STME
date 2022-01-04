function [] = DrawCircle(x,y, r,label)
% 绘制圆
%   x:中心x
%   y:中心y
%   r:半径
%   label:用于区分画不同的颜色

theta=0:2*pi/3600:2*pi;
Circle1=x+r*cos(theta);
Circle2=y+r*sin(theta);
if label == 1
    plot(Circle1,Circle2,'r','Linewidth',1);
else
     plot(Circle1,Circle2,'b','Linewidth',1);
end

hold on;

axis equal
end

