function [] = DrawCircle(x,y, r,label)
% ����Բ
%   x:����x
%   y:����y
%   r:�뾶
%   label:�������ֻ���ͬ����ɫ

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

