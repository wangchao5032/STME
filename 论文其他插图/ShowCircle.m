function [showlegend] = ShowCircle(Circles, Radius)
%Circles：圆心列表 
%Radius：半径列表
%legend：需要绘制图例的图形
%   绘制背景的圆形

MinVal = -1;
MaxVal = 1;
MaxRadius = 0.5;
t = 0 : .1 : 2 * pi;

% Circles = [8,0;0,0;1,0]; % 圆心坐标
% Radius = [8; 4]; % 半径
[nCircles,Dimension] = size(Circles); % 圆的数目，圆的维度

% cmap = hsv(nCircles); %// define colors. You could change `hsv` to `jet`, `cool`, ...
% colorTable = hsv2rgb([linspace(0, 0.9, nCircles)' linspace(0.99, 0.99, nCircles)' linspace(0.99, 0.99, nCircles)']);
colorTable = ['g','y','g','m','b'];
alpha = .3; % 透明度

showlegend = [];
for i = 1 : length(Radius)
    x = Radius(i) * cos(t) + Circles(i,1);
    y = Radius(i) * sin(t) + Circles(i,2);
    %     patch(x, y, cmap(i,:), 'facealpha', alpha, 'edgecolor', 'none'); %// plot filled circle with transparency
    showlegend(i)=patch(x, y, colorTable(i), 'facealpha', alpha, 'edgecolor', 'none'); %// plot filled circle with transparency
end

% axis equal; %// same aspect ratio in both axes
% axis([-5,17,-9,9]);
% set(gcf,'WindowStyle','normal');
% grid on;


end

