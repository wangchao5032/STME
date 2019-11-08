function [showlegend] = ShowCircle(Circles, Radius)
%Circles��Բ���б� 
%Radius���뾶�б�
%legend����Ҫ����ͼ����ͼ��
%   ���Ʊ�����Բ��

MinVal = -1;
MaxVal = 1;
MaxRadius = 0.5;
t = 0 : .1 : 2 * pi;

% Circles = [8,0;0,0;1,0]; % Բ������
% Radius = [8; 4]; % �뾶
[nCircles,Dimension] = size(Circles); % Բ����Ŀ��Բ��ά��

% cmap = hsv(nCircles); %// define colors. You could change `hsv` to `jet`, `cool`, ...
% colorTable = hsv2rgb([linspace(0, 0.9, nCircles)' linspace(0.99, 0.99, nCircles)' linspace(0.99, 0.99, nCircles)']);
colorTable = ['g','y','g','m','b'];
alpha = .3; % ͸����

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

