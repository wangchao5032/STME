function [ showlegend ] = ShowCluster2Dimension( D,real_clusterLabel )
%% ���������� ��ʾ����Ч��
%   D�����ݼ�
% figure();
%% ��ɫ�б�
%colorTable = ['b','r','y','y','m','c','r','y','g','m','b','r','y','g','m','b','r','y','g','m','b','r','y','g','m','b','r','y','g','m','b','r','y','g','m'];
%colorTable = ['b','r','y','g','m','c','b','r','y','g','m','c','b','r','y','g','m','c','b','r','y','g','m','c','b','r','y','g','m','c','b','r','y','g','m','c','b','r','y','g','m','c','b','r','y','g','m','c','b','r','y','g','m','c','b','r','y','g','m','c','b','r','y','g','m','c','b','r','y','g','m','c','b','r','y','g','m','c','b','r','y','g','m','c','b','r','y','g','m','c'];
% colorTable = ['y','b','r','g','m','c','r','c','g','r','b','m','c','r','r','y','g','m','c','b','c','y','m','m','r','b','c','g','g','m','r','b','r','y','g','m','c','b','r','y','g','m','c','b','r','y','g','m','c','b','r','y','g','m','c','b','r','y','g','m','c','b','r','y','g','m','c','b','r','y','g','m','c','b','r','y','g','m','c','b','r','y','g','m','c','b','r','y','g','m','c'];
colorTable = hsv2rgb([linspace(0, 0.9, real_clusterLabel+1)' linspace(0.99, 0.99, real_clusterLabel+1)' linspace(0.99, 0.99, real_clusterLabel+1)']);
% colorTable = hsv2rgb([linspace(0, 0.9, clusterNum+1)' linspace(0.99, 0.99, clusterNum+1)' linspace(0.99, 0.99, clusterNum+1)']);
%% �ص���Ŀ
clusterNum=max(D(:,5));  % �ϲ���С�Ĵ�֮ǰ�ص���Ŀ
real_clusterNum = 0; % ��ʵ�Ĵص���Ŀ
showlegend=[];% ��ʾͼ����ͼ��

%% �����㻭СԲ��
cluster = D((D(:, 5) <= 0), :);  
% scatter(cluster(:, 2),cluster(:,3), 15,'k','.');
showlegend = [showlegend, scatter(cluster(:, 2),cluster(:,3), 50,'k','.')]; %[0.5 0.5 0.5]
hold on

%% ���ƾ����� 
for i=1:clusterNum
    cluster = D((D(:, 5) == i), :);
    if ( length(cluster) == 0 )
        continue;
    end
    real_clusterNum = real_clusterNum + 1;
    cluster1 = cluster((cluster(:,4) == 0),:);  %��i�����е�A����
    cluster2 = cluster((cluster(:,4) == 1),:);  %��i�����е�B����
    cluster3 = cluster((cluster(:,4) == 2),:);  %��i�����е�C����
    cluster4 = cluster((cluster(:,4) == 3),:);  %��i�����е�C����
    % ------------------�ù���ɫ����--------------------
    showlegend = [showlegend,scatter(cluster1(:, 2), cluster1(:, 3), 50,colorTable(real_clusterNum, :),'.')]; % +
    scatter(cluster2(:, 2), cluster2(:, 3), 50,colorTable(real_clusterNum, :),'.'); % <
    scatter(cluster3(:, 2), cluster3(:, 3), 50,colorTable(real_clusterNum, :),'.'); % o
    scatter(cluster4(:, 2), cluster4(:, 3), 50,colorTable(real_clusterNum, :),'.'); % o
    % ------------------��ָ����ɫ����--------------------
%     showlegend = [showlegend,scatter(cluster1(:, 2), cluster1(:, 3), 30,[colorTable(real_clusterNum),'.'])]; % +
%     scatter(cluster2(:, 2), cluster2(:, 3), 50,[colorTable(real_clusterNum),'.']); % <
%     scatter(cluster3(:, 2), cluster3(:, 3), 50,[colorTable(real_clusterNum),'.']); % <
%     scatter(cluster4(:, 2), cluster4(:, 3), 50,[colorTable(real_clusterNum),'.']); % <
    hold on;
    
    %  �ڴصĵ�һ�����Ա߱�����֣�ע�⣬���ﲢ���ǵ�һ������صĵ㣩
%     if ~isempty(cluster)
%         text(cluster(1, 2), cluster(1, 3), ['c-' int2str(real_clusterNum)]);
%     end
end

% grid on;    % ��������
% axis([0,100,0,100]);
% axis([-5,17,-8,8]);
% set(gcf,'WindowStyle','normal');
% axis equal 
end

