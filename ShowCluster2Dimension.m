function [ showlegend ] = ShowCluster2Dimension( D,real_clusterLabel )
%% 函数描述： 显示聚类效果
%   D：数据集
% figure();
%% 颜色列表
%colorTable = ['b','r','y','y','m','c','r','y','g','m','b','r','y','g','m','b','r','y','g','m','b','r','y','g','m','b','r','y','g','m','b','r','y','g','m'];
%colorTable = ['b','r','y','g','m','c','b','r','y','g','m','c','b','r','y','g','m','c','b','r','y','g','m','c','b','r','y','g','m','c','b','r','y','g','m','c','b','r','y','g','m','c','b','r','y','g','m','c','b','r','y','g','m','c','b','r','y','g','m','c','b','r','y','g','m','c','b','r','y','g','m','c','b','r','y','g','m','c','b','r','y','g','m','c','b','r','y','g','m','c'];
% colorTable = ['y','b','r','g','m','c','r','c','g','r','b','m','c','r','r','y','g','m','c','b','c','y','m','m','r','b','c','g','g','m','r','b','r','y','g','m','c','b','r','y','g','m','c','b','r','y','g','m','c','b','r','y','g','m','c','b','r','y','g','m','c','b','r','y','g','m','c','b','r','y','g','m','c','b','r','y','g','m','c','b','r','y','g','m','c','b','r','y','g','m','c'];
colorTable = hsv2rgb([linspace(0, 0.9, real_clusterLabel+1)' linspace(0.99, 0.99, real_clusterLabel+1)' linspace(0.99, 0.99, real_clusterLabel+1)']);
% colorTable = hsv2rgb([linspace(0, 0.9, clusterNum+1)' linspace(0.99, 0.99, clusterNum+1)' linspace(0.99, 0.99, clusterNum+1)']);
%% 簇的数目
clusterNum=max(D(:,5));  % 合并碎小的簇之前簇的数目
real_clusterNum = 0; % 真实的簇的数目
showlegend=[];% 显示图例的图形

%% 噪声点画小圆点
cluster = D((D(:, 5) <= 0), :);  
% scatter(cluster(:, 2),cluster(:,3), 15,'k','.');
showlegend = [showlegend, scatter(cluster(:, 2),cluster(:,3), 50,'k','.')]; %[0.5 0.5 0.5]
hold on

%% 绘制聚类结果 
for i=1:clusterNum
    cluster = D((D(:, 5) == i), :);
    if ( length(cluster) == 0 )
        continue;
    end
    real_clusterNum = real_clusterNum + 1;
    cluster1 = cluster((cluster(:,4) == 0),:);  %第i个簇中的A类型
    cluster2 = cluster((cluster(:,4) == 1),:);  %第i个簇中的B类型
    cluster3 = cluster((cluster(:,4) == 2),:);  %第i个簇中的C类型
    cluster4 = cluster((cluster(:,4) == 3),:);  %第i个簇中的C类型
    % ------------------用过渡色绘制--------------------
    showlegend = [showlegend,scatter(cluster1(:, 2), cluster1(:, 3), 50,colorTable(real_clusterNum, :),'.')]; % +
    scatter(cluster2(:, 2), cluster2(:, 3), 50,colorTable(real_clusterNum, :),'.'); % <
    scatter(cluster3(:, 2), cluster3(:, 3), 50,colorTable(real_clusterNum, :),'.'); % o
    scatter(cluster4(:, 2), cluster4(:, 3), 50,colorTable(real_clusterNum, :),'.'); % o
    % ------------------用指定颜色绘制--------------------
%     showlegend = [showlegend,scatter(cluster1(:, 2), cluster1(:, 3), 30,[colorTable(real_clusterNum),'.'])]; % +
%     scatter(cluster2(:, 2), cluster2(:, 3), 50,[colorTable(real_clusterNum),'.']); % <
%     scatter(cluster3(:, 2), cluster3(:, 3), 50,[colorTable(real_clusterNum),'.']); % <
%     scatter(cluster4(:, 2), cluster4(:, 3), 50,[colorTable(real_clusterNum),'.']); % <
    hold on;
    
    %  在簇的第一个点旁边标记文字（注意，这里并不是第一被加入簇的点）
%     if ~isempty(cluster)
%         text(cluster(1, 2), cluster(1, 3), ['c-' int2str(real_clusterNum)]);
%     end
end

% grid on;    % 画出网格
% axis([0,100,0,100]);
% axis([-5,17,-8,8]);
% set(gcf,'WindowStyle','normal');
% axis equal 
end

