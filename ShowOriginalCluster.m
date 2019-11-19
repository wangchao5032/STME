function [ showlegend] = ShowOriginalCluster(D)
%% 函数描述： 显示原始数据聚类效果
%   D：数据集

%% 函数主体
figure();
showlegend=[];% 显示图例的图形

clusterA = D((D(:, 4) == 0), :);
showlegend = [showlegend,scatter3(clusterA(:, 2),clusterA(:,3),clusterA(:, 1), 'r','.')];
hold on;
clusterB = D((D(:, 4) == 1), :);
showlegend = [showlegend,scatter3(clusterB(:, 2),clusterB(:,3),clusterB(:, 1), 'b','.')];
% hold on;
% clusterC = D((D(:, 4) == 2), :);
% scatter3(clusterC(:, 2),clusterC(:,3),clusterC(:, 1), 'b','*');  
% scatter3(D(:, 2),D(:,3), D(:, 1), 'k','.');  

grid on;    % 画出网格
end

