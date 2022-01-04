function [showlegend ] = ShowOriginalCluster2Dimension( D )
%% 函数描述： 显示原始数据聚类效果
%   D：数据集

%% 函数主体
figure();
showlegend=[];% 显示图例的图形

clusterA = D((D(:, 4) == 0), :);
showlegend = [showlegend,scatter(clusterA(:, 2),clusterA(:,3), 'r','.')];  
hold on;
clusterB = D((D(:, 4) == 1), :);
showlegend = [showlegend,scatter(clusterB(:, 2),clusterB(:,3), 'b','.')];  
hold on;
clusterC = D((D(:, 4) == 2), :);
showlegend = [showlegend,scatter(clusterC(:, 2),clusterC(:,3), 'k','.')];  

% scatter(D(:, 2),D(:,3), 'k','.');  
grid on;    % 画出网格

end

