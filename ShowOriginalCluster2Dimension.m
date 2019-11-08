function [ ] = ShowOriginalCluster2Dimension( D )
%% 函数描述： 显示原始数据聚类效果
%   D：数据集

%% 函数主体
figure();

for j = 1 : length(D(:,1))
    if D(j,4) == 0    %如果第j个点是A类型的
        scatter(D(j, 2), D(j, 3),['.', 'r']);
    else    %如果第j个点是B类型的
        scatter(D(j, 2), D(j, 3),['.', 'b']);
    end
    hold on;
end

% clusterA = D((D(:, 4) == 0), :);
% scatter(clusterA(:, 2),clusterA(:,3), 'r','.');  
% hold on;
% clusterB = D((D(:, 4) == 1), :);
% scatter(clusterB(:, 2),clusterB(:,3), 'b','.');  

% scatter(D(:, 2),D(:,3), 'k','.');  
grid on;    % 画出网格
% axis([0,100,0,100]);
end

