function [ showlegend] = ShowCluster( D,real_clusterLabel )
%% 函数描述： 显示聚类效果
%   D：数据集
%% 函数主体
figure();
shapeTable = ['o','+','*','x','s','d','^','v','>','<','p','h','o','+','*','x','s','d','^','v','>','<','p','h','o','+','*','x','s','d','^','v','>','<','p','h','o','+','*','x','s','d','^','v','>','<','p','h','o','+','*','x','s','d','^','v','>','<','p','h','o','+','*','x','s','d','^','v','>','<','p','h','o','+','*','x','s','d','^','v','>','<','p','h','o','+','*','x','s','d','^','v','>','<','p','h'];
colorTable = hsv2rgb([linspace(0, 0.9, real_clusterLabel+1)' linspace(0.99, 0.99, real_clusterLabel+1)' linspace(0.99, 0.99, real_clusterLabel+1)']);
%% 簇的数目
clusterNum=max(D(:,5)); % 合并碎小的簇之前簇的数目
real_clusterNum = 0; % 真实的簇的数目
showlegend=[];% 显示图例的图形
%% 噪声点画小圆点
cluster = D((D(:, 5) == -1), :);
showlegend = [showlegend,scatter3(cluster(:,2),cluster(:, 3),cluster(:,1),50,'.','k')];  
hold on;

% 画出已经成功聚类的点，颜色按照colorTable，形状按照shapeTable
for i=1:clusterNum
    cluster = D((D(:, 5) == i), :);
    if ( length(cluster) == 0 )
        continue;
    end
    real_clusterNum = real_clusterNum + 1;
    clusterA = cluster(cluster(:, 4) == 0, :);
    clusterB = cluster(cluster(:, 4) == 1, :);
    clusterC = cluster(cluster(:, 4) == 2, :);
    clusterD = cluster(cluster(:, 4) == 3, :);
    showlegend = [showlegend,scatter3(clusterA(:, 2), clusterA(:, 3), clusterA(:, 1),50,colorTable(real_clusterNum, :),'.')];
    scatter3(clusterB(:, 2), clusterB(:, 3), clusterB(:, 1),50,colorTable(real_clusterNum, :),'.');
    scatter3(clusterC(:, 2), clusterC(:, 3), clusterC(:, 1),50,colorTable(real_clusterNum, :),'.');
    scatter3(clusterD(:, 2), clusterD(:, 3), clusterD(:, 1),50,colorTable(real_clusterNum, :),'.');
    
%     if ~isempty(cluster)
%         text(cluster(1, 2), cluster(1, 3), cluster(1, 1), ['c-' int2str(real_clusterNum)]);
%     end
end


% for i=1:clusterNum
%     cluster = D((D(:, 5) == i), :);     %第i个簇
%     scatter3(cluster(:, 2),cluster(:, 3),cluster(:, 1),shapeTable(i),colorTable(i));
%         hold on;
% end
% axis([0,100,0,100,0,100]);
% grid on;    % 画出网格
% set(gcf,'WindowStyle','normal');

end

