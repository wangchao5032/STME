function [ ] = ShowClusterXZ( D )
%% 函数描述： 显示聚类效果
%   D：数据集
%% 函数主体
figure();
% colorTable = ['b','r','y','g','m','c','b','r','y','g','m','c','b','r','y','g','m','c','b','r','y','g','m','c','b','r','y','g','m','c','b','r','y','g','m','c','b','r','y','g','m','c','b','r','y','g','m','c','b','r','y','g','m','c','b','r','y','g','m','c','b','r','y','g','m','c','b','r','y','g','m','c','b','r','y','g','m','c','b','r','y','g','m','c','b','r','y','g','m','c'];
% shapeTable = ['o','+','*','x','s','d','^','v','>','<','p','h','o','+','*','x','s','d','^','v','>','<','p','h','o','+','*','x','s','d','^','v','>','<','p','h','o','+','*','x','s','d','^','v','>','<','p','h','o','+','*','x','s','d','^','v','>','<','p','h','o','+','*','x','s','d','^','v','>','<','p','h','o','+','*','x','s','d','^','v','>','<','p','h','o','+','*','x','s','d','^','v','>','<','p','h'];
colorTable = ['b','r','b','r','m','c','r','g','m','m','r','c','b','g','c','r','m','c','b','r','y','g','m','c','b','r','y','g','m','c','b','r','y','g','m','c','b','r','y','g','m','c','b','r','y','g','m','c','b','r','y','g','m','c','b','r','y','g','m','c','b','r','y','g','m','c','b','r','y','g','m','c','b','r','y','g','m','c','b','r','y','g','m','c','b','r','y','g','m','c'];
shapeTable = ['o','+','*','x','s','d','^','v','>','<','p','h','+','p','d','s','d','p','*','*','p','+','<','s','o','+','*','x','s','d','^','v','>','<','p','h','o','+','*','x','s','d','^','v','>','<','p','h','o','+','*','x','s','d','^','v','>','<','p','h','o','+','*','x','s','d','^','v','>','<','p','h','o','+','*','x','s','d','^','v','>','<','p','h','o','+','*','x','s','d','^','v','>','<','p','h'];
% 噪声点画小圆点
cluster = D((D(:, 5) == -1), :);
scatter(cluster(:,2),cluster(:,1),20,'.','k');  
hold on;

% 画出已经成功聚类的点，颜色按照colorTable，形状按照shapeTable
clusterNum=max(D(:,5));

for i=1:clusterNum
    cluster = D((D(:, 5) == i), :);
    for j = 1 : length(cluster(:,1))
        if cluster(j,4) == 0    %如果第j个点是A类型的
            scatter3(cluster(j, 2), cluster(j, 1),cluster(j, 1),10,'+',colorTable(i));
        else    %如果第j个点是B类型的
            scatter3(cluster(j, 2), cluster(j, 1),cluster(j, 1),10,'<', colorTable(i));
        end
    end
    hold on;
end


% for i=1:clusterNum
%     cluster = D((D(:, 5) == i), :);     %第i个簇
%     scatter(cluster(:, 2),cluster(:, 1),shapeTable(i),colorTable(i));
%         hold on;
% end
% axis([0,100,0,100,0,100]);
grid on;    % 画出网格

end

