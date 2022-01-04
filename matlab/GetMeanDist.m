function [ distmean ] = GetMeanDist( Data )
%% 函数说明 计算簇内对象的平均距离
%   Data:簇数据
%% 函数内容
distList = [];
for i = 1 : length(Data(:,1))
    x_i = Data(i,2);
    y_i = Data(i,3);
    for j = 1 : length(Data(:,1))
        if(i == j)
            continue;
        end
        x_j = Data(j,2);
        y_j = Data(j,3);
        dist_i_j = sqrt((x_i - x_j)^2 + (y_i - y_j)^2);
        distList = [distList,dist_i_j];
    end
end
distmean = mean(distList);

% distList = [];
% xCenter = mean(Data(:,2));
% yCenter = mean(Data(:,3));
% % xCenter = (max(Data(:,2))+ min(Data(:,2)))/2;
% % yCenter = (max(Data(:,3))+ min(Data(:,3)))/2;
% for i=1 : (length(Data(:,1)))
%     
%     dist = sqrt((Data(i, 2) - xCenter)^2 + (Data(i, 3) - yCenter)^2);   %计算第i个点到中心点的空间距离
%     distList = [distList,dist];
% end
% distmean = mean(distList);
end

