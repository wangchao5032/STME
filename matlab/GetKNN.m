function [ distK_all, distKUpID,distKList] = GetKNN( D,k,t )
%% 函数说明 计算所有点的最近的K个邻居
%   D:数据
%   k：第k个邻居
%   t:时间窗口
%   distK_all:每个点的最近的k个邻居
%% 函数主体
distK_all = {}; %所有点的k个邻居

% 2019-11--------------------------
types_num = max(D(:, 4)) + 1;  % 类型的编号从0开始，为0 1 2，所以这里要多加一个1，即一共有3种类型  
% distKList 从一个1*n的向量变成types_num*n的矩阵
% 第一行为不考虑类型的k距离，第2至第4行为只考虑A类型，B类型C类型的k距离
distKList = zeros(types_num + 1, length(D(:, 1)));   
% --------------------------


for i=1 : (length(D(:,1)))
    temp = inf * ones(1, length(D(:, 1)));
    tempID = [];
    dist = [];
    for j = 1 : (length(D(:,1)))
        if i~=j && (abs(D(i,1)-D(j,1)) <= t)    %i不等于j，且两者时间间隔小于等于t
             dist_i_j = sqrt((D(i, 2) - D(j,2))^2 + (D(i, 3) - D(j,3))^2);   %计算i和j的空间距离
        else
            dist_i_j = +inf;    %i和j的空间距离为正无穷
        end
        temp(j) = dist_i_j;
    end
    [dist,tempID] = sort(temp);     %对第i个点距离所有点的距离排序
    distK_all{i} = tempID(1:k);   %保存距离最近的k个邻居的ID
    distKList(1, i) = dist(k);     %第i个点第k个最近的距离
    
    % 2019-11 ----------------------------------
    % 计算考虑类型的k距离，存储在distKList的后3行
    for q = 0:types_num-1
        temp_typet = temp(D(:, 4) == q);
        dist_t = mink(temp_typet, k);   % mink() 函数可以计算temp_typet向量中第k小的值
        if length(dist_t) == k
            distKList(q+2, i) = dist_t(k); 
        else
            distKList(q+2, i) = Inf;
        end
    end
    % ----------------------------------
    
end
[distKUp, distKUpID] = sort(distKList(1, :)); %对所有点第k个最近的距离进行排序

end

