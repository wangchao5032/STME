clc;clear;
format long ;   %[胡伟健添加]保留小数位数增加，便于判断
%% 从excel表格中读取数据
%D:模拟数据， 字段：时间，x,y,类型（0代表A类型，1代表B类型），簇标号（默认为0），可能还有点的理想标号
addpath('论文其他插图') 
% D = xlsread('simple_data.xlsx');
D = xlsread('北京市数据\博物馆-剧场-酒吧\博物馆-剧场-酒吧-data.xlsx');

% 2019-11 ------------------------------
% 生成新的随机的模拟数据
% D3 = []; D2 = D3; D1 = D2;
% D1 = generate_ring_data(0, 0,  0, 4, 0.3, 0.5, 0);
% D2 = generate_ring_data(6, 8,  0, 10, 0.9, 0.5, 1);
% D3 = generate_ring_data(16, -8, 0, 18, 2,   0.5, 2);
% D = [D1; D2; D3];
% --------------------------------------------
%% 绘制原始数据
% showlegend=[];
% showlegend = [showlegend, ShowOriginalCluster(D)];
% % 4、绘制图例
% legend(showlegend,'A-Type','B-Type');
% axis([-8,8,-8,8]);
% set(gcf,'WindowStyle','normal');
% grid on;    % 画出网格
% xlabel('x');
% ylabel('y');
% zlabel('t');
% ShowOriginalCluster2Dimension(D);
%% 输入参数
%d = Ripley_L(D); %聚类大小的直径 因为本算法并不固定簇的直径，所以不用计算
k = 13;
tWindow = 5;
kt = 8; % ceil(k*0.5);
MinPts = 8; %ceil(k*0.5);
distK_sigma_multi = 1; % 点数最多的类型可以放宽k距离的阈值
distK_sigma_times = 3; 
%% 计算每个点最近的K个邻居 和 所有点第k个邻居的升序排序后的距离
[knnList, distKUpList,distKList] = GetKNN(D,k,tWindow);  
distKList_typet = distKList(2:end, :);  % 考虑类型的k距离
distKList = distKList(1, :); % 不考虑类型的k距离
%% 计算每个点的共享近邻
RSTNList = {};   %共享近邻列表
for i = 1 : length(D(:,1))
    RSTN = [];      %共享近邻
    for j = 1 : k     %遍历k个邻居
        compareObject = knnList{i}(j);  %k个邻居中的第j个
        if i == compareObject
            continue;
        end
        %如果第j个k邻居和i的共享近邻个数大于kt，且i也是第j个k邻居的K邻居
        if length(intersect(knnList{i},knnList{compareObject})) >= kt && ismember(i,knnList{compareObject})==1  
             RSTN = [RSTN,compareObject];    %将该邻居加入共享近邻队列
        end
    end
    RSTNList{i} = RSTN;
end
%% 聚类过程
clusterLabel = 0;   %初始化簇标号
clusterDensityList = [];    %簇有效密度列表
clusterDensityListA = [];   %簇中A类型点的有效密度列表
clusterDensityListB = [];   %簇中B类型点的有效密度列表
clusterDensityListC = [];   %簇中C类型点的有效密度列表
clusterNumList = []; % 簇中点的总数 

clusterCoreList = []; % 簇的第一个核心点列表

labelKList = []; %按照标记的顺序记录被标记点的k距离
labelList = []; %按照标记的顺序记录被标记的簇号
labelFirstList = []; %簇中第一次被标记的点 是 所有数据中第几个被标记的点
for ii = 1 : length(D(:,1))     %按照第k个邻居距离的升序进行枚举，也就是第i个点一定是目前所有点里k邻域范围最小的。
    i = distKUpList(ii);    %第ii个“所有点第k个邻居距离的升序排序”的标号
    if D(i,5) ~= 0  %第i个点已经被处理过
        continue;
    end
    Y = knnList{i};     % 第i个点的k个邻居
    X = RSTNList{i};    % 第i个点的共享邻居/直接可达点
    RSTN_label = D(X,5); % 第i个点的共享近邻的标签
    RSTN_free_pts = nnz((RSTN_label == 0 | RSTN_label == -1)); % 第i个点的共享近邻中标签为0或-1的点
    RSTN_most_label = mode(RSTN_label); % 第i个点的共享近邻中出现最多的类型
    if length(X) < MinPts  % 如果第i个点是噪声，不是核心点
        D(i, 5) = -1; 
    elseif RSTN_free_pts < MinPts &&  RSTN_most_label ~=0
%         RSTN_most_label = mode(RSTN_label);
        D(i,5) = RSTN_most_label;
    else % 如果第i个点是核心点，则开始一个新的簇
        clusterLabel = clusterLabel + 1;    % 簇标记+1
        clusterValidList = [];  % 参与计算当前簇有效体积的集合，存放当前簇的所有核心点。和type_buf_core的作用一样，其实可以删除的。
        clusterCoreList = [clusterCoreList,i]; % 当前簇的第一个核心点
        labelFirstList = [labelFirstList,length(labelList)+1]; % 产生当前簇的第一个核心点时，已经有length(labelList)个点被标记过了，所以当前簇的第一个核心点即将是第length(labelList)+1个被标记的点
         % 构造初始队列
        queue = [i, RSTNList{i}];               %将第i个点及其RSTN加入队列
        distK_buf = distKList_typet(:, queue); % 队列中所有点的考虑类型的k距离
        type_buf_init = D(queue,4); %  簇的有效点1:簇的初始队列中的所有点
        type_buf_core=[]; %  簇的有效点2:后续加入队列的点中的核心点
        
        % 计算初始队列中所有点的均值和标准差
        distK_mu = mean(distK_buf, 2);                                          %计算“队列中所有点的考虑类型的k距离”中点的 均值
        distK_sigma_init = sqrt(var(distK_buf, 1, 2));                               %计算“队列中所有点的考虑类型的k距离”中点的 标准差
        
        cluster1_num = 0;
        % 遍历队列
        while ~isempty(queue) 
            %  为了绘图，临时加的，正式跑时记得删除
%             cluster1_num = cluster1_num+1;
%             if (clusterLabel==1 && cluster1_num==80) 
%                 break;
%             end
            % -----------------1、统计当前对列中各类型点的数目
            type_buf = [type_buf_init; type_buf_core]; % 簇的有效点：初始队列中所有点+后续加入队列的点中的核心点
            
            type_hist = [nnz(type_buf == 0) nnz(type_buf == 1) nnz(type_buf == 2)]; % 统计当前队列中有效点的所有类型
            type_valid = (type_hist >= MinPts);                   % 当前队列的有效点中，如果某种类型的点数目大于MinPts，则认为这些类型是有影响力的
            [~, dominate_type] = max(type_hist);                                    % 找出当前队列的有效点中，影响力最大的类型（文章中忽略）
            distK_sigma = distK_sigma_init;
            distK_sigma(dominate_type) = distK_sigma(dominate_type) * distK_sigma_multi;          % 影响力最大的类型可以放宽k距离的阈值（为了使模拟数据边界密度较低的点包含进来）（文章中忽略） 
            % -----------------2、取出当前点并标记
            ptCurrent = queue(1);   % 取出队列中的第一个点作为当前点
            queue(1) = []; % 把当前点移出队列
            if D(ptCurrent, 5) > 0  % 如果当前点被标记过，这个点只可能是作为共享近邻被加入队列了。
                continue;
            end 
            % 队列中的点有以下情况:
            % (1)初始簇形成时依赖的第1个点，一定是没有被标记也不是噪声的核心点。
            % (2)初始簇形成时依赖的第1个点的共享近邻，不一定是核心点，可能被标记过。会在后文进一步筛选是否是核心点（if length(X2) >= MinPts），如果是再将其k邻居加入队列
            % (3)后期加入队列中的点，一定是没有被标记过的点或噪声，不一定是核心点。会在后文进一步筛选是否是核心点（if length(X2) >= MinPts），如果是再将其k邻居加入队列
            D(ptCurrent, 5) = clusterLabel;     
            labelKList = [labelKList,distKList(ptCurrent)];
            labelList = [labelList,clusterLabel]; 
            % -----------------3、遍历当前点的k邻居，加入队列
            Y2 = knnList{ptCurrent};     % 队列中当前点的k个邻居
            X2 = RSTNList{ptCurrent};    % 队列中当前点的共享邻居/直接可达点
            if length(X2) >= MinPts  % 如果当前点是核心点
                type_buf_core = [type_buf_core; D(ptCurrent,4)];
                clusterValidList = [clusterValidList,ptCurrent];             
                for j = 1 : length(Y2)   % 遍历当前点的k邻居
                    distK_diff = abs(distKList_typet(:, Y2(j)) - distK_mu);   % 计算当前点第j个邻居的考虑类型的k距离（1*3） 和 初始簇考虑类型的k距离均值（1*3） 的差值（1*3）                 
                    % 如果当前点的第j个邻居(1)是噪声或还没有被标记过;(2)不在对列中;
                    % (3)在只考虑有影响力的类型下，当前点第j个邻居的考虑类型的k距离（1*2） 和 初始簇考虑类型的k距离均值（1*2）的差值distK_diff(type_valid)（1*2），若差值在初始簇的3倍标准差以内
                    % D(Y2(j),5) == -1 || 
                    if ( D(Y2(j),5) == -1 ||D(Y2(j),5) == 0 )  && isempty(find(queue == Y2(j), 1)) &&  all(distK_diff(type_valid) < distK_sigma(type_valid) * distK_sigma_times) 
                        queue = [queue, Y2(j)] ;    % 将当前点的第j个邻居加入队列
                    end
                end
            end
        end
        % -----------------4、计算当前簇的点数
        cluster = D((D(:, 5) == clusterLabel), :);     %当前簇
        clusterNumList = [clusterNumList, length(find(D(:, 5) == clusterLabel))]; % 当前簇的点的总数
        clusterNumA = nnz(cluster(:,4)==0); %当前簇中A类型点的数目
        clusterNumB = nnz(cluster(:,4)==1); %当前簇中B类型点的数目
        clusterNumC = nnz(cluster(:,4)==2); %当前簇中C类型点的数目
        
        % -----------------5、计算当前簇有效密度
        %计算三维的凸包
%           [K,V] = convhull(cluster(:,2),cluster(:,3),cluster(:,1));   %当前簇的最小凸包
        %计算二维的凸包
        [K,V] = convhull(cluster(:,2),cluster(:,3));   %当前簇的最小凸包
        
        clusterDensity =  length(cluster) / V ;    %当前簇的密度=当前簇的数量/当前簇的体积
        clusterDensityList = [clusterDensityList;clusterDensity];
        clusterDensityA = clusterNumA /V; %当前簇的A类型点的密度 = 当前簇中A类型点的数量/当前簇的体积    
        clusterDensityListA = [clusterDensityListA,clusterDensityA];
        clusterDensityB = clusterNumB /V; %当前簇的B类型点的密度 = 当前簇中B类型点的数量/当前簇的体积
        clusterDensityListB = [clusterDensityListB,clusterDensityB];
        clusterDensityC = clusterNumC /V; %当前簇的C类型点的密度 = 当前簇中C类型点的数量/当前簇的体积
        clusterDensityListC = [clusterDensityListC,clusterDensityC];
       
%         clusterA = []; %保存当前簇的所有核心点中的A类型
%         clusterB = []; %保存当前簇的所有核心点中的B类型
%         cluster = []; %保存当前簇的所有核心点
%         for p =1:length(clusterValidList)
%             if( D(clusterValidList(p),4)==0)
%                 clusterA = [clusterA;D(clusterValidList(p),:)];   
%             else
%                 clusterB = [clusterB;D(clusterValidList(p),:)];   
%             end
%         end
%         cluster = [clusterA;clusterB];
        
    end
end

%% 将碎小的簇合并至其他的簇
% merge_num = 0;
% for i = 1:length(clusterNumList)
%     if (clusterNumList(i) < MinPts) % 如果簇中点的总数比MinPts还小
%         [id,id_] = find(D(:, 5) == i); 
%         max_num_list = zeros(1,length(clusterNumList)); % 存放该簇的所有点的共享近邻中，各个簇号出现的频率
%         for j = 1:length(id) % 遍历该簇中的所有点，统计簇中大多数点最近的簇
%             RSTN_j = RSTNList{id(j)}; % 簇的第j个点的共享近邻
%             labels = D(RSTN_j,5); % 簇的第j个点的共享近邻的所有簇号
%             % 出现次数最多的簇号
%             temp = 1:max(labels); % 数组中可能存在的整数的集合
%             c = histc(labels,temp); % 各数字出现的次数
%             [max_num, max_index] = max(c); % 出现次数最多的元素出现的次数max_num，和该元素在b中的下标
%             max_num_list(max_index) =max_num_list(max_index) + max_num;
%         end
%         [max_num, max_index] = max(max_num_list); % 出现次数最多的元素出现的次数max_num，和该元素在b中的下标
%         D(id,5) = max_index;
%         merge_num = merge_num+1;
%     end
% end
%% 计算事件之间的空间相关系数
% 1、计算每个簇的中心点坐标
cluster_point = zeros(clusterLabel,2); % 存放每个簇的中心点坐标
for i = 1:clusterLabel
    cluster = D((D(:, 5) == i), :); % 第i个簇的所有点
    x = mean(cluster(:,2)); % 第i个簇的所有点的x坐标均值
    y = mean(cluster(:,3)); % 第i个簇的所有点的y坐标均值
    cluster_point(i,:) = [x y];
end
% 2、计算簇与簇的中心点坐标之间的距离（即簇之间的距离）矩阵
cluster_dist = zeros(clusterLabel,clusterLabel); % 存放每个簇的中心点坐标之间的距离矩阵
for i = 1 : clusterLabel
    for j = 1 : clusterLabel
        if i==j
            continue;
        end
        cluster_dist(i,j) = sqrt((cluster_point(i, 1) - cluster_point(j,1))^2 + (cluster_point(i, 2) - cluster_point(j,2))^2);   %计算i和j的空间距离
    end
end
% 让对角线的值等于所在行的倒数第二个最小值（即除了0之外的最小值）
for i = 1: length(cluster_dist)
     cluster_dist(i,i) =  min(cluster_dist(find(cluster_dist-min(cluster_dist))));
end
% 对距离矩阵的倒数进行归一化:距离矩阵的每个值都除以所在行的值之和
cluster_dist_sum = sum(cluster_dist); % cluster_dist的每行的累加值
for i = 1: length(cluster_dist)
    cluster_dist(:,i) =  cluster_dist(:,i)/cluster_dist_sum(i);    
end
% 3、计算簇的归一化密度，即有效相对密度 (归一化 X'=(X-MIN)/(MAX-MIN)) 
cluster_density = [clusterDensityListA;clusterDensityListB;clusterDensityListC]'; % 所有簇的混合密度
for i = 1:  size(cluster_density,2)
    cluster_density(:,i) = (cluster_density(:,i)-min(cluster_density(:,i))) / (max(cluster_density(:,i))-min(cluster_density(:,i)));
end
% 4、计算事件之间的空间相关系数
X = cluster_density(:,1);
Y = cluster_density(:,2);
Z = cluster_density(:,3);
W = cluster_dist;
% 计算A和B类型事件的相关系数
fenzi = 1/clusterLabel * (X - W * X)'* (Y-W*Y );
fenmu = sqrt(1/clusterLabel * (X - W * X)' * (X - W * X)) * sqrt(1/clusterLabel *(Y-W*Y )' * (Y-W*Y ));
sAB = fenzi/fenmu;
% 计算B和C类型事件的相关系数
fenzi = 1/clusterLabel * (Y - W * Y)'* (Z-W*Z );
fenmu = sqrt(1/clusterLabel * (Y - W * Y)' * (Y - W * Y)) * sqrt(1/clusterLabel *(Z-W*Z )' * (Z-W*Z ));
sBC = fenzi/fenmu;
% 计算A和C类型事件的相关系数
fenzi = 1/clusterLabel * (X - W * X)'* (Z-W*Z );
fenmu = sqrt(1/clusterLabel * (X - W * X)' * (X - W * X)) * sqrt(1/clusterLabel *(Z-W*Z )' * (Z-W*Z ));
sAC = fenzi/fenmu;

%% 计算簇内两种类型点的比例
% totalNumA = nnz(D(:,4)==0); %数据集中A类型点的数目
% totalNumB = nnz(D(:,4)==1); %数据集中A类型点的数目
% ratioList = [];     %簇中A类型点数量/B类型点数量列表
% % clusterNumList = [];    %簇中所有点的数目
% for i=1:clusterLabel
%     cluster = D((D(:, 5) == i), :);     %第i个簇
% %     clusterNumList = [clusterNumList;length(cluster(:,1))];
%     cluster1 = cluster((cluster(:,4) == 0),:);  %第i个簇中的A类型
%     cluster2 = cluster((cluster(:,4) == 1),:);  %第i个簇中的B类型
% 
%     ratio = (length(cluster1)/totalNumA) / (length(cluster2)/totalNumB);
%     ratioList = [ratioList; [i,ratio]];
% end

%% 绘图论文其他的插图
% 绘制柱状图
figure();
% temp = [clusterDensityListA', clusterDensityListB'];
% bar(temp(1:4,:),'stacked');
temp = [clusterDensityListA', clusterDensityListB',clusterDensityListC'];
bar(temp(1:6,:),'stacked');

% ids =1:6;
% [AX,H1,H2]=plotyy(ids,clusterDensityList',ids,clusterRelativeDensityList','bar','plot');
% set(AX(1),'ylim',[0 15]);
% set(AX(2),'ylim',[0 2])

% showlegend=[];
% % 1、绘制背景的圆
% Circles = [8,0;0,0;1,0]; % 圆心坐标
% Radius = [8; 4]; % 半径
% figure();
% showlegend = [showlegend,ShowCircle(Circles,Radius)];
% hold on;
% 2、绘制每个簇的第一个核心点位置
% colorTable = ['r','y','m'];
% for i = 1 : 1 % length(clusterCoreList)
%      showlegend = [showlegend, scatter(D(clusterCoreList(i), 2), D(clusterCoreList(i), 3),400,'.',colorTable(i))];
%      hold on;
% end

% 3、绘制带图例的二维聚类结果
% figure();
showlegend=[];
merge_num = 0;
real_clusterLabel = clusterLabel - merge_num;
showlegend = [showlegend, ShowCluster2Dimension(D,real_clusterLabel)];
% showlegend = [showlegend, ShowCluster(D,real_clusterLabel)];
% 4、绘制图例
legend(showlegend,'Noise','Cluster1','Cluster2','Cluster3','Cluster4','Cluster5','Cluster6','Cluster7','Cluster8','Cluster9','Cluster10','Cluster11','Cluster12','Cluster13','Cluster14','Cluster15','Cluster16','Cluster17','Cluster18','Cluster19','Cluster20','Cluster21','Cluster22','Cluster23','Cluster24','Cluster25');
% axis([-20,20,-15,15,-0.5,1.5]);
axis([-8,8,-8,8]);
set(gcf,'WindowStyle','normal');
grid on;    % 画出网格
xlabel('x');
ylabel('y');
zlabel('t');
% axis equal;
% set(gca,'YTicklabel',[] ); % 不显示坐标轴刻度数字
% set(gca,'ytick',[]); % 不显示坐标轴刻度
% set(gca,'xtick',[]);
% axis off;
% box on; % 加边框
% ShowCluster(D,real_clusterLabel);
%  ShowClusterXZ(D);
%% 绘制聚类结果
% figure();
% merge_num = 0;
% real_clusterLabel = clusterLabel - merge_num;
% ShowCluster(D,real_clusterLabel);
% ShowCluster2Dimension(D,real_clusterLabel);
% axis([-12,12,-9,9]);
%  axis equal;
%% 对所有簇内所有点的k近邻距离按照升序排序
% start = 1;
% stop = 0;
% prex = 0;
% prey = 0;
% for i = 1:clusterLabel
%     stop = stop+sum(labelList==i);
%     temp = labelKList(start:stop);
%     y = [prey];
%     tempID = [];
%     [y2,tempID] = sort(temp);     %对第i个点距离所有点的距离排序
%     y = [y,y2];
%     x = [prex,start:stop];
%     plot(x,y);
%     hold on;
%     start = stop+1;
%     prex = x(length(x));
%     prey = y(length(y));
% end
%% 对所有簇内的所有点的k近邻距离按照被标记（加入簇）的顺序进行绘图
% a = 1:length(labelList);
% plot(a,labelKList,'k');
% hold on;
% x=[];
% y=[];
% for i = 1 : length(labelFirstList)
%     x = [x,a(labelFirstList(i))];
%         y = [y,labelKList(labelFirstList(i))];
%     scatter(a(labelFirstList(i)),labelKList(labelFirstList(i)),['k', 'd']);
%     hold on;
% end
% plot(x,y,'r--');
% axis([0,1114,0,16]);
