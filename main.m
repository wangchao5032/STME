clc;clear;
format long ;   %[��ΰ������]����С��λ�����ӣ������ж�
%% ��excel�����ж�ȡ����
%D:ģ�����ݣ� �ֶΣ�ʱ�䣬x,y,���ͣ�0����A���ͣ�1����B���ͣ����ر�ţ�Ĭ��Ϊ0�������ܻ��е��������
addpath('����������ͼ') 
% D = xlsread('simple_data.xlsx');
D = xlsread('ģ������11\data.xlsx');
 
% 2019-11 ------------------------------
% �����µ������ģ������
% D3 = []; D2 = D3; D1 = D2;
% D1 = generate_ring_data(0, 0,  0, 4, 0.3, 0.5, 0);
% D2 = generate_ring_data(6, 8,  0, 10, 0.9, 0.5, 1);
% D3 = generate_ring_data(16, -8, 0, 18, 2,   0.5, 2);
% D = [D1; D2; D3];
% --------------------------------------------

% ����ԭʼ����
% ShowOriginalCluster(D);
% ShowOriginalCluster2Dimension(D);
%% �������
%d = Ripley_L(D); %�����С��ֱ�� ��Ϊ���㷨�����̶��ص�ֱ�������Բ��ü���
k = 20;
tWindow = 5;
kt = 10; % ceil(k*0.5);
MinPts = 10; % ceil(k*0.5);
distK_sigma_multi = 10; % �����������Ϳ��Էſ�k�������ֵ
distK_sigma_times = 3; 
%% ����ÿ���������K���ھ� �� ���е��k���ھӵ����������ľ���
[knnList, distKUpList,distKList] = GetKNN(D,k,tWindow);  
distKList_typet = distKList(2:end, :);  % �������͵�k����
distKList = distKList(1, :); % ���������͵�k����
%% ����ÿ����Ĺ�������
RSTNList = {};   %���������б�
for i = 1 : length(D(:,1))
    RSTN = [];      %��������
    for j = 1 : k     %����k���ھ�
        compareObject = knnList{i}(j);  %k���ھ��еĵ�j��
        if i == compareObject
            continue;
        end
        %�����j��k�ھӺ�i�Ĺ������ڸ�������kt����iҲ�ǵ�j��k�ھӵ�K�ھ�
        if length(intersect(knnList{i},knnList{compareObject})) >= kt && ismember(i,knnList{compareObject})==1  
             RSTN = [RSTN,compareObject];    %�����ھӼ��빲�����ڶ���
        end
    end
    RSTNList{i} = RSTN;
end
%% �������
clusterLabel = 0;   %��ʼ���ر��
clusterDensityList = [];    %����Ч�ܶ��б�
clusterDensityListA = [];   %����A���͵����Ч�ܶ��б�
clusterDensityListB = [];   %����B���͵����Ч�ܶ��б�
clusterDensityListC = [];   %����C���͵����Ч�ܶ��б�
clusterNumList = []; % ���е������ 

clusterCoreList = []; % �صĵ�һ�����ĵ��б�

labelKList = []; %���ձ�ǵ�˳���¼����ǵ��k����
labelList = []; %���ձ�ǵ�˳���¼����ǵĴغ�
labelFirstList = []; %���е�һ�α���ǵĵ� �� ���������еڼ�������ǵĵ�
for ii = 1 : length(D(:,1))     %���յ�k���ھӾ�����������ö�٣�Ҳ���ǵ�i����һ����Ŀǰ���е���k����Χ��С�ġ�
    i = distKUpList(ii);    %��ii�������е��k���ھӾ�����������򡱵ı��
    if D(i,5) ~= 0  %��i�����Ѿ���������
        continue;
    end
    Y = knnList{i};     % ��i�����k���ھ�
    X = RSTNList{i};    % ��i����Ĺ����ھ�/ֱ�ӿɴ��
    RSTN_label = D(X,5); % ��i����Ĺ������ڵı�ǩ
    RSTN_free_pts = nnz((RSTN_label == 0 | RSTN_label == -1)); % ��i����Ĺ��������б�ǩΪ0��-1�ĵ�
    if length(X) < MinPts  % �����i���������������Ǻ��ĵ�
        D(i, 5) = -1; 
    elseif RSTN_free_pts < MinPts
        RSTN_most_label = mode(RSTN_label);
        D(i,5) = RSTN_most_label;
    else % �����i�����Ǻ��ĵ㣬��ʼһ���µĴ�
        clusterLabel = clusterLabel + 1;    % �ر��+1
        clusterValidList = [];  % ������㵱ǰ����Ч����ļ��ϣ���ŵ�ǰ�ص����к��ĵ�
        clusterCoreList = [clusterCoreList,i]; % ��ǰ�صĵ�һ�����ĵ�
        labelFirstList = [labelFirstList,length(labelList)+1]; % ������ǰ�صĵ�һ�����ĵ�ʱ���Ѿ���length(labelList)���㱻��ǹ��ˣ����Ե�ǰ�صĵ�һ�����ĵ㼴���ǵ�length(labelList)+1������ǵĵ�
         % �����ʼ����
        queue = [i, RSTNList{i}];               %����i���㼰��RSTN�������
        distK_buf = distKList_typet(:, queue); % ���������е�Ŀ������͵�k����
        type_buf = D(queue, 4); % �ص���Ч�㣺��ʼ���������е�+����������еĵ��еĺ��ĵ�
        % �����ʼ���������е�ľ�ֵ�ͱ�׼��
        distK_mu = mean(distK_buf, 2);                                          %���㡰���������е�Ŀ������͵�k���롱�е�� ��ֵ
        distK_sigma_init = sqrt(var(distK_buf, 1, 2));                               %���㡰���������е�Ŀ������͵�k���롱�е�� ��׼��
        
        cluster1_num = 0;
        % ��������
        while ~isempty(queue) 
            %  Ϊ�˻�ͼ����ʱ�ӵģ���ʽ��ʱ�ǵ�ɾ��
%             cluster1_num = cluster1_num+1;
%             if (clusterLabel==1 && cluster1_num==80) 
%                 break;
%             end
            % -----------------1��ͳ�Ƶ�ǰ�����и����͵����Ŀ
            type_hist = [nnz(type_buf == 0) nnz(type_buf == 1) nnz(type_buf == 2)]; % ͳ�Ƶ�ǰ��������Ч�����������
            type_valid = (type_hist >= MinPts);                   % ��ǰ���е���Ч���У����ĳ�����͵ĵ���Ŀ����MinPts������Ϊ��Щ��������Ӱ������
            [~, dominate_type] = max(type_hist);                                    % �ҳ���ǰ���е���Ч���У�Ӱ�����������ͣ������к��ԣ�
            distK_sigma = distK_sigma_init;
            distK_sigma(dominate_type) = distK_sigma(dominate_type) * distK_sigma_multi;          % Ӱ�����������Ϳ��Էſ�k�������ֵ��Ϊ��ʹģ�����ݱ߽��ܶȽϵ͵ĵ�����������������к��ԣ� 
            % -----------------2��ȡ����ǰ�㲢���
            ptCurrent = queue(1);   % ȡ�������еĵ�һ������Ϊ��ǰ��
            queue(1) = []; % �ѵ�ǰ���Ƴ�����
            if D(ptCurrent, 5) > 0  % �����ǰ�㱻��ǹ��������ֻ��������Ϊ�������ڱ���������ˡ�
                continue;
            end 
            % �����еĵ����������:
            % (1)��ʼ���γ�ʱ�����ĵ�1���㣬һ����û�б����Ҳ���������ĺ��ĵ㡣
            % (2)��ʼ���γ�ʱ�����ĵ�1����Ĺ������ڣ���һ���Ǻ��ĵ㣬���ܱ���ǹ������ں��Ľ�һ��ɸѡ�Ƿ��Ǻ��ĵ㣨if length(X2) >= MinPts����������ٽ���k�ھӼ������
            % (3)���ڼ�������еĵ㣬һ����û�б���ǹ��ĵ����������һ���Ǻ��ĵ㡣���ں��Ľ�һ��ɸѡ�Ƿ��Ǻ��ĵ㣨if length(X2) >= MinPts����������ٽ���k�ھӼ������
            D(ptCurrent, 5) = clusterLabel;     
            labelKList = [labelKList,distKList(ptCurrent)];
            labelList = [labelList,clusterLabel]; 
            % -----------------3��������ǰ���k�ھӣ��������
            Y2 = knnList{ptCurrent};     % �����е�ǰ���k���ھ�
            X2 = RSTNList{ptCurrent};    % �����е�ǰ��Ĺ����ھ�/ֱ�ӿɴ��
            if length(X2) >= MinPts  % �����ǰ���Ǻ��ĵ�
                type_buf = [type_buf; D(ptCurrent, 4)]; % �ѵ�ǰ����뵱ǰ���е���Ч���б���������bug���еĵ�ᱻ�ظ����룡������
                clusterValidList = [clusterValidList,ptCurrent];             
                for j = 1 : length(Y2)   % ������ǰ���k�ھ�
                    distK_diff = abs(distKList_typet(:, Y2(j)) - distK_mu);   % ���㵱ǰ���j���ھӵĿ������͵�k���루1*3�� �� ��ʼ�ؿ������͵�k�����ֵ��1*3�� �Ĳ�ֵ��1*3��                 
                    % �����ǰ��ĵ�j���ھ�(1)��������û�б���ǹ�;(2)���ڶ�����;
                    % (3)��ֻ������Ӱ�����������£���ǰ���j���ھӵĿ������͵�k���루1*2�� �� ��ʼ�ؿ������͵�k�����ֵ��1*2���Ĳ�ֵdistK_diff(type_valid)��1*2��������ֵ�ڳ�ʼ�ص�3����׼������
                    % D(Y2(j),5) == -1 || 
                    if ( D(Y2(j),5) == -1 ||D(Y2(j),5) == 0 )  && isempty(find(queue == Y2(j), 1)) &&  all(distK_diff(type_valid) < distK_sigma(type_valid) * distK_sigma_times) 
                        queue = [queue, Y2(j)] ;    % ����ǰ��ĵ�j���ھӼ������
                    end
                end
            end
        end
        % -----------------4�����㵱ǰ�صĵ���
        cluster = D((D(:, 5) == clusterLabel), :);     %��ǰ��
        clusterNumList = [clusterNumList, length(find(D(:, 5) == clusterLabel))]; % ��ǰ�صĵ������
        clusterNumA = nnz(cluster(:,4)==0); %��ǰ����A���͵����Ŀ
        clusterNumB = nnz(cluster(:,4)==1); %��ǰ����B���͵����Ŀ
        clusterNumC = nnz(cluster(:,4)==2); %��ǰ����C���͵����Ŀ
        
        % -----------------5�����㵱ǰ����Ч�ܶ�
        %������ά��͹��
%           [K,V] = convhull(cluster(:,2),cluster(:,3),cluster(:,1));   %��ǰ�ص���С͹��
        %�����ά��͹��
        [K,V] = convhull(cluster(:,2),cluster(:,3));   %��ǰ�ص���С͹��
        
        clusterDensity =  length(cluster) / V ;    %��ǰ�ص��ܶ�=��ǰ�ص�����/��ǰ�ص����
        clusterDensityList = [clusterDensityList;clusterDensity];
        clusterDensityA = clusterNumA /V; %��ǰ�ص�A���͵���ܶ� = ��ǰ����A���͵������/��ǰ�ص����    
        clusterDensityListA = [clusterDensityListA,clusterDensityA];
        clusterDensityB = clusterNumB /V; %��ǰ�ص�B���͵���ܶ� = ��ǰ����B���͵������/��ǰ�ص����
        clusterDensityListB = [clusterDensityListB,clusterDensityB];
        clusterDensityC = clusterNumC /V; %��ǰ�ص�C���͵���ܶ� = ��ǰ����C���͵������/��ǰ�ص����
        clusterDensityListC = [clusterDensityListC,clusterDensityC];
       
%         clusterA = []; %���浱ǰ�ص����к��ĵ��е�A����
%         clusterB = []; %���浱ǰ�ص����к��ĵ��е�B����
%         cluster = []; %���浱ǰ�ص����к��ĵ�
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
%% ����С�Ĵغϲ��������Ĵ�
% merge_num = 0;
% for i = 1:length(clusterNumList)
%     if (clusterNumList(i) < MinPts) % ������е��������MinPts��С
%         [id,id_] = find(D(:, 5) == i); 
%         max_num_list = zeros(1,length(clusterNumList)); % ��Ÿôص����е�Ĺ��������У������غų��ֵ�Ƶ��
%         for j = 1:length(id) % �����ô��е����е㣬ͳ�ƴ��д����������Ĵ�
%             RSTN_j = RSTNList{id(j)}; % �صĵ�j����Ĺ�������
%             labels = D(RSTN_j,5); % �صĵ�j����Ĺ������ڵ����дغ�
%             % ���ִ������Ĵغ�
%             temp = 1:max(labels); % �����п��ܴ��ڵ������ļ���
%             c = histc(labels,temp); % �����ֳ��ֵĴ���
%             [max_num, max_index] = max(c); % ���ִ�������Ԫ�س��ֵĴ���max_num���͸�Ԫ����b�е��±�
%             max_num_list(max_index) =max_num_list(max_index) + max_num;
%         end
%         [max_num, max_index] = max(max_num_list); % ���ִ�������Ԫ�س��ֵĴ���max_num���͸�Ԫ����b�е��±�
%         D(id,5) = max_index;
%         merge_num = merge_num+1;
%     end
% end
%% ����ص���Ч����ܶ�(��һ�� X'=(X-MIN)/(MAX-MIN))
% 
clusterRelativeDensityList =[]; %�ص���Ч����ܶ�
clusterRelativeDensityListA = []; %����A���͵����Ч����ܶ�
clusterRelativeDensityListB = []; %����B���͵����Ч����ܶ�
% clusterRelativeDensityListC = []; %����C���͵����Ч����ܶ�
for i = 1 : clusterLabel
    clusterRelativeDensityList = [clusterRelativeDensityList;(clusterDensityList(i)-min(clusterDensityList))/(max(clusterDensityList)-min(clusterDensityList))];
    clusterRelativeDensityListA = [clusterRelativeDensityListA;(clusterDensityListA(i)-min(clusterDensityListA))/(max(clusterDensityListA)-min(clusterDensityListA))];
    clusterRelativeDensityListB = [clusterRelativeDensityListB;(clusterDensityListB(i)-min(clusterDensityListB))/(max(clusterDensityListB)-min(clusterDensityListB))];
%     clusterRelativeDensityListC = [clusterRelativeDensityListC;(clusterDensityListC(i)-min(clusterDensityListC))/(max(clusterDensityListC)-min(clusterDensityListC))];
end
%% ���㲻ͬ���͵�֮������ϵ�� 
%{
sAB = corrcoef([clusterRelativeDensityListA,clusterRelativeDensityListB]);
corAB = sAB(1,2);
% sAC = corrcoef([clusterRelativeDensityListA,clusterRelativeDensityListC]);
% corAC = sAC(1,2);
% sBC = corrcoef([clusterRelativeDensityListB,clusterRelativeDensityListC]);
% corBC = sBC(1,2);
%}
%% ��������������͵�ı���
% totalNumA = nnz(D(:,4)==0); %���ݼ���A���͵����Ŀ
% totalNumB = nnz(D(:,4)==1); %���ݼ���A���͵����Ŀ
% ratioList = [];     %����A���͵�����/B���͵������б�
% % clusterNumList = [];    %�������е����Ŀ
% for i=1:clusterLabel
%     cluster = D((D(:, 5) == i), :);     %��i����
% %     clusterNumList = [clusterNumList;length(cluster(:,1))];
%     cluster1 = cluster((cluster(:,4) == 0),:);  %��i�����е�A����
%     cluster2 = cluster((cluster(:,4) == 1),:);  %��i�����е�B����
% 
%     ratio = (length(cluster1)/totalNumA) / (length(cluster2)/totalNumB);
%     ratioList = [ratioList; [i,ratio]];
% end
%% ���� ���������
% clusterResult = D(:,5);
% clusterTrue = D(:,6);
% falseNum = length(find((clusterResult-clusterTrue)~=0)) ;
%% ��ͼ���������Ĳ�ͼ
% ������״ͼ
figure();
temp = [clusterDensityListA', clusterDensityListB',clusterDensityListC'];
bar(temp(1:6,:),'stacked');

% ids =1:6;
% [AX,H1,H2]=plotyy(ids,clusterDensityList',ids,clusterRelativeDensityList','bar','plot');
% set(AX(1),'ylim',[0 15]);
% set(AX(2),'ylim',[0 2])

% showlegend=[];
% % 1�����Ʊ�����Բ
% Circles = [8,0;0,0;1,0]; % Բ������
% Radius = [8; 4]; % �뾶
% figure();
% showlegend = [showlegend,ShowCircle(Circles,Radius)];
% hold on;
% 2������ÿ���صĵ�һ�����ĵ�λ��
% colorTable = ['r','y','m'];
% for i = 1 : 1 % length(clusterCoreList)
%      showlegend = [showlegend, scatter(D(clusterCoreList(i), 2), D(clusterCoreList(i), 3),400,'.',colorTable(i))];
%      hold on;
% end

% 3�����ƴ�ͼ���Ķ�ά������
% figure();
showlegend=[];
merge_num = 0;
real_clusterLabel = clusterLabel - merge_num;
% showlegend = [showlegend, ShowCluster2Dimension(D,real_clusterLabel)];
showlegend = [showlegend, ShowCluster(D,real_clusterLabel)];
% 4������ͼ��
legend(showlegend,'Noise','Cluster1','Cluster2','Cluster3','Cluster4','Cluster5','Cluster6');
axis([-20,20,-15,15,-0.5,1.5]);
% axis([-20,20,-15,15]);
set(gcf,'WindowStyle','normal');
grid on;    % ��������
xlabel('x');
ylabel('y');
zlabel('t');
% axis equal;
% set(gca,'YTicklabel',[] ); % ����ʾ������̶�����
% set(gca,'ytick',[]); % ����ʾ������̶�
% set(gca,'xtick',[]);
% axis off;
% box on; % �ӱ߿�
% ShowCluster(D,real_clusterLabel);
%  ShowClusterXZ(D);
%% ���ƾ�����
% figure();
% merge_num = 0;
% real_clusterLabel = clusterLabel - merge_num;
% ShowCluster(D,real_clusterLabel);
% ShowCluster2Dimension(D,real_clusterLabel);
% axis([-12,12,-9,9]);
%  axis equal;
%% �����д������е��k���ھ��밴����������
% start = 1;
% stop = 0;
% prex = 0;
% prey = 0;
% for i = 1:clusterLabel
%     stop = stop+sum(labelList==i);
%     temp = labelKList(start:stop);
%     y = [prey];
%     tempID = [];
%     [y2,tempID] = sort(temp);     %�Ե�i����������е�ľ�������
%     y = [y,y2];
%     x = [prex,start:stop];
%     plot(x,y);
%     hold on;
%     start = stop+1;
%     prex = x(length(x));
%     prey = y(length(y));
% end
%% �����д��ڵ����е��k���ھ��밴�ձ���ǣ�����أ���˳����л�ͼ
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
%% ͳ��ÿ�����ڻ�ϵ����Ŀ
% clusterNumList = [];
% for i = 1 : clusterLabel
%     clusterNumList = [clusterNumList,length(D(D(:,5)==i))];
% end
%% ����ÿ������AB�������͵�ı���
% ratioList = [];
% for i = 1 : clusterLabel
%     ratioList = [ratioList;clusterDensityListA(i)/clusterDensityListB(i)];
% end