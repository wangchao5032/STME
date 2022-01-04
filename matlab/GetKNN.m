function [ distK_all, distKUpID,distKList] = GetKNN( D,k,t )
%% ����˵�� �������е�������K���ھ�
%   D:����
%   k����k���ھ�
%   t:ʱ�䴰��
%   distK_all:ÿ����������k���ھ�
%% ��������
distK_all = {}; %���е��k���ھ�

% 2019-11--------------------------
types_num = max(D(:, 4)) + 1;  % ���͵ı�Ŵ�0��ʼ��Ϊ0 1 2����������Ҫ���һ��1����һ����3������  
% distKList ��һ��1*n���������types_num*n�ľ���
% ��һ��Ϊ���������͵�k���룬��2����4��Ϊֻ����A���ͣ�B����C���͵�k����
distKList = zeros(types_num + 1, length(D(:, 1)));   
% --------------------------


for i=1 : (length(D(:,1)))
    temp = inf * ones(1, length(D(:, 1)));
    tempID = [];
    dist = [];
    for j = 1 : (length(D(:,1)))
        if i~=j && (abs(D(i,1)-D(j,1)) <= t)    %i������j��������ʱ����С�ڵ���t
             dist_i_j = sqrt((D(i, 2) - D(j,2))^2 + (D(i, 3) - D(j,3))^2);   %����i��j�Ŀռ����
        else
            dist_i_j = +inf;    %i��j�Ŀռ����Ϊ������
        end
        temp(j) = dist_i_j;
    end
    [dist,tempID] = sort(temp);     %�Ե�i����������е�ľ�������
    distK_all{i} = tempID(1:k);   %������������k���ھӵ�ID
    distKList(1, i) = dist(k);     %��i�����k������ľ���
    
    % 2019-11 ----------------------------------
    % ���㿼�����͵�k���룬�洢��distKList�ĺ�3��
    for q = 0:types_num-1
        temp_typet = temp(D(:, 4) == q);
        dist_t = mink(temp_typet, k);   % mink() �������Լ���temp_typet�����е�kС��ֵ
        if length(dist_t) == k
            distKList(q+2, i) = dist_t(k); 
        else
            distKList(q+2, i) = Inf;
        end
    end
    % ----------------------------------
    
end
[distKUp, distKUpID] = sort(distKList(1, :)); %�����е��k������ľ����������

end

