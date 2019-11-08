function [ d_result ] = Ripley_L( D )
%% 函数说明：计算L'(d)取最小值对应的距离d
%%   读取数据，计算研究区面积
xMin = min(D(:,2));
xMax = max(D(:,2));
yMin = min(D(:,3));
yMax = max(D(:,3));
A=(xMax-xMin)*(yMax-yMin); %研究区面积
%% 计算L(d)
n = length(D(:,1));
dist_i_j_List = zeros(n,n); %第i个点和第j个点的距离集合
LdList=[]; %L(d)集合
q = 50; %

for i = 1:n
    for j = 1:n
        dist_i_j = sqrt((D(i, 2) - D(j,2))^2 + (D(i, 3) - D(j,3))^2);   %计算i和j的空间距离
        dist_i_j_List(i,j) = dist_i_j;
        %dist_i_j_List=[dist_i_j_List,dist_i_j];
    end
 end


for d=0 : q
    total =0;
    for i = 1:n
        for j = 1:n
             dist_i_j = dist_i_j_List(i,j);
             %dist_i_j = sqrt((D(i, 2) - D(j,2))^2 + (D(i, 3) - D(j,3))^2);   %计算i和j的空间距离
             %dist_i_j_List=[dist_i_j_List,dist_i_j];
             if dist_i_j<=d
                 delta = 1;  
             else
                 delta = 0;
             end
             total =  total + delta/(n^2);
        end
    end
    Kd = A*total; %K(d)
    Ld = sqrt(Kd/pi)-d;  %L(d)
    LdList = [LdList,Ld];
end

figure();
plot(LdList); 

%% 求L(d)的导数 
derList=[]; %L(d)导数集合
for i=1: length(LdList(1,:))-1
  der = LdList(i+1)-LdList(i);
  derList = [derList,der];
end
figure();
plot(derList);
%L'(d)的最小值和L'(d)最小值所在的位置
[minDer,minDerLocate] = min(derList);
 %L'(d)取最小值时对应的距离d
 d_result = dist_i_j_List(minDerLocate);
end

