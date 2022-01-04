function [ d_result ] = Ripley_L( D )
%% ����˵��������L'(d)ȡ��Сֵ��Ӧ�ľ���d
%%   ��ȡ���ݣ������о������
xMin = min(D(:,2));
xMax = max(D(:,2));
yMin = min(D(:,3));
yMax = max(D(:,3));
A=(xMax-xMin)*(yMax-yMin); %�о������
%% ����L(d)
n = length(D(:,1));
dist_i_j_List = zeros(n,n); %��i����͵�j����ľ��뼯��
LdList=[]; %L(d)����
q = 50; %

for i = 1:n
    for j = 1:n
        dist_i_j = sqrt((D(i, 2) - D(j,2))^2 + (D(i, 3) - D(j,3))^2);   %����i��j�Ŀռ����
        dist_i_j_List(i,j) = dist_i_j;
        %dist_i_j_List=[dist_i_j_List,dist_i_j];
    end
 end


for d=0 : q
    total =0;
    for i = 1:n
        for j = 1:n
             dist_i_j = dist_i_j_List(i,j);
             %dist_i_j = sqrt((D(i, 2) - D(j,2))^2 + (D(i, 3) - D(j,3))^2);   %����i��j�Ŀռ����
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

%% ��L(d)�ĵ��� 
derList=[]; %L(d)��������
for i=1: length(LdList(1,:))-1
  der = LdList(i+1)-LdList(i);
  derList = [derList,der];
end
figure();
plot(derList);
%L'(d)����Сֵ��L'(d)��Сֵ���ڵ�λ��
[minDer,minDerLocate] = min(derList);
 %L'(d)ȡ��Сֵʱ��Ӧ�ľ���d
 d_result = dist_i_j_List(minDerLocate);
end

