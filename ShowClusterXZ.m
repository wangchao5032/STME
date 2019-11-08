function [ ] = ShowClusterXZ( D )
%% ���������� ��ʾ����Ч��
%   D�����ݼ�
%% ��������
figure();
% colorTable = ['b','r','y','g','m','c','b','r','y','g','m','c','b','r','y','g','m','c','b','r','y','g','m','c','b','r','y','g','m','c','b','r','y','g','m','c','b','r','y','g','m','c','b','r','y','g','m','c','b','r','y','g','m','c','b','r','y','g','m','c','b','r','y','g','m','c','b','r','y','g','m','c','b','r','y','g','m','c','b','r','y','g','m','c','b','r','y','g','m','c'];
% shapeTable = ['o','+','*','x','s','d','^','v','>','<','p','h','o','+','*','x','s','d','^','v','>','<','p','h','o','+','*','x','s','d','^','v','>','<','p','h','o','+','*','x','s','d','^','v','>','<','p','h','o','+','*','x','s','d','^','v','>','<','p','h','o','+','*','x','s','d','^','v','>','<','p','h','o','+','*','x','s','d','^','v','>','<','p','h','o','+','*','x','s','d','^','v','>','<','p','h'];
colorTable = ['b','r','b','r','m','c','r','g','m','m','r','c','b','g','c','r','m','c','b','r','y','g','m','c','b','r','y','g','m','c','b','r','y','g','m','c','b','r','y','g','m','c','b','r','y','g','m','c','b','r','y','g','m','c','b','r','y','g','m','c','b','r','y','g','m','c','b','r','y','g','m','c','b','r','y','g','m','c','b','r','y','g','m','c','b','r','y','g','m','c'];
shapeTable = ['o','+','*','x','s','d','^','v','>','<','p','h','+','p','d','s','d','p','*','*','p','+','<','s','o','+','*','x','s','d','^','v','>','<','p','h','o','+','*','x','s','d','^','v','>','<','p','h','o','+','*','x','s','d','^','v','>','<','p','h','o','+','*','x','s','d','^','v','>','<','p','h','o','+','*','x','s','d','^','v','>','<','p','h','o','+','*','x','s','d','^','v','>','<','p','h'];
% �����㻭СԲ��
cluster = D((D(:, 5) == -1), :);
scatter(cluster(:,2),cluster(:,1),20,'.','k');  
hold on;

% �����Ѿ��ɹ�����ĵ㣬��ɫ����colorTable����״����shapeTable
clusterNum=max(D(:,5));

for i=1:clusterNum
    cluster = D((D(:, 5) == i), :);
    for j = 1 : length(cluster(:,1))
        if cluster(j,4) == 0    %�����j������A���͵�
            scatter3(cluster(j, 2), cluster(j, 1),cluster(j, 1),10,'+',colorTable(i));
        else    %�����j������B���͵�
            scatter3(cluster(j, 2), cluster(j, 1),cluster(j, 1),10,'<', colorTable(i));
        end
    end
    hold on;
end


% for i=1:clusterNum
%     cluster = D((D(:, 5) == i), :);     %��i����
%     scatter(cluster(:, 2),cluster(:, 1),shapeTable(i),colorTable(i));
%         hold on;
% end
% axis([0,100,0,100,0,100]);
grid on;    % ��������

end

