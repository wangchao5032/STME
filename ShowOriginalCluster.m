function [ showlegend] = ShowOriginalCluster(D)
%% ���������� ��ʾԭʼ���ݾ���Ч��
%   D�����ݼ�

%% ��������
figure();
showlegend=[];% ��ʾͼ����ͼ��

clusterA = D((D(:, 4) == 0), :);
showlegend = [showlegend,scatter3(clusterA(:, 2),clusterA(:,3),clusterA(:, 1), 'r','.')];
hold on;
clusterB = D((D(:, 4) == 1), :);
showlegend = [showlegend,scatter3(clusterB(:, 2),clusterB(:,3),clusterB(:, 1), 'b','.')];
% hold on;
% clusterC = D((D(:, 4) == 2), :);
% scatter3(clusterC(:, 2),clusterC(:,3),clusterC(:, 1), 'b','*');  
% scatter3(D(:, 2),D(:,3), D(:, 1), 'k','.');  

grid on;    % ��������
end

