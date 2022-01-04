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
hold on;
clusterC = D((D(:, 4) == 2), :);
showlegend = [showlegend,scatter3(clusterC(:, 2),clusterC(:,3),clusterC(:, 1), 'g','.')];

grid on;    % ��������
end

