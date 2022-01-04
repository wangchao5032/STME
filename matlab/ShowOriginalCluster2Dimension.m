function [showlegend ] = ShowOriginalCluster2Dimension( D )
%% ���������� ��ʾԭʼ���ݾ���Ч��
%   D�����ݼ�

%% ��������
figure();
showlegend=[];% ��ʾͼ����ͼ��

clusterA = D((D(:, 4) == 0), :);
showlegend = [showlegend,scatter(clusterA(:, 2),clusterA(:,3), 'r','.')];  
hold on;
clusterB = D((D(:, 4) == 1), :);
showlegend = [showlegend,scatter(clusterB(:, 2),clusterB(:,3), 'b','.')];  
hold on;
clusterC = D((D(:, 4) == 2), :);
showlegend = [showlegend,scatter(clusterC(:, 2),clusterC(:,3), 'k','.')];  

% scatter(D(:, 2),D(:,3), 'k','.');  
grid on;    % ��������

end

