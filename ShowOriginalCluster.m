function [ ] = ShowOriginalCluster( D )
%% ���������� ��ʾԭʼ���ݾ���Ч��
%   D�����ݼ�

%% ��������
figure();

for j = 1 : length(D(:,1))
    if D(j,4) == 0    %�����j������A���͵�
        scatter3(D(j, 2), D(j, 3),D(j, 1),15,['.', 'r']);
    else    %�����j������B���͵�
        scatter3(D(j, 2), D(j, 3),D(j, 1),15,['.', 'b']);
    end
    hold on;
end

% clusterA = D((D(:, 4) == 0), :);
% scatter3(clusterA(:, 2),clusterA(:,3),clusterA(:, 1), 'r','.');  
% hold on;
% clusterB = D((D(:, 4) == 1), :);
% scatter3(clusterB(:, 2),clusterB(:,3),clusterB(:, 1), 'k','*');  

% scatter3(D(:, 2),D(:,3), D(:, 1), 'k','.');  

grid on;    % ��������
% axis([0,100,0,100,0,100]);
end

