function [ ] = ShowOriginalCluster2Dimension( D )
%% ���������� ��ʾԭʼ���ݾ���Ч��
%   D�����ݼ�

%% ��������
figure();

for j = 1 : length(D(:,1))
    if D(j,4) == 0    %�����j������A���͵�
        scatter(D(j, 2), D(j, 3),['.', 'r']);
    else    %�����j������B���͵�
        scatter(D(j, 2), D(j, 3),['.', 'b']);
    end
    hold on;
end

% clusterA = D((D(:, 4) == 0), :);
% scatter(clusterA(:, 2),clusterA(:,3), 'r','.');  
% hold on;
% clusterB = D((D(:, 4) == 1), :);
% scatter(clusterB(:, 2),clusterB(:,3), 'b','.');  

% scatter(D(:, 2),D(:,3), 'k','.');  
grid on;    % ��������
% axis([0,100,0,100]);
end

