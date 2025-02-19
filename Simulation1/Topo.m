clear all,close all ,clc
% ��������ڵ�����
rng(1); % ������������Ա���һ����
numNodes = 16;
nodeRadius = 0.3; % �ڵ�뾶

% ���ɽڵ�����
x = (rand(1, numNodes) - 0.5) * 2; % ��(-1, 1)��Χ���������x����
y = (rand(1, numNodes) - 0.5) * 2; % ��(-1, 1)��Χ���������y����

% ���Žڵ�����
x = x * nodeRadius;

x(1)=x(1)+0.04;
x(3)=x(3)-0.02;
x(8)=x(8)+0.07;
x(7)=x(7)+0.17;
x(12)=x(12)+0.06;
x(13)=x(13)-0.015;
x(5)=x(5)+0.04;
x(16)=x(16)+0.08;


y = y * nodeRadius;
y(14)=y(14)+0.05;
y(1)=y(1)+0.01;
y(3)=y(3)+0.2;
y(8)=y(8)-0.01;
y(7)=y(7)-0.005;
y(4)=y(4)+0.03;
y(16)=y(16)-0.12;
y(10)=y(10)+0.03;
% ���ƽڵ�ͱ��
plot(x, y, 'o', 'MarkerSize', 15, 'MarkerFaceColor', 'y');
hold on;
for i = 1:numNodes


    if i == 14
        text(x(i)-0.000, y(i), num2str(i), 'Color', 'k', 'FontSize', 12, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    elseif i==1
        text(x(i), y(i)+0.000, num2str(i), 'Color', 'k', 'FontSize', 12, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    elseif i==16
        text(x(i)-0.000, y(i)+0.000, num2str(i), 'Color', 'k', 'FontSize', 12, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    elseif i==18
        text(x(i)-0.000, y(i)+0.000, num2str(i), 'Color', 'k', 'FontSize', 12, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    elseif i==3
        text(x(i)-0.000, y(i)+0.000, num2str(i), 'Color', 'k', 'FontSize', 12, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    elseif i==10
        text(x(i)-0.000, y(i)+0.000, num2str(i), 'Color', 'k', 'FontSize', 12, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    elseif i==11
        text(x(i)-0.000, y(i)+0.000, num2str(i), 'Color', 'k', 'FontSize', 12, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    elseif i==13
        text(x(i)+0.000, y(i)+0.000, num2str(i), 'Color', 'k', 'FontSize', 12, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    elseif i==19
        text(x(i)-0.000, y(i)+0.000, num2str(i), 'Color', 'k', 'FontSize', 12, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    elseif i==5
        text(x(i)+0.000, y(i)+0.000, num2str(i), 'Color', 'k', 'FontSize', 12, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    elseif i==17
        text(x(i)-0.000, y(i)+0.000, num2str(i), 'Color', 'k', 'FontSize', 12, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    elseif i==6
        text(x(i), y(i)+0.000, num2str(i), 'Color', 'k', 'FontSize', 12, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    elseif i==12
        text(x(i)+0.000, y(i)+0.000, num2str(i), 'Color', 'k', 'FontSize', 12, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    elseif i==7
        text(x(i)+0.000, y(i)+0.000, num2str(i), 'Color', 'k', 'FontSize', 12, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    elseif i==2
        text(x(i)-0.000, y(i)+0.000, num2str(i), 'Color', 'k', 'FontSize', 12, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    elseif i==9
        text(x(i)-0.000, y(i)+0.000, num2str(i), 'Color', 'k', 'FontSize', 12, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    elseif i==8
        text(x(i)-0.000, y(i)+0.000, num2str(i), 'Color', 'k', 'FontSize', 12, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    else
        text(x(i)+0.000, y(i)+0.000, num2str(i), 'Color', 'k', 'FontSize', 12, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    end



end
% ���ɶԳƵ��ڽӾ���
adjMatrix = zeros(numNodes);
for i = 1:numNodes
    numConnections = randi([2, 6]);  % �����������������Χ��3��5֮��
    connectedNodes = randperm(numNodes, numConnections);  % ���ѡ�����ӽڵ㣬��������ڵ�
    adjMatrix(i, connectedNodes) = 1; % �����ӽڵ���Ϊ1
end
adjMatrix = adjMatrix | adjMatrix'; % ���ɶԳ��ڽӾ���

% ���Խ���Ԫ����Ϊ1����ʾ�ڵ�������
adjMatrix = adjMatrix | eye(numNodes);

% ȷ��ÿһ�е�Ԫ�ز�С��2�Ҳ�����6
for i = 1:numNodes
    rowSum = sum(adjMatrix(i, :));
    if rowSum < 2
        while rowSum < 2
            colIdx = find(adjMatrix(i, :) == 0); % �ҵ�δ���ӵĽڵ�
            nodeIdx = randsample(colIdx, 2 - rowSum); % ���ѡ��ڵ��������
            adjMatrix(i, nodeIdx) = 1;
            adjMatrix(nodeIdx, i) = 1;
            rowSum = sum(adjMatrix(i, :));
        end
    elseif rowSum > 6
        while rowSum > 6
            colIdx = find(adjMatrix(i, :) == 1); % �ҵ������ӵĽڵ�
            nodeIdx = randsample(colIdx, rowSum - 5); % ���ѡ��ڵ�Ͽ�����
            adjMatrix(i, nodeIdx) = 0;
            adjMatrix(nodeIdx, i) = 0;
            rowSum = sum(adjMatrix(i, :));
        end
    end
end

% ���Խ���Ԫ����Ϊ1
adjMatrix = adjMatrix | eye(numNodes);

% ����������
for i = 1:numNodes
    for j = i+1:numNodes
        if adjMatrix(i, j) == 1
            x_start = x(i);
            y_start = y(i);
            x_end = x(j);
            y_end = y(j);

            x_line = [x_start, x_end];
            y_line = [y_start, y_end];
            line(x_line, y_line, 'Color', 'b','LineWidth',1.5);
        end
    end
end

% ����ͼ������
axis equal;
grid on;
xlabel('X-coordinate');
ylabel('Y-coordinate');
title('Complex Network Topology');
xlim([-0.4,0.3]);
ylim([-0.29,0.29]);

% ��ÿһ�й�һ��
% ע�⣺���ڽӾ�����������У���һ�����ܲ��Ǻ������壬��ΪԪ��ͨ����0��1
% ����������������󣬿���������
normalizedAdjacencyMatrix = bsxfun(@rdivide, adjMatrix, sum(adjMatrix, 2));
% ���ڿ��ܴ���ȫ0�У���Щ���ڹ�һ����ᵼ��NaN����Ҫ����
% �������ǿ���ѡ����ЩNaN�滻Ϊ0������ֵ
normalizedAdjacencyMatrix(isnan(normalizedAdjacencyMatrix)) = 0;
% ��ʾ���ɵ��ڽӾ���͹�һ����ľ���
disp('����ڽӾ���:');
disp(adjMatrix);
disp('��һ������ڽӾ���:');
disp(normalizedAdjacencyMatrix);

save('adjacencyMatrix.mat', 'normalizedAdjacencyMatrix');