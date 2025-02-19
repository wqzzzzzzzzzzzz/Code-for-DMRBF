clear all,close all ,clc
% 生成随机节点坐标
rng(1); % 设置随机种子以保持一致性
numNodes = 16;
nodeRadius = 0.3; % 节点半径

% 生成节点坐标
x = (rand(1, numNodes) - 0.5) * 2; % 在(-1, 1)范围内随机生成x坐标
y = (rand(1, numNodes) - 0.5) * 2; % 在(-1, 1)范围内随机生成y坐标

% 缩放节点坐标
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
% 绘制节点和标号
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
% 生成对称的邻接矩阵
adjMatrix = zeros(numNodes);
for i = 1:numNodes
    numConnections = randi([2, 6]);  % 随机生成连接数，范围在3到5之间
    connectedNodes = randperm(numNodes, numConnections);  % 随机选择连接节点，包括自身节点
    adjMatrix(i, connectedNodes) = 1; % 将连接节点标记为1
end
adjMatrix = adjMatrix | adjMatrix'; % 生成对称邻接矩阵

% 将对角线元素设为1，表示节点自连接
adjMatrix = adjMatrix | eye(numNodes);

% 确保每一行的元素不小于2且不大于6
for i = 1:numNodes
    rowSum = sum(adjMatrix(i, :));
    if rowSum < 2
        while rowSum < 2
            colIdx = find(adjMatrix(i, :) == 0); % 找到未连接的节点
            nodeIdx = randsample(colIdx, 2 - rowSum); % 随机选择节点进行连接
            adjMatrix(i, nodeIdx) = 1;
            adjMatrix(nodeIdx, i) = 1;
            rowSum = sum(adjMatrix(i, :));
        end
    elseif rowSum > 6
        while rowSum > 6
            colIdx = find(adjMatrix(i, :) == 1); % 找到已连接的节点
            nodeIdx = randsample(colIdx, rowSum - 5); % 随机选择节点断开连接
            adjMatrix(i, nodeIdx) = 0;
            adjMatrix(nodeIdx, i) = 0;
            rowSum = sum(adjMatrix(i, :));
        end
    end
end

% 将对角线元素设为1
adjMatrix = adjMatrix | eye(numNodes);

% 绘制连接线
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

% 设置图形属性
axis equal;
grid on;
xlabel('X-coordinate');
ylabel('Y-coordinate');
title('Complex Network Topology');
xlim([-0.4,0.3]);
ylim([-0.29,0.29]);

% 将每一行归一化
% 注意：在邻接矩阵的上下文中，归一化可能不是很有意义，因为元素通常是0或1
% 但如果你有特殊需求，可以这样做
normalizedAdjacencyMatrix = bsxfun(@rdivide, adjMatrix, sum(adjMatrix, 2));
% 由于可能存在全0行，这些行在归一化后会导致NaN，需要处理
% 这里我们可以选择将这些NaN替换为0或其他值
normalizedAdjacencyMatrix(isnan(normalizedAdjacencyMatrix)) = 0;
% 显示生成的邻接矩阵和归一化后的矩阵
disp('随机邻接矩阵:');
disp(adjMatrix);
disp('归一化后的邻接矩阵:');
disp(normalizedAdjacencyMatrix);

save('adjacencyMatrix.mat', 'normalizedAdjacencyMatrix');