clc
clear
close all

%% 数据准备
% 优化后数据
load('GA_results.mat');
a_ = results.population;
r_optimized = a_(end,:);

% 均匀阵列数据
Na = 16; % 辐条数量
Nm = 8;  % 每个辐条的阵元数
M = Na * Nm; % 总阵元数量

radiusMin = 0.1; % 最小半径
radiusMax = 0.7; % 最大半径
frequency = 1000; % 入射波频率
c = 343; % 声速

% 计算波数
k = 2*pi*frequency/c;

% 初始化阵元位置向量数组，准备存储所有阵元的位置
rm = zeros(3,M);
rm_optimized = zeros(3,M);

% 为每个辐条生成阵元位置
for a = 1:Na
    % 计算该辐条上的阵元角度间隔
    thetaStep = 2*pi/Nm;
    thetaStart = (a-1)*2*pi/Na;
    
    for m = 1:Nm
        theta = thetaStart + (m-1)*thetaStep;
        r = radiusMin + (radiusMax - radiusMin)/(Nm)*(m-1); % 简化假设，线性分布，实际可能需要更复杂的分布逻辑
        rm(:, (a-1)*Nm + m) = [r*cos(theta), r*sin(theta), 0];
    end
end
rm = rm';

for a = 1:Na
    % 计算该辐条上的阵元角度
    thetaStep = 2*pi/Nm;
    thetaStart = (a-1)*2*pi/Na;
    theta = thetaStart + (a-1)*thetaStep;
    
    for m = 1:Nm
        % 对于所有辐条，都使用第一行的相同r值序列
        r = r_optimized(m);  % 从第一行提取r值，循环使用
        rm_optimized(:, (a-1)*Nm + m) = [r*cos(theta), r*sin(theta), 0];
    end
end
rm_optimized = rm_optimized';

%% V_db

load('kx_ky_KZ.mat');

% 计算方向图矩阵V
V = zeros(1, length(kx));
V_optimized = zeros(1, length(kx));
for i = 1:length(kx)

        % 验证 kx 平方和是否小于1
        if kx(i)^2  < 1
            % 计算三维单位向量 \kappa
            kappa = [kx(i), 0, KZ(i, length(kx)/2+1)];
            V(i) = sum(exp(1j*k*rm*kappa.')); % 注意这里的转置，保证维度匹配
            V_optimized(i) = sum(exp(1j*k*rm_optimized*kappa.'));
        else
            % 如果不满足条件，将该位置的值设为0
            V(i) = 0;
            V_optimized(i) = 0;
        end

end

    V_abs = abs(V);
    maxVal = max(max(V_abs(:)));
    V_db = 20*log10(V_abs/maxVal);

    V_optimized_abs = abs(V_optimized);
    maxVal = max(max(V_optimized_abs(:)));
    V_optimized_db = 20*log10(V_optimized_abs/maxVal);

%% Direction_Graph
% 定义颜色
colorOriginal = [47, 127, 193]/255; % 不饱和蓝色
colorOptimized = [255, 127, 111]/255; % 不饱和红色

% 绘制原始数据曲线
plot(kx, V_db, 'color', colorOriginal, 'LineWidth', 1.5);
hold on; % 保持图形以便在同一图上绘制另一条线

% 绘制优化后数据曲线
plot(kx, V_optimized_db, 'color', colorOptimized, 'LineWidth', 1.5);

% 设置图例
legend('Original', 'Optimized', 'Location', 'best');

% 设置横纵坐标标签
xlabel('kx (Normalized Distance)', 'FontSize', 12);
ylabel('Amplitude (dB)', 'FontSize', 12);

% 设置图的标题（可选）
title('Comparison of Original and Optimized Amplitude Responses', 'FontSize', 14);

% 确保横坐标是对齐的，如果kx不是等间距的，可能需要调整 xlim
xlim([min(kx), max(kx)]);

% 调整网格线（可选）
grid on;

% 优化图表布局
legend('boxoff'); % 移除图例边框
set(gca, 'XTickLabel', {}, 'YTickLabel', {}); % 清除刻度标签，如果需要刻度则注释掉这行
set(gcf, 'Color', 'w'); % 设置背景为白色

% 显示图表
hold off; % 结束保持状态