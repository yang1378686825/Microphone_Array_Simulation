clc
clear
close all
%% 参数定义

% 参数定义
M = 16; % 麦克风阵元数量
radius = 0.1; % 阵列半径
frequency = 1000; % 入射波频率
c = 343; % 声速
resolution_azimuth = 360; % 方位角分辨率
resolution_elevation = 360; % 俯仰角分辨率

% 计算波数
k = 2*pi*frequency/c;

% 初始化阵元位置向量（仅x-y平面，但为了计算方便，添加z=0）
theta = linspace(0, 2*pi, M);
rm = radius*[cos(theta)', sin(theta)', zeros(1,M)']; % 阵元相对于中心的位置向量（三维）

% 绘制阵列形状
figure;
hold on;
plot(rm(:,1), rm(:,2), 'o-', 'MarkerSize', 10, 'LineWidth', 2);
plot([rm(1,1), rm(end,1)], [rm(1,2), rm(end,2)], 'k-', 'LineWidth', 2); % 连接首尾形成圆环
axis equal; % 等比例缩放
xlim([-radius*1.5 radius*1.5]);
ylim([-radius*1.5 radius*1.5]);
title('Circular Microphone Array Layout');
xlabel('X-axis');
ylabel('Y-axis');

%%  计算方向图

% 方位角和俯仰角范围
azimuth = linspace(0, 2*pi, resolution_azimuth); % 转换为弧度
elevation = linspace(0, pi, resolution_elevation); % 俯仰角范围，转换为弧度

% 初始化方向图矩阵
V = zeros(length(azimuth), length(elevation));

for az_idx = 1:length(azimuth)
    for el_idx = 1:length(elevation)
        % 计算三维单位向量 \kappa
        kappa = [cos(azimuth(az_idx))*sin(elevation(el_idx)), sin(azimuth(az_idx))*sin(elevation(el_idx)), cos(elevation(el_idx))];
        V(az_idx, el_idx) = sum(exp(1j*k*rm*kappa.')); % 注意这里的转置，保证维度匹配
    end
end


figure;
imagesc( elevation./pi, azimuth./pi, abs(V));
xlabel('Elevation (pi)');
ylabel('Azimuth (pi)');
title('Beam Pattern Magnitude of Circular Microphone Array');
colorbar;
set(gca, 'YDir', 'normal'); % 确保俯仰角从下到上递增
caxis([0, max(abs(V(:)))]); % 色标范围


