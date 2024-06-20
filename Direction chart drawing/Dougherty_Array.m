clear all;
close all;
clc;

frequency = 1000; % 入射波频率
c = 343; % 声速
% 计算波数
k = 2*pi*frequency/c;


% 螺旋参数
r0 = 0.25;         % 最小半径 (m)
rMax = 0.75;       % 最大半径 (m)
dTheta = 2*pi/48;  % 每个阵元间的螺旋角度增量
d_r = 0.005;

% 初始化
theta = 0;          % 初始角度
rm = [];            % 用于存储阵元位置的矩阵

% 生成阵列位置直到达到最大半径
while r0 <= rMax
    x = r0 * cos(theta);   % 计算x坐标
    y = r0 * sin(theta);   % 计算y坐标
    rm = [rm; x, y, 0];       % 添加到阵列
    r0 = r0 + d_r;      % 更新半径以生成螺旋效果
    theta = theta + dTheta/(2.25*r0); % 增加角度
end

figure;
hold on;

% 绘制螺旋阵列
scatter(rm(:,1), rm(:,2),  'SizeData', 75, 'LineWidth', 2); % 点的大小设为50，填充样式

% 设置X轴和Y轴的极限
xlim([-rMax*1.1 rMax*1.1]);
ylim([-rMax*1.1 rMax*1.1]);

title('Dougherty Array Spiral Microphone Array Layout');
xlabel('X-axis (m)');
ylabel('Y-axis (m)');

grid off;
hold off; % 结束hold状态


%% 初始化kx, ky网格
kx = linspace(-1, 1, 720);
ky = linspace(-1, 1, 720);
[KX,KY] = meshgrid(kx, ky);

% 创建逻辑掩码，保留kx^2 + ky^2 <= 1的点
valid_points_mask = (KX.^2 + KY.^2) <= 1;

% 应用掩码到KX和KY，移除无效点
KX = KX .* valid_points_mask;
KY = KY .* valid_points_mask;

% 计算对应的kz，现在所有点都满足kx^2 + ky^2 <= 1
KZ = sqrt(1 - KX.^2 - KY.^2); % 此处不会出现复数，因为已预先筛选
% 将KZ中接近1的值（考虑到浮点运算误差，这里设定一个很小的容差eps）置为0
eps = 1e-8; % 设定一个很小的正数作为容差
KZ(abs(KZ - 1) < eps) = 0;

%% 计算方向图
V = zeros(length(kx), length(ky));
for i = 1:length(kx)
    for j = 1:length(ky)
        % 验证 kx 和 ky 的平方和是否小于1
        if (kx(i)^2 + ky(j)^2) < 1
            % 计算三维单位向量 \kappa
            kappa = [kx(i), ky(j), KZ(i, j)];
            V(i, j) = sum(exp(1j*k*rm*kappa.')); % 注意这里的转置，保证维度匹配
        else
            % 如果不满足条件，将该位置的值设为0
            V(i, j) = 0;
        end
    end
end

V_abs = abs(V);
maxVal = max(max(V_abs(:)));
V_db = 20*log10(V_abs/maxVal);

%% 绘制方向图
figure;
imagesc(kx, ky, abs(V));
xlabel('kx');
ylabel('ky');
title('Direction Diagram--Circular Array ');
colormap('jet');
colorbar;
caxis('auto'); % 色标范围
% caxis([0, max(max(V_abs))]);
caxis('auto');


figure;
imagesc(kx, ky, V_db);
xlabel('kx');
ylabel('ky');
title('Direction Diagram--Circular Array ');
% 使用更明显的colormap，例如'jet'，并自动调整色标范围
colormap('jet'); % 或选择其他你喜欢的colormap，如'parula'
colorbar;
% 自动调整色标范围以覆盖数据的全部动态范围
caxis([-30, max(max(V_db))]);
% caxis('auto');
