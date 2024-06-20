clc
clear
close all

%% 参数定义

M = 64; % 麦克风阵元数量
radius = 0.4; % 阵列半径
frequency = 1000; % 入射波频率
c = 343; % 声速

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

%% 
% 初始化kx, ky网格
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

