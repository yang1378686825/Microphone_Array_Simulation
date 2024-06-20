clear all;
close all;
clc;

frequency = 1000; % 入射波频率
c = 343; % 声速
% 计算波数
k = 2*pi*frequency/c;

% Dougherty多臂阵列设计参数
r_max = 0.7;         % 最大螺旋半径 (m)
r_min = 0.1;          % 最小螺旋半径 (m)
v = 5*pi/16;          % 螺旋切线与半径的夹角 (rad)
Nm = 8;               % 每条臂上的阵元数
Na = 16;              % 总臂数

% 中间变量计算函数
ln_func = @(n) ((n-1)/(Nm-1)) * (r_max * sqrt(1+cot(v)^2)) / cot(v);

% 初始化阵列位置矩阵
theta_1_n = zeros(1, Nm);   % 第一条臂上的阵元角度
r_1_n = zeros(1, Nm);       % 第一条臂上的阵元半径
theta_m_n = zeros(Na, Nm);  % 多臂阵列中所有阵元的角度
r_m_n = zeros(Na, Nm);      % 多臂阵列中所有阵元的半径

% 计算第一条臂上每个阵元的位置
for n = 1:Nm
    l_n = ln_func(n);
    theta_1_n(n) = 1/cot(v) * log(1 + cot(v)*l_n/(r_min*sqrt(1+cot(v)^2)));
    r_1_n(n) = r_min * exp(cot(v)*theta_1_n(n));
end

% 扩展到所有臂
for m = 1:Na
    for n = 1:Nm
        theta_m_n(m, n) = theta_1_n(n) + (m-1)/Na * 2*pi;
        r_m_n(m, n) = r_1_n(n);  % 在多臂阵列中，所有臂的r相同，只改变角度
    end
end

% 计算x,y坐标
[x_1d, y_1d] = pol2cart(theta_m_n, r_m_n);
x_1d = x_1d(:);
y_1d = y_1d(:);
z_1d = zeros(length(x_1d),1);
rm = [x_1d, y_1d, z_1d];

figure;
hold on;
for m = 1:Na
    scatter(x_1d, y_1d, 50,[0, 0.4470, 0.7410], 'filled');
end
axis equal;
xlim([-r_max*1.1 r_max*1.1]);
ylim([-r_max*1.1 r_max*1.1]);

title('Dougherty Multi-Arm Spiral Array Layout');
xlabel('X-axis (m)');
ylabel('Y-axis (m)');
grid off;

hold off;


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


