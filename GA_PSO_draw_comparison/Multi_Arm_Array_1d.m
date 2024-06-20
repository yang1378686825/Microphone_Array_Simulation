clc
clear
close all

%% 参数定义

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

% 为每个辐条生成阵元位置
for a = 1:Na
    % 计算该辐条上的阵元角度间隔
    thetaStep = 2*pi/Nm;
    thetaStart = (a-1)*2*pi/Na;
    theta = thetaStart + (a-1)*thetaStep;
    
    for m = 1:Nm
        r = radiusMin + (radiusMax - radiusMin)/(Nm)*(m-1); % 简化假设，线性分布，实际可能需要更复杂的分布逻辑
        rm(:, (a-1)*Nm + m) = [r*cos(theta), r*sin(theta), 0];
    end
end

rm = rm';

figure;
hold on;

% 绘制每个阵元的位置
scatter(rm(:,1), rm(:,2), 'SizeData', 75, 'LineWidth', 2); % 散点图展示每个阵元


xlim([-radiusMax*1.1 radiusMax*1.1]);
ylim([-radiusMax*1.1 radiusMax*1.1]);

title('Multi-Arm Spoke Microphone Array Layout (2D)');
xlabel('X-axis');
ylabel('Y-axis');

grid off; % 关闭网格以减少视觉杂乱


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
V = zeros(1, length(kx));

for i = 1:length(kx)

        % 验证 kx 平方和是否小于1
        if kx(i)^2  < 1
            % 计算三维单位向量 \kappa
            kappa = [kx(i), 0, KZ(i, length(kx)/2+1)];
            V(i) = sum(exp(1j*k*rm*kappa.')); % 注意这里的转置，保证维度匹配
        else
            % 如果不满足条件，将该位置的值设为0
            V(i) = 0;
        end

end

V_abs = abs(V);
maxVal = max(max(V_abs(:)));
V_db = 20*log10(V_abs/maxVal);




