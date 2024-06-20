function fitnessValue = fitnessFunction(RValues, K)
    % 调用analyzeArrayGeometry函数获取B_{3dB}和M_{side}
    [threeDBBandwidth, maxSidelobeVal] = analyzeArrayGeometry(RValues);
    
    % 根据给定的公式计算适应度值
    fitnessValue = (1/threeDBBandwidth) + (K/maxSidelobeVal) + 5;
    % 检查fitnessValue是否为NaN，并处理
    if isnan(fitnessValue)
        fitnessValue = 0; % 如果是NaN，则置为0
    end
    
    % 注意：由于我们要最大化适应度值，而原公式会导致越小的值适应度“越高”，
    % 实际上我们可能需要对公式稍作调整，或是在遗传算法框架内正确处理这种最大化问题。
end


function [threeDBBandwidth, maxSidelobeVal] = analyzeArrayGeometry(RValues)
% Function to analyze microphone array geometry and compute 3dB bandwidth and max sidelobe level

Na = 16; % 辐条数量
Nm = length(RValues);  % 每个辐条的阵元数
M = Na * Nm; % 总阵元数量
frequency = 1000; % 入射波频率
c = 343; % 声速

% 计算波数
k = 2*pi*frequency/c;

% 初始化阵元位置向量数组，准备存储所有阵元的位置
rm = zeros(3,M);

%% 为每个辐条生成阵元位置
for a = 1:Na
    % 计算该辐条上的阵元角度
    thetaStep = 2*pi/Nm;
    thetaStart = (a-1)*2*pi/Na;
    theta = thetaStart + (0:(Nm-1))*thetaStep;
    
    for m = 1:Nm
        % 对于所有辐条，使用给定的r值序列
        r = RValues(m);  % 从输入参数提取r值
        rm(:, (a-1)*Nm + m) = [r*cos(theta(m)), r*sin(theta(m)), 0];
    end
end

rm = rm';

%% 加载kx, ky网格
% 指定文件名
filename = 'kx_ky_KZ.mat';

% 尝试加载已保存的变量
if exist(filename, 'file') == 2 % 检查文件是否存在
    load(filename);
%     fprintf('Loaded existing kx, ky, and KZ from %s.\n', filename);
else
    % 初始化kx, ky网格
    kx = linspace(-1, 1, 360);
    ky = linspace(-1, 1, 360);
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
    
    % 保存结果
    save(filename, 'kx', 'ky', 'KZ');
%     fprintf('Computed and saved kx, ky, and KZ to %s.\n', filename);
end

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

% 3dB带宽计算
row_of_interest = V_db;
[maxVal, maxValPos] = max(row_of_interest);
three_dB_down = maxVal - 3;
lowerBound = find(row_of_interest(1:maxValPos-1) >= three_dB_down, 1, 'first'); 
if isempty(lowerBound)
    lowerBound = 1; % 如果没有找到满足条件的点，设为序列起始
end
upperBound = find(row_of_interest(maxValPos:end) >= three_dB_down, 1, 'last'); 
if isempty(upperBound)
    upperBound = length(row_of_interest); % 如果没有找到满足条件的点，设为序列结束
else
    upperBound = maxValPos + upperBound - 1; % 调整索引以反映原序列的位置
end
if lowerBound > upperBound
    fprintf('No distinct 3dB boundaries found around the peak.\n');
else
    threeDBBandwidth = abs(kx(upperBound) - kx(lowerBound)); % 假设kx是等间距的频率轴
end

% 最大旁瓣级
[pks, ~] = findpeaks(row_of_interest);
[~, idxMaxPeak] = max(pks);
pksWithoutMax = pks;
pksWithoutMax(idxMaxPeak) = NaN;
[maxSidelobeVal, ~] = max(pksWithoutMax);

end


    
    
    
    