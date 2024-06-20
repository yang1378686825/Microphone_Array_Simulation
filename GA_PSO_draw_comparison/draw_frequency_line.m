clc
clear
close all

% 主函数调用
frequency_line = [500, 100, 1500];

    % 加载必要的数据
    filename = 'GA_results.mat';
    load(filename);
    load('kx_ky_KZ.mat');
    plotDirectionGraph(filename)

    % 初始化存储变量
    threeDBBandwidths = [];
    maxSidelobeVals = [];
    threeDBBandwidthsOptimized = [];
    maxSidelobeValsOptimized = [];

    for frequency = frequency_line(1):frequency_line(2):frequency_line(3)
        disp(['Processing frequency: ', num2str(frequency)]);
        % 更新与频率相关的参数
        c = 343; % 声速
        k = 2*pi*frequency/c;
        [rm, rm_optimized] = calculateArrayPositions(results);
        [V, V_optimized, kx] = calculateBeamPatterns(rm, rm_optimized, kx, ky, KZ, k);
        [threeDBBandwidth, maxSidelobeVal, threeDBBandwidth_optimized, maxSidelobeVal_optimized] = evaluatePerformance(V, V_optimized, kx);
            
        % 存储结果
        threeDBBandwidths = [threeDBBandwidths, threeDBBandwidth];
        maxSidelobeVals = [maxSidelobeVals, maxSidelobeVal];
        threeDBBandwidthsOptimized = [threeDBBandwidthsOptimized, threeDBBandwidth_optimized];
        maxSidelobeValsOptimized = [maxSidelobeValsOptimized, maxSidelobeVal_optimized];
        
    end

    % 定义颜色  
    colorOriginal = [47, 127, 193]/255; % 不饱和蓝色  
    colorOptimized = [255, 127, 111]/255; % 不饱和红色  
      
    % 绘制比较图 
    x = frequency_line(1): frequency_line(2): frequency_line(3);
    figure;  
      
    % 第一个子图：3dB带宽比较  
    subplot(2,1,1);  
    plot(x, threeDBBandwidths, 'LineStyle', '-', 'Color', colorOriginal, 'Marker', 'o', 'LineWidth', 2, 'MarkerSize', 8);  
    hold on;  
    plot(x, threeDBBandwidthsOptimized, 'LineStyle', '-', 'Color', colorOptimized, 'Marker', 's', 'LineWidth', 2, 'MarkerSize', 8);  
    title('3dB Bandwidth Comparison');  
    legend('Original', 'Optimized','Location', 'best');  
    xlabel('Frequency (Hz)');  
    ylabel('3dB Bandwidth');  
      
    % 第二个子图：最大旁瓣值比较  
    subplot(2,1,2);  
    plot(x, maxSidelobeVals, 'LineStyle', '-', 'Color', colorOriginal, 'Marker', 'o', 'LineWidth', 2, 'MarkerSize', 8);  
    hold on;  
    plot(x, maxSidelobeValsOptimized, 'LineStyle', '-', 'Color', colorOptimized, 'Marker', 's', 'LineWidth', 2, 'MarkerSize', 8);  
    title('Max Sidelobe Value Comparison');  
    legend('Original', 'Optimized','Location', 'best');  
    xlabel('Frequency (Hz)');  
    ylabel('Max Sidelobe Value (dB)');  
    hold off;



function [rm, rm_optimized] = calculateArrayPositions(results)
% 根据频率计算原始和优化的阵列位置

    % 优化后阵列数据
    a_ = results.population;
    r_optimized = a_(end,:);

    % 均匀阵列数据
    Na = 16; % 辐条数量
    Nm = 8;  % 每个辐条的阵元数
    M = Na * Nm; % 总阵元数量
    radiusMin = 0.1; % 最小半径
    radiusMax = 0.7; % 最大半径

    % 初始化阵元位置向量数组，准备存储所有阵元的位置
    rm = zeros(M,3);
    rm_optimized = zeros(M,3);

    for a = 1:Na
        % 计算该辐条上的阵元角度间隔
        thetaStep = 2*pi/Nm;
        thetaStart = (a-1)*2*pi/Na;
        theta = thetaStart + (a-1)*thetaStep;

        for m = 1:Nm  
            r = radiusMin + (radiusMax - radiusMin)/(Nm)*(m-1); % r均匀分布
            rm((a-1)*Nm + m, :) = [r*cos(theta), r*sin(theta), 0];
        end
    end

    for a = 1:Na
        % 计算该辐条上的阵元角度
        thetaStep = 2*pi/Nm;
        thetaStart = (a-1)*2*pi/Na;
        theta = thetaStart + (a-1)*thetaStep;

        for m = 1:Nm
            r = r_optimized(m);  % 提取第m个阵元的r值
            rm_optimized((a-1)*Nm + m, :) = [r*cos(theta), r*sin(theta), 0];
        end
    end

end


function [V_db, V_optimized_db, kx] = calculateBeamPatterns(rm, rm_optimized, kx, ky, KZ, k)

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
end


function [threeDBBandwidth, maxSidelobeVal, threeDBBandwidth_optimized, maxSidelobeVal_optimized] = evaluatePerformance(V_db, V_optimized_db, kx)

    % 提取第181行
    row_of_interest = V_db;
    row_of_interest_optimized = V_optimized_db;
    
    % 调用calculateMetrics函数
    [threeDBBandwidth, maxSidelobeVal] = calculateMetrics(row_of_interest);
    [threeDBBandwidth_optimized, maxSidelobeVal_optimized] = calculateMetrics(row_of_interest_optimized);

end


function [threeDBBandwidth, maxSidelobeVal] = calculateMetrics(row_of_interest)
% 找到峰值位置（主瓣）
    [maxVal, maxValPos] = max(row_of_interest);

    % 定义3dB下降点的dB值（相对于峰值）
    three_dB_down = maxVal - 3;

    % 从峰值点开始向前搜索
    lowerBound = find(row_of_interest(1:maxValPos-1) >= three_dB_down, 1, 'first'); 
    if isempty(lowerBound)
        lowerBound = 1; % 如果没有找到满足条件的点，设为序列起始
    end

    % 从峰值点开始向后搜索
    upperBound = find(row_of_interest(maxValPos:end) >= three_dB_down, 1, 'last'); 
    if isempty(upperBound)
        upperBound = length(row_of_interest); % 如果没有找到满足条件的点，设为序列结束
    else
        upperBound = maxValPos + upperBound - 1; % 调整索引以反映原序列的位置
    end

    % 计算3dB带宽
    if lowerBound > upperBound
        fprintf('No distinct 3dB boundaries found around the peak.\n');
    else
        % 假设kx是等间距的频率轴，此处需要确保row_of_interest与kx的对应关系
        % 由于kx未在本函数中未给出，需确保调用时外部已知或传递对应关系
        threeDBBandwidth = abs(upperBound - lowerBound);
    end

    % 使用findpeaks函数找到峰值
    [pks, ~] = findpeaks(row_of_interest);

    % 找到最大峰值的索引和值
    [~, idxMaxPeak] = max(pks);

    % 移除最大峰值，找到次大峰值（第一旁瓣）
    pksWithoutMax = pks;
    pksWithoutMax(idxMaxPeak) = NaN; % 或者直接移除最大值
    [~, locMaxSidelobe] = find(~isnan(pksWithoutMax), 1, 'first'); % 找到第一个非NaN值的位置
    maxSidelobeVal = pks(locMaxSidelobe); % 第一旁瓣的幅值
end



