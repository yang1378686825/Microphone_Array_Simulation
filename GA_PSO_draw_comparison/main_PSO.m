clc
clear
close all

% PSO参数设置
filename = 'PSO_results.mat';
particleCount = 40;       % 粒子数量
dimensions = 8;           % 搜索空间维度，对应于8个位置的选择
minVal = 0.1;             % 解空间的最小值
maxVal = 0.7;             % 解空间的最大值
maxVelocity = 0.4;        % 速度的最大值，控制探索范围
inertiaWeight = 0.9;      % 惯性权重，控制速度的继承程度
cognitiveCoefficient = 2;   % 认知系数，影响粒子向个人最优解移动的程度
socialCoefficient = 2;    % 社会系数，影响粒子向全局最优解移动的程度
maxEpochs = 50;           % 最大迭代次数
K_factor = 25;             % 适应度函数设置：(1/threeDBBandwidth) + (K/maxSidelobeVal)

% 新增参数
mutationProbability = 0.5; % 每个粒子变异的概率
resetProbability = 0.2; % 每个粒子重置的概率
cognitiveDecay = 0.99; % 认知系数衰减速率
socialGrowth = 1.01; % 社会系数增长速率
inertiaWeightStart = 0.9; % 初始惯性权重
inertiaWeightEnd = 0.4; % 最终惯性权重
inertiaWeightLinearDecayRate = (inertiaWeightStart - inertiaWeightEnd) / maxEpochs; % 线性递减率
noImprovementThreshold = 5; % 若干代无显著改进的阈值
lastImprovementEpoch = 0; % 上一次改进的代数记录

% 初始化粒子群
particles = Random_Numbers(particleCount, dimensions, minVal, maxVal, 0.01);
velocities = randn(particleCount, dimensions) * maxVelocity; % 初始化速度矩阵
velocities = round(velocities / 0.01) * 0.01;      % 调整速度，确保变化步长为0.01的整数倍
personalBestPositions = particles; % 初始化每个粒子的个人最优位置
personalBestFitness = zeros(particleCount, 1); % 初始化每个粒子的个人最优适应度

% 全局最优初始化
globalBestFitness = -Inf;
globalBestPosition = [];
lastGlobalBestFitness = -Inf;

% 初始化保存结构体
results = struct('population', {}, 'best_fitness', {});



%% PSO主循环

for epoch = 1:maxEpochs
    
    % 计算适应度
    fitnessValues = zeros(particleCount, 1);
    for i = 1:particleCount
        RValues = particles(i, :);
        fitnessValues(i) = fitnessFunction(RValues, K_factor); 
    end
    
    % 更新个人最优解
    betterFitness = fitnessValues > personalBestFitness;
    personalBestFitness(betterFitness) = fitnessValues(betterFitness);
    personalBestPositions(betterFitness,:) = particles(betterFitness,:);
    
    % 更新全局最优解
    if max(fitnessValues) > globalBestFitness
        globalBestFitness = max(fitnessValues);
        [~, idx] = max(fitnessValues);
        globalBestPosition = particles(idx, :);
    end
    
    % 输出当前代信息（可选）
    fprintf('Epoch %d, Global Best Fitness: %f\n', epoch-1, globalBestFitness);
    results(end+1).population = particles;  % 保存当前代的population
    results(end).best_fitness = globalBestFitness;  % 保存当前代的最佳fitness
 
    % 判断连续无显著改进的情况
    if (globalBestFitness == lastGlobalBestFitness) && (epoch - lastImprovementEpoch >= noImprovementThreshold)
        % 识别并标记精英粒子
        eliteMask = fitnessValues == max(fitnessValues);
        nonEliteMask = ~eliteMask; % 非精英粒子的掩码

        % % 对非精英粒子应用自适应变异
        % mutationProbabilityAdaptive = mutationProbability * (epoch - lastImprovementEpoch); % 自适应变异概率，随未改进代数增加
        % mutationMask = nonEliteMask & (rand(particleCount, 1) < mutationProbabilityAdaptive); % 变异粒子掩码，只对非精英粒子应用变异
        % 
        % mutationAmount = 0.75 * (maxVal - minVal) * (2 * rand(particleCount, dimensions) - 1); % 计算实际需要变异的粒子的变异量
        % mutationAmount = round(mutationAmount / 0.01) * 0.01;
        % particles(mutationMask,:) = particles(mutationMask,:) + mutationAmount(mutationMask,:); % 应用变异

        % 重置非精英粒子的逻辑可以保留在这里，但根据您的要求，我们主要关注变异
        % 如果确实需要重置而非变异，请在这里添加相应的代码逻辑
        newPositionsNonElite = Random_Numbers(particleCount - sum(eliteMask), dimensions, minVal, maxVal, 0.01);
        particles(nonEliteMask, :) = newPositionsNonElite;

        
        particles(particles < minVal) = minVal; % 确保位置在范围内
        particles(particles > maxVal) = maxVal;
        
        % 检查并处理每一行的重复值
        for i = 1:particleCount
            if any(diff(sort(particles(i,:)))==0) % 检查是否有重复值
                particles(i,:) = handleDuplicates(particles(i,:), minVal, maxVal); % 使用处理函数
            end
        end
        
    else
        if globalBestFitness > lastGlobalBestFitness
            lastImprovementEpoch = epoch-1; % 更新上一次改进的代数
        end
        lastGlobalBestFitness = globalBestFitness; % 更新历史最优适应度
    end
    
%     % 重置粒子位置的条件判断
%     if (globalBestFitness == lastGlobalBestFitness) && (epoch - lastImprovementEpoch >= noImprovementThreshold)
%     % 若连续若干代无显著改进，重置非精英粒子
%         nonEliteMask = fitnessValues ~= max(fitnessValues);
%         particles(nonEliteMask,:) = Random_Numbers(sum(nonEliteMask), dimensions, minVal, maxVal, 0.01);
%         lastImprovementEpoch = epoch-1; % 更新上一次改进的代数
%     else
%         if globalBestFitness > lastGlobalBestFitness
%             lastImprovementEpoch = epoch-1; % 更新上一次改进的代数
%         end
%         lastGlobalBestFitness = globalBestFitness; % 更新历史最优适应度
%     end
    
    
    % 动态调整惯性权重
    inertiaWeight = max(inertiaWeightStart - inertiaWeightLinearDecayRate * epoch, inertiaWeightEnd);
    
    % 更新速度和位置之前动态调整认知与社会系数
    cognitiveCoefficient = cognitiveCoefficient * cognitiveDecay;
    socialCoefficient = socialCoefficient * socialGrowth;
    
    % 更新速度和位置
    r1 = rand(particleCount, dimensions);
    r2 = rand(particleCount, dimensions);
    velocities = inertiaWeight * velocities + ...
                  cognitiveCoefficient * r1 .* (personalBestPositions - particles) + ...
                  socialCoefficient * r2 .* (globalBestPosition - particles);
    velocities(velocities > maxVelocity) = maxVelocity; % 限制速度
    velocities(velocities < -maxVelocity) = -maxVelocity;
    
    % 调整速度，确保变化步长为0.01的整数倍
    velocities = round(velocities / 0.01) * 0.01;
    
    % 位置更新
    particles = particles + velocities;
    
    
    % 重置部分粒子到随机位置
    resetMask = rand(particleCount, 1) < resetProbability;
    particles(resetMask,:) = Random_Numbers(sum(resetMask), dimensions, minVal, maxVal, 0.01);

    particles(particles < minVal) = minVal; % 确保位置在范围内
    particles(particles > maxVal) = maxVal;
    
    % 检查并处理每一行的重复值
    for i = 1:particleCount
        if any(diff(sort(particles(i,:)))==0) % 检查是否有重复值
            particles(i,:) = handleDuplicates(particles(i,:), minVal, maxVal); % 使用处理函数
        end
    end
     
    
end



%% 输出最终最佳解
% 计算适应度
fitnessValues = zeros(particleCount, 1);
for i = 1:particleCount
    RValues = particles(i, :);
    fitnessValues(i) = fitnessFunction(RValues, K_factor); 
end

% 更新个人最优解
betterFitness = fitnessValues > personalBestFitness;
personalBestFitness(betterFitness) = fitnessValues(betterFitness);
personalBestPositions(betterFitness,:) = particles(betterFitness,:);

% 更新全局最优解
if max(fitnessValues) > globalBestFitness
    globalBestFitness = max(fitnessValues);
    [~, idx] = max(fitnessValues);
    globalBestPosition = particles(idx, :);
end
fprintf('Best solution found: \n');
disp(globalBestPosition);
fprintf('With fitness value: %f\n', globalBestFitness);



%% 在循环结束后，保存results到.mat文件
results(end+1).population = particles;  % 保存当前代的population
results(end).best_fitness = globalBestFitness;  % 保存当前代的最佳fitness
save(filename, 'results');
plotDirectionGraph(filename);


