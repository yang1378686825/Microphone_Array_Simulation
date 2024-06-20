clc
clear
close all

% 参数设置
populationSize = 40;      % 种群大小
chromosomeLength = 8;     % 染色体长度
minVal = 0.1;                % 基因取值最小范围
maxVal = 0.7;              % 基因取值最大范围
mutationRate = 0.5;       % 变异率
epoch = 100;               % 进化代数
K_factor = 50;            % 适应度函数设置：(1/threeDBBandwidth) + (K/maxSidelobeVal)
filename = 'GA_results.mat';

% 调用遗传算法
geneticAlgorithm(populationSize, chromosomeLength, minVal, maxVal, mutationRate, epoch, K_factor, filename);
plotDirectionGraph(filename);


function geneticAlgorithm(populationSize, chromosomeLength, minVal, maxVal, mutationRate, epoch, K_factor, filename)

% 初始化种群
population = Random_Numbers(populationSize, chromosomeLength, minVal, maxVal, 0.01);
% % 初始化上一代最佳适应度
% prevBestFitness = -Inf;
% 初始化保存结构体
results = struct('population', {}, 'best_fitness', {});

for generation = 1:epoch
    
    % 计算适应度
    fitnessValues = zeros(populationSize, 1);
    for i = 1:populationSize
        RValues = population(i, :);
        fitnessValues(i) = fitnessFunction(RValues, K_factor); 
    end
    
    % 输出并保存当前代信息（可选）
    fprintf('Generation %d, Best Fitness: %f\n', generation-1, max(fitnessValues));
    results(end+1).population = population;  % 保存当前代的population
    results(end).best_fitness = max(fitnessValues);  % 保存当前代的最佳fitness
    
    % 精英保留策略：保留当前代的前10%最优个体
    eliteCount = round(populationSize * 0.1);
    [~, sortedIndices] = sort(fitnessValues, 'descend');
    elites = population(sortedIndices(1:eliteCount,:), :);
    
%     % 检查停止条件：适应度差值是否足够小
%     currentBestFitness = fitnessValues(sortedIndices(1));
%     if abs(currentBestFitness - prevBestFitness) < 0.01
%         fprintf('Termination condition met: Improvement in best fitness is less than 0.01.\n');
%         break; % 跳出循环
%     end
%     % 更新上一代最佳适应度
%     prevBestFitness = currentBestFitness;
    

    % 选择操作：轮盘赌选择法
    parentIndices = rouletteWheelSelection(fitnessValues, populationSize);
    
    % 交叉操作
    newPopulation = zeros(populationSize - eliteCount, chromosomeLength);
    for i = 1:2:(populationSize - eliteCount - 1)
        crossoverPoint = randi([1, chromosomeLength - 1]);
        parent1 = population(parentIndices(i), :);
        parent2 = population(parentIndices(i+1), :);
        offspring1 = [parent1(1:crossoverPoint), parent2(crossoverPoint+1:end)];
        offspring2 = [parent2(1:crossoverPoint), parent1(crossoverPoint+1:end)];
        offspring1 = handleDuplicates(offspring1, minVal, maxVal);
        offspring2 = handleDuplicates(offspring2, minVal, maxVal);
        newPopulation(i, :) = offspring1;
        newPopulation(i+1, :) = offspring2;
    end

    % 变异操作
    for i = 1:populationSize - eliteCount
        if rand < mutationRate
            mutationPoint = randi([1, chromosomeLength]);
            newPopulation(i, mutationPoint) = round((rand() * (maxVal - minVal) + minVal) * 100) / 100;
            newPopulation(i,:) = handleDuplicates(newPopulation(i,:), minVal, maxVal);
        end
    end
      
    % 合并精英与新生成的个体
    population = [elites; newPopulation];
   
    
end
    
% 首先，我们需要重新计算合并精英和新生成个体后的种群的整体适应度
fitnessValuesFinal = zeros(populationSize, 1);
for i = 1:populationSize
    RValuesFinal = population(i, :);
    fitnessValuesFinal(i) = fitnessFunction(RValuesFinal, K_factor); 
end

% 然后，找到具有最高适应度的个体的索引
[~, bestFinalIndex] = max(fitnessValuesFinal);

fprintf('Best solution found: \n');
disp(population(bestFinalIndex, :)); % 输出最佳个体
fprintf('With fitness value: %f\n', fitnessValuesFinal(bestFinalIndex)); % 输出最佳适应度

% 在循环结束后，保存results到.mat文件
results(end+1).population = population;  % 保存当前代的population
results(end).best_fitness = fitnessValuesFinal(bestFinalIndex);  % 保存当前代的最佳fitness
save(filename, 'results');


end



