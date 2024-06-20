function randomNumbers = Random_Numbers(numSets, numbersPerSet, minValue, maxValue,interval)
% generateRandomNumbers 生成指定数量的随机数数组
%
%   randomNumbers = Random_Numbers(numSets, numbersPerSet, minValue, maxValue)
%
%   输入参数:
%       numSets:      生成数据集的数量(遗传算法的种群中个体的数量)
%       numbersPerSet: 每个数据集中随机数的数量
%       minValue:     随机数可取的最小值（包含）
%       maxValue:     随机数可取的最大值（包含），注意这里的值应该能够被适当缩放以满足后续计算需求
%       interval:     随机数间隔
%
%   输出参数:
%       randomNumbers: 一个numSets x numbersPerSet的矩阵，每行代表一组随机数


% 初始化矩阵来存储结果
randomNumbers = zeros(numSets, numbersPerSet);

% 循环生成指定数量的数据集
for i = 1:numSets
    % 生成不重复的索引，这里先生成0到(最大值-最小值)*10的整数，然后除以10并加上最小值
    idx = randperm(((maxValue - minValue)/interval)-1,numbersPerSet)*interval+ minValue;
    % 直接使用idx，因为我们已经调整了范围，不需要额外加最小值
    randomNumbers(i,:) = idx;
end


end


% % 调用函数
% randomNumbersExample = Random_Numbers(50, 8, 0.1, 0.7, 0.01);

