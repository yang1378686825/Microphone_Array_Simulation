
function indices = rouletteWheelSelection(fitnessValues, populationSize)

%             fitnessValues, % 各个体的适应度值数组
%             populationSize % 当前种群的大小


    % 轮盘赌选择法
    cumulativeFitness = cumsum(fitnessValues) / sum(fitnessValues);
    randomNumbers = rand(1, populationSize);
    indices = zeros(1, populationSize);
    for i = 1:populationSize
        
        
        try
            % 尝试执行的代码
            indices(i) = find(cumulativeFitness >= randomNumbers(i), 1, 'first');
        catch exception
            % 如果出错，则执行这里的代码
            disp('错误发生了：')
            disp(exception.message) % 打印出错信息
            disp('cumulativeFitness:')
            disp(cumulativeFitness) % 显示cumulativeFitness的内容
            disp(fitnessValues) % 显示cumulativeFitness的内容
            fprintf('randomNumbers(i): %f\n', randomNumbers(i)) % 显示randomNumbers(i)的值
            try
                fprintf('find()的结果：')
                disp(find(cumulativeFitness >= randomNumbers(i), 1, 'first')) % 尝试单独执行find并打印结果
            catch innerException
                disp('find()也出错了：')
                disp(innerException.message) % 如果find也出错，打印其错误信息
            end
            disp('错误前的indices:')
            disp(indices) % 显示错误发生前indices的状态
        end
        
          
        
    end
end

