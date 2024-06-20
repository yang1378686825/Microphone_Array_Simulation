
function processedOffspring = handleDuplicates(offspring, minValue, maxValue)
    % 获取去重后的数组C和每个唯一值的首次出现位置ia
    [~, ia] = unique(offspring, 'stable');

    % 初始化一个新的数组来存储处理过的唯一值
    processedOffspring = zeros(size(offspring));

    % 遍历ia，根据ia的索引将原数组中的唯一值复制到新数组中
    for k = 1:numel(ia)
        processedOffspring(ia(k)) = offspring(ia(k));
    end

    % 如果有重复值，它们在processedOffspring中仍为0，现在对这些位置生成新值
    missingIndices = find(processedOffspring == 0);

    for idx = missingIndices
        while true
            newValue = round((rand() * (maxValue - minValue) + minValue) * 100) / 100;
            if ~ismember(newValue, processedOffspring)
                processedOffspring(idx) = newValue;
                break;
            end
        end
    end
end