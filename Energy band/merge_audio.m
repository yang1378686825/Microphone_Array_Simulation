clc
clear
close all

% 定义音频文件的基本名称和预期的总数
baseFileName = 'Bebop_';
totalFiles = 6;
outputFileName = 'Bebop_merged_audio.wav';

% 初始化一个空的音频数据容器
allAudioData = [];

% 循环遍历每个音频文件
for i = 1:totalFiles
    % 构建当前文件的完整路径
    fileName = [baseFileName num2str(i) '.wav'];
    
    % 读取音频文件
    [audioData, sampleRate] = audioread(fileName);
    
    % 将当前音频数据添加到总数据中
    allAudioData = [allAudioData; audioData];
end

% 写入合并后的音频文件
audiowrite(outputFileName, allAudioData, sampleRate);

fprintf('所有音频文件已成功合并至：%s\n', outputFileName);