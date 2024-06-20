clc
clear
close all
%%
filename = '_daji_mini_9_fixed.wav';
[y, Fs] = audioread(filename); 
if size(y, 2) == 2 % 检查是否为双通道
    y = mean(y, 2); % 转换为单声道
end

% 获取_noise
startSample = round(length(y)- 0.5* Fs); 
endSample = length(y); 
y_clipped = y(startSample+1:endSample);
audiowrite('noise.wav', y_clipped, Fs); 

% 获取音频段落
startTime = 0.5;
endTime = 9;
timeInterval = 0.5;

startSample = round(startTime * Fs);
endSample = round(endTime * Fs);
numSegments = ceil((endTime - startTime) / timeInterval);

% 分割音频段落
y_clipped = zeros(numSegments, timeInterval*Fs);
for i = 1:numSegments
    % 计算当前段落的起始样本点
    segmentStartSample = startSample + (i-1) * round(timeInterval * Fs);
    % 确保最后一个段落不会超出总长度
    segmentEndSample = min(segmentStartSample + round(timeInterval * Fs), endSample);
    % 提取当前段落并存储
    y_clipped(i, :) = y(segmentStartSample+1:segmentEndSample);
end

%% 

y_clipped = y_clipped';
results = []; % 初始化一个空的行向量来存储每列处理后的validLocs

for col = 1:size(y_clipped, 2)
    y = y_clipped(:, col); % 取当前列的数据
    y = y - mean(y); % 减去均值

    % 使用highpass函数进行滤波
    cutoffFreq = 250;
    [y_filtered, d] = highpass(y, cutoffFreq, Fs, Steepness=0.8); % 'steepness'参数可调整滤波器陡峭度，默认值为0.85，可根据需要调整
    
    % 加载噪声
    noise_part = audioread('noise.wav');
    signal_part = y_filtered; % 全部信号
    
    % 计算噪声功率谱
    NFFT = 2^nextpow2(length(noise_part)); % 使用合适的FFT点数
    P_noise = abs(fft(noise_part,NFFT)).^2 / NFFT;
    P_noise =  P_noise / max(P_noise);
    
    % 计算信号功率谱
    P_signal = abs(fft(signal_part,NFFT)).^2 / NFFT;
    P_signal =  P_signal / max(P_signal);
    
    % 谱减法
    P_clean = P_signal - 2.*P_noise; 
    P_clean(P_clean < 0) = 0; % 确保P_clean中的值非负
    
    % 逆变换回时域
    y_clean = real(ifft(sqrt(P_clean),NFFT));
    
    % 7. 裁剪至合适长度并重设幅度
    y_clean = y_clean(1:length(signal_part));
    y_clean = y_clean / max(abs(y_clean)); % 重设最大幅值
    
    % 周期图法计算信号的功率谱
    [pxx_y_clean, f_y_clean] = pwelch(y_filtered, [], [], [], Fs);
    pxx_y_clean = pxx_y_clean / max(pxx_y_clean);
    
    % 平滑功率谱
    pxx_y_smoothed = movmean(pxx_y_clean, 5, 'omitmissing'); % 移动平均，窗口大小为5

    % 使用 findpeaks 函数找到波峰
    coordinate = find(f_y_clean == Fs/8);
    pxx_y_smoothed_db = 10*log10(pxx_y_smoothed);
    pxx_y_smoothed_db = pxx_y_smoothed_db(1: coordinate);
    [peaks, locs] = findpeaks(pxx_y_smoothed_db, 'SortStr', 'descend');

    % 筛选出波峰值大于等于-10的波峰
    validPeaks = peaks(peaks >= -10);
    
    % 对应的横坐标（波峰位置）
    validLocs = f_y_clean(locs(peaks >= -10));

    % 将当前列处理得到的validLocs添加
    results = [results; validLocs]; 
end

% %% 导出到Excel
% filename = 'band_peaks_results.xlsx';
% % 读取现有的Excel文件数据
% data = readtable(filename);
% 
% % 检查数据表的最后一列是否为空或确定要追加的位置
% % 这里直接追加到末尾
% VariableNames_list = {'Daji'; 'Bebop'; 'Membo'};
% numColsBefore = size(data, 2); % 获取当前列数
% data = [data array2table(results, 'VariableNames', VariableNames_list(2))]; % 追加新列
% 
% % 将更新后的数据写回Excel文件
% writetable(data, filename);

%% 




