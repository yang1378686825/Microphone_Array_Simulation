clc
clear
close all
%%
% [y, Fs] = audioread('_daji_mini_9_fixed.wav'); 
% if size(y, 2) == 2 % 检查是否为双通道
%     y = mean(y, 2); % 转换为单声道
% end

% startSample = round(length(y)- 0.5* Fs); 
% endSample = length(y); 
% y_clipped = y(startSample+1:endSample);
% audiowrite('noise.wav', y_clipped, Fs); 
%%
% [y, Fs] = audioread('_daji_mini_9_fixed.wav'); % 替换为您的输入文件名
% if size(y, 2) == 2 % 检查是否为双通道
%     y = mean(y, 2); % 转换为单声道
% end
% 
% startSample = round(1.5* Fs); 
% endSample = round(2* Fs); 
% y_clipped = y(startSample+1:endSample);


[y, Fs] = audioread('Membo_merged_audio.wav'); 
% 检查并转换为单声道（如果需要）
if size(y, 2) == 2
    y = mean(y, 2);
end

% 定义起始时间和结束时间（秒）
startTime = 0.5;
endTime = 6;

% 定义时间间隔（秒）
timeInterval = 0.5;

% 计算起始和结束的样本点
startSample = round(startTime * Fs);
endSample = round(endTime * Fs);

% 计算总共能分割的段落数
numSegments = ceil((endTime - startTime) / timeInterval);

% 初始化用于存储每个段落的矩阵
y_clipped = zeros(numSegments, timeInterval*Fs);

% 分割音频段落
for i = 1:numSegments
    % 计算当前段落的起始样本点
    segmentStartSample = startSample + (i-1) * round(timeInterval * Fs);
    % 确保最后一个段落不会超出总长度
    segmentEndSample = min(segmentStartSample + round(timeInterval * Fs), endSample);
    
    % 提取当前段落并存储
    y_clipped(i, :) = y(segmentStartSample:segmentEndSample-1);
end

%% 

y_clipped = y_clipped';
y = y_clipped(:, 6);
y = y - mean(y);

% 设定高通滤波器的截止频率为250Hz
cutoffFreq = 250;

% 使用highpass函数进行滤波
[y_filtered, d] = highpass(y, cutoffFreq, Fs, Steepness=0.8); % 'steepness'参数可调整滤波器陡峭度，默认值为0.85，可根据需要调整


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

% 确定窗口长度，这里以信号长度的1/4为例，实际根据需求调整
windowLength = length(y_clean)-1;
if mod(windowLength, 2) == 0 % 确保窗口长度为奇数，以避免相位失真
    windowLength = windowLength + 1;
end
hammingWindow = blackman(windowLength);

% 将窗函数应用于信号
y_window = y_clean(1:windowLength) .* hammingWindow;


% 计算原始信号的功率谱
[pxx_y, f_y] = pwelch(y, [], [], [], Fs);
[pxx_y_filtered, f_y_filtered] = pwelch(y_filtered, [], [], [], Fs);
[pxx_y_clean, f_y_clean] = pwelch(y_clean, [], [], [], Fs);
[pxx_y_window, f_y_window] = pwelch(y_window, [], [], [], Fs);

pxx_y = pxx_y / max(pxx_y);
pxx_y_filtered = pxx_y_filtered / max(pxx_y_filtered);
pxx_y_clean = pxx_y_clean / max(pxx_y_clean);
pxx_y_window = pxx_y_window / max(pxx_y_window);

% 平滑功率谱
pxx_y_smoothed = movmean(pxx_y_filtered, 5, 'omitmissing'); % 移动平均，窗口大小为10

% 绘制功率谱对比图
figure;
% plot(f_y, 10*log10(pxx_y), 'g.', 'LineWidth', 1.5);
hold on;
plot(f_y_filtered, 10*log10(pxx_y_filtered), 'yo', 'LineWidth', 1.5);
plot(f_y_clean, 10*log10(pxx_y_clean), 'b-', 'LineWidth', 1.5);
plot(f_y_clean, 10*log10(pxx_y_smoothed), 'r-', 'LineWidth', 1.5);
% legend('Original Signal Power Spectrum', 'Filtered Signal Power Spectrum', Location='best');
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');
title('Power Spectral Density Comparison');
xlim([0, Fs/4]); % 仅展示一半的频谱，因为是实信号且采样率为Fs
grid on;

%% 
% 使用 findpeaks 函数找到波峰
coordinate = find(f_y_clean == Fs/4);
pxx_y_smoothed_db = 10*log10(pxx_y_smoothed);
pxx_y_smoothed_db = pxx_y_smoothed_db(1: coordinate);
[peaks, locs] = findpeaks(pxx_y_smoothed_db, 'SortStr', 'descend');

% 筛选出波峰值大于等于-10的波峰
validPeaks = peaks(peaks >= -10);

% 对应的横坐标（波峰位置）
validLocs = f_y_window(locs(peaks >= -10));



