clc
clear
close all
%%
[y, Fs] = audioread('Bebop_merged_audio.wav'); 
y = y - mean(y);

% 设定高通滤波器的截止频率为250Hz
cutoffFreq = 300;

% 使用highpass函数进行滤波
[y_filtered, d] = highpass(y, cutoffFreq, Fs, Steepness=0.9); % 'steepness'参数可调整滤波器陡峭度，默认值为0.85，可根据需要调整

% 计算信号的功率谱
[pxx_y, f_y] = pwelch(y, [], [], [], Fs);
[pxx_y_filtered, f_y_filtered] = pwelch(y_filtered, [], [], [], Fs);

pxx_y = pxx_y / max(pxx_y);
pxx_y_filtered = pxx_y_filtered / max(pxx_y_filtered);

% 平滑功率谱
pxx_y_filtered = movmean(pxx_y_filtered, 10, 'omitmissing'); % 移动平均，窗口大小为10

% 绘制功率谱对比图
figure;
plot(f_y, 10*log10(pxx_y), 'g.', 'LineWidth', 1.5);
hold on;
plot(f_y_filtered, 10*log10(pxx_y_filtered), 'b-', 'LineWidth', 1.5);
legend('Original Signal Power Spectrum', 'Filtered Signal Power Spectrum', Location='best');
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');
title('Power Spectral Density Comparison');
xlim([0, Fs/8]); % 仅展示一半的频谱，因为是实信号且采样率为Fs
grid on;




