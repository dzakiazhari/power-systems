close all; clc; clf; clear all; 

%% Input Variabel
% Read Audio x(t)
Fs      = 2;            % Sampling Frequency
[x, Fs] = audioread('spesifikasic.wav');

% Parameter Variabel
L   = length(x);       % Signal Length
T   = 1/Fs;            % Sampling Period
t   = (0:L-1)*T;       % Time Vector
t_plt= (0:2*(L-1))*T;
N   = length(t);
n   = (0:1:N-1);

%% TUGAS 1
% Spesifikasi C: Susutan 72%
h1 = 0.28*(n==0);
y1 = conv(x,h1);

% PLOT 1 : X(t) Isyarat Masukan Audio
figure(1)
subplot(2,1,1);
plot(t,x);
title('Transmitted Signal Time Domain');
xlabel('Time')
ylabel('Magnitude')
[axis, X]=plot_freq(x,L,Fs);
title('Transmitted Signal Freq. Domain');


% PLOT 2 : H(t) Tanggapan Impuls Susutan
figure(2)
subplot(2,1,1);
stem(h1);
title('Impulse Response Time Domain');
xlabel('Time')
ylabel('Magnitude')
[axis, H1]=plot_freq(h1,L,Fs);
title('Impulse Response Freq. Domain');

% PLOT 3 : Y(t) Output Hasil Susutan
figure(3)
subplot(2,1,1);
plot(t_plt,y1);
title('Received Signal Time Domain');
xlabel('Time')
ylabel('Magnitude')
[axis, Y1]=plot_freq(y1,2*L,Fs);
title('Received Signal Freq. Domain');

filename = 'ytugas1.wav';
audiowrite(filename,y1,Fs);

%% TUGAS 2
[H,channel_axis,H_freq]=wireless_channel(L,Fs);

y2 = conv(x,H);

% PLOT 4 : X(t) Isyarat Masukan Audio
figure(4)
subplot(2,1,1);
plot(t,x);
title('Transmitted Signal Time Domain');
xlabel('Time ')
ylabel('Magnitude')
[axis, X]=plot_freq(x,L,Fs);
title('Transmitted Signal Freq. Domain');


% PLOT 5 : H(t) Tanggapan Impuls Susutan
figure(5);
title('Impulse Response');
subplot(2,1,1);
stem(H);
title('Impulse Response Time Domain');
subplot(2,1,2);
xlabel('Time')
ylabel('Magnitude')
plot(channel_axis, H_freq);
title('Impulse Response Freq. Domain');


% PLOT 6 : Y(t) Output Hasil Susutan
figure(6)
subplot(2,1,1);
plot(t_plt,y2);
title('Received Signal Time Domain');
xlabel('Time')
ylabel('Magnitude')
[axis, Y2]=plot_freq(y2,2*L,Fs);
title('Received Signal Freq. Domain');


