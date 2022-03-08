clear; clc; clf; close all;

%% Input variabel
Fs  = 1000;            % Sampling Frequency
L   = 1000;            % Signal Length
T   = 1/Fs;            % Sampling Period
t   = (0:L-1)*T;       % Time Vector
N   = length(t);
n   = (0:1:N-1);

% parameter sinyal
A_m = 20;
A_c = 1;
f_m = 10;
f_c = 450;
m_t = A_m*cos(2*pi*f_m*t); %Sinyal Isyarat
c_t = A_c*cos(2*pi*f_c*t); %SInyal Carier
u_t = m_t .* c_t;          %Sinyal Modulasi AM
f1 = f_c + f_m;
f2 = f_c - f_m;
f3 = 60;
theta1 = (2*pi*f1)/Fs;
theta2 = (2*pi*f2)/Fs;
theta3 = (2*pi*f3)/Fs;

%% GAMBAR 1 : PLOT m(t)
figure(1)
subplot(2,1,1)
plot(t,m_t)
xlabel('Time (seconds)')
ylabel('Amplitude (Volt or Ampere)')

m_fft=fft(m_t); 
m_fft_magnitude=abs(m_fft)./L; 
m_fft_single=m_fft_magnitude(1:(0.5*L)+1); 
m_fft_single(2:(0.5*L)+1)=2*m_fft_single(2:(0.5*L)+1); 
axis_single=(0:L/2)*(Fs/L);
subplot(2,1,2)
plot(axis_single,m_fft_single);
xlabel('Frequency (Hz)')
ylabel('Magnitude (Volt or Ampere)')

%% GAMBAR 2 : PLOT c(t)
figure(2)
subplot(2,1,1)
plot(t,c_t)
xlabel('Time (seconds)')
ylabel('Amplitude (Volt or Ampere)')

c_fft=fft(c_t); 
c_fft_magnitude=abs(c_fft)./L; 
c_fft_single=c_fft_magnitude(1:(0.5*L)+1); 
c_fft_single(2:(0.5*L)+1)=2*c_fft_single(2:(0.5*L)+1); 
axis_single=(0:L/2)*(Fs/L); 

subplot(2,1,2)
plot(axis_single,c_fft_single);
xlabel('Frequency (Hz)')
ylabel('Magnitude (Volt or Ampere)')

%% GAMBAR 3 : PLOT u(t)
figure(3)
subplot(2,1,1)
plot(t,u_t)
xlabel('Time (seconds)')
ylabel('Amplitude (Volt or Ampere)')

u_fft=fft(u_t); 
u_fft_magnitude=abs(u_fft)./L; 
u_fft_single=u_fft_magnitude(1:(0.5*L)+1); 
u_fft_single(2:(0.5*L)+1)=2*u_fft_single(2:(0.5*L)+1); 
axis_single=(0:L/2)*(Fs/L); 

subplot(2,1,2)
plot(axis_single,u_fft_single)
xlabel('Frequency (Hz)')
ylabel('Magnitude (Volt or Ampere)')

%% USB Filter
A   = (theta2/pi).*sinc(theta2.*(n-0.5*N)/pi);
B   = (theta1/pi).*sinc(theta1.*(n-0.5*N)/pi);
h_n = A - B;
USB = conv(u_t,h_n,'same');

%% GAMBAR 4 : PLOT frequency response filter h(n)
figure(4);
[h,w]=freqz(h_n,1);
plot(abs(h));

%% GAMBAR 5 : PLOT Hasil filter 
figure(5)
subplot(2,1,1)
plot(t,USB);
xlabel('Time (seconds)')
ylabel('Amplitude (Volt or Ampere)')

USB_fft=fft(USB); 
USB_fft_magnitude=abs(USB_fft)./L; 
USB_fft_single=USB_fft_magnitude(1:(0.5*L)+1); 
USB_fft_single(2:(0.5*L)+1)=2*USB_fft_single(2:(0.5*L)+1); 
axis_single=(0:L/2)*(Fs/L); 

subplot(2,1,2)
plot(axis_single,USB_fft_single)
xlabel('Frequency (Hz)')
ylabel('Magnitude (Volt or Ampere)')

%% Amplification 1
amp = 10;
AMP = USB * amp;
%% GAMBAR 6 : PLOT Hasil filter + Ampli
figure(6)
subplot(2,1,1)
plot(t,AMP);
xlabel('Time (seconds)')
ylabel('Amplitude (Volt or Ampere)')

AMP_fft=fft(USB); 
AMP_fft_magnitude=abs(AMP_fft)./L; 
AMP_fft_single=AMP_fft_magnitude(1:(0.5*L)+1); 
AMP_fft_single(2:(0.5*L)+1)=2*AMP_fft_single(2:(0.5*L)+1); 
axis_single=(0:L/2)*(Fs/L); 

subplot(2,1,2)
plot(axis_single,AMP_fft_single)
xlabel('Frequency (Hz)')
ylabel('Magnitude (Volt or Ampere)')
%% Transmission Line
% Variables
RG59_atten  = 17.38*0.1; 
RG6_atten   = 13.97*0.1;
RG11_atten  = 8.69*0.1;
L1 = 12;
L2 = 4;
L3 = 2;
ZL1 = 33;
ZL2 = 30;
ZL3 = 25;
Z0  = 10;
%% Formula losses
total_loss_cable = 0.5*RG59_atten*L1 + 0.25*RG6_atten*L2 + 0.25*RG11_atten*L2;
Ref  = abs((0.5*ZL1+0.25*ZL2+0.25*ZL3-Z0)/(0.5*ZL1+0.25*ZL2+0.25*ZL3+Z0));
Ref_dB = 10*log10(1/(1-(abs(Ref)^2)));
total_loss_transmission = 10^(-0.1*(total_loss_cable+Ref_dB));
%Formula isyarat setelah lewat transmisi line
v_t = AMP;
delf(n+1) = 0;
delf(1) = 1;
y_t = conv(v_t, total_loss_transmission.*delf(n+1),'same');
%% GAMBAR 7 : PLOT setelah melewati transmission line
figure(7)
subplot(2,1,1);
plot(t,y_t);
xlabel('Time (seconds)')
ylabel('Amplitude (Volt or Ampere)')

y_fft=fft(y_t); 
y_fft_magnitude=abs(y_fft)./L; 
y_fft_single=y_fft_magnitude(1:(0.5*L)+1); 
y_fft_single(2:(0.5*L)+1)=2*y_fft_single(2:(0.5*L)+1); 
axis_single=(0:L/2)*(Fs/L); 

subplot(2,1,2)
plot(axis_single,y_fft_single);
xlabel('Frequency (Hz)')
ylabel('Magnitude (Volt or Ampere)')
%% Proses demodulasi
mod_t = y_t.*c_t;
%% GAMBAR 8 : PLOT Demodulasi
figure(8)
subplot(2,1,1)
plot(t,mod_t)
xlabel('Time (seconds)')
ylabel('Amplitude (Volt or Ampere)')

mod_fft=fft(mod_t); 
mod_fft_magnitude=abs(mod_fft)./L; 
mod_fft_single=mod_fft_magnitude(1:(0.5*L)+1); 
mod_fft_single(2:(0.5*L)+1)=2*mod_fft_single(2:(0.5*L)+1); 
axis_single=(0:L/2)*(Fs/L);
subplot(2,1,2)
plot(axis_single,mod_fft_single);
xlabel('Frequency (Hz)')
ylabel('Magnitude (Volt or Ampere)')
%% Low Pass Filter
A2   = (theta3/pi).*sinc(theta3.*(n-0.5*N)/pi);
B2   = (0/pi).*sinc(0.*(n-0.5*N)/pi);
lpf_n = A2 - B2;
LPF = conv(mod_t,lpf_n,'same');

%% GAMBAR 9 : PLOT frequency response filter h(n)
figure(9);
[h,w]=freqz(lpf_n,1);
plot(abs(h));
xlabel('Frequency (Hz)')
ylabel('Magnitude (Volt or Ampere)')
%% GAMBAR 10 : PLOT Hasil filter 
figure(10)
subplot(2,1,1)
plot(t,LPF);
xlabel('Time (seconds)')
ylabel('Amplitude (Volt or Ampere)')

LPF_fft=fft(LPF); 
LPF_fft_magnitude=abs(LPF_fft)./L; 
LPF_fft_single=LPF_fft_magnitude(1:(0.5*L)+1); 
LPF_fft_single(2:(0.5*L)+1)=2*LPF_fft_single(2:(0.5*L)+1); 
axis_single=(0:L/2)*(Fs/L); 

subplot(2,1,2)
plot(axis_single,LPF_fft_single)
xlabel('Frequency (Hz)')
ylabel('Magnitude (Volt or Ampere)')
%% Amplification 2
amp = 5;
AMP2 = LPF * amp;
%% GAMBAR 11 : PLOT Hasil filter + Ampli
figure(11)
subplot(2,1,1)
plot(t,AMP2);
xlabel('Time (seconds)')
ylabel('Amplitude (Volt or Ampere)')

AMP2_fft=fft(AMP2); 
AMP2_fft_magnitude=abs(AMP2_fft)./L; 
AMP2_fft_single=AMP2_fft_magnitude(1:(0.5*L)+1); 
AMP2_fft_single(2:(0.5*L)+1)=2*AMP2_fft_single(2:(0.5*L)+1); 
axis_single=(0:L/2)*(Fs/L); 

subplot(2,1,2)
plot(axis_single,AMP2_fft_single)
xlabel('Frequency (Hz)')
ylabel('Magnitude (Volt or Ampere)')