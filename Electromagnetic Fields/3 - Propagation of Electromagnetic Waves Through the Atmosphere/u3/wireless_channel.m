%% Code by Gabriel, Electrical Engineering Department
%% Gadjah Mada University  
%% Simple Wireless Channel Function (Exp decay)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [H,channel_axis,H_freq]=wireless_channel(L,Fs)
channel_length=1:L;
channel_axis=(0:L/2)*(Fs/L);
total_path_loss=0.5;
loss_coeff=total_path_loss./(sum(exp(-channel_length/8)));
path_loss=loss_coeff*exp(-channel_length/8); %path loss in power, exponential decay model
H=zeros(1,length(channel_length)); %channel impulse response
        for j=1:length(channel_length)
              H(j)=path_loss(j).*(randn(1,1));
        end
H=abs(H);
H_freq=abs(fft(H))/L;
H_freq=H_freq(1:(0.5*L)+1);
H_freq(2:(0.5*L)+1)=2*H_freq(2:(0.5*L)+1);