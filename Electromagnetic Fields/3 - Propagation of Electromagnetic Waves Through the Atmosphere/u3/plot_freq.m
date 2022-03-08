function [axis_single, X_fft_single]=plot_freq(x,L,Fs)
 
X_fft_magnitude=abs(fft(x))./L; 
X_fft_single=X_fft_magnitude(1:(0.5*L)+1); 
X_fft_single(2:(0.5*L)+1)=2*X_fft_single(2:(0.5*L)+1); 
axis_single=(0:L/2)*(Fs/L);
subplot(2,1,2)
gambar = plot(axis_single,X_fft_single);
xlabel('Frequency (Hz)')
ylabel('Magnitude')