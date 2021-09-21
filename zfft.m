function [y] = zfft(signal, f0, f1, fs, N)
f0_norm = 2*f0/fs;
f1_norm = 2*f1/fs;
y = czt(signal, N, exp(-1i*pi*(f1_norm-f0_norm)/(N-1)), exp(1i*pi*f0_norm)); 
end