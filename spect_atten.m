function [PS, f] = spect_atten(varargin)
if ischar(varargin{1})
    name = varargin{1};
    [pa0, Ts]= Load_Snapshot(name, 'l');
    if nargin == 3
        flo = varargin{2};
        fhi = varargin{3};
        flag_cutoff = 0;
    else
        freq_cutoff = varargin{2};
        flag_cutoff = 1;
    end
else
    pa0 = varargin{1};
    Ts = varargin{2};
    if nargin == 4
        flo = varargin{3};
        fhi = varargin{4};
        flag_cutoff = 0;
    else
        freq_cutoff = varargin{3};
        flag_cutoff = 1;
    end 
end

pa = pa0;
pb = pa0;
N = 8192;   
df = 1/Ts/N;    
PS_allfreq = abs(fft(pa,N));  
PS_allfreq_b = abs(fft(pb,N));  %
PS_half = PS_allfreq(1:round(end/3)); 
PS_half_b = PS_allfreq_b(1:round(end/3));%
if flag_cutoff
    C = 10^(freq_cutoff/20);
    [ps_max, i_max] = max(PS_half);
    f1 = (find(PS_half > ps_max*C, 1, 'first')-1)*df;
    f2 = (find(PS_half > ps_max*C, 1, 'last')-1)*df;
    f1 = f1*1e-6;
    f2 = f2*1e-6;
else
    f1 = flo;    f2 = fhi;  % MHz
end

% use czt below
N = 128*8;
f = linspace(f1, f2, N);  f = f.';
fs = 1/Ts*1e-6;             %MHz
PS = abs(zfft(pa,f1, f2, fs, N));  
end
