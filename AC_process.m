%%%%%%%%%%%%%%% Modify these fields and run the script %%%%%%%%%%%

clc; clear all; close all;

refName = 'Txdc_10M_F4_3D_ref_out.h5';
attName1_top = 'Txdc_10M_F4_3D_att_top_out.h5';
sosName = 'Txdc_10M_F4_3D_att_out.h5';
attName = sosName;

% Refer to Figure 2 in paper
sos_tstart = 10;            % [us] start of echoes - E1+E2+E3
sos_tend = 49;              % [us] end of echoes - E1+E2+E3
att_tstart = 47;            % [us] start of echo E3  
att_tend = 49;              % [us] end of echo E3
ref_tstart = 48;            % [us] start of echo Er
ref_tend = 51;              % [us] end of echo Er
att_top_tstart = ref_tstart;% [us] start of echo Em
att_top_tend = ref_tend;    % [us] end of echo Em

%% System parameters
cutoff = -30; 
temperature = 20;       % Temperature of the setup
gain_ref = 0;           % Gain of receive chain for Setup1
gain_att = 0;           % Gain of receive chain for Setup2
gain_att_top = 0;       % Gain of receive chain for Setup3

load('Txdc_10M_3D_ref_param.mat');

f0 = 10e6;              % [Hz]
bw = 0.5;               % bw*f0 = bandwidth
bwr = -3;               % [dB] Freq cutoff 
transmit_gain = 1e6/100;
receive_gain = 100/1e8;
tt  = -5/f0 : t_step : 5/f0;
g_pulse = gauspuls(tt,f0,bw,bwr);

t_step = dx/speed_plexi/4;
f_nyq = speed_water/dx/1e6;     % [MHz] Spatial-nyquist frequency for the simulation setup
freq_low = f_nyq/20;            % [MHz] Lower threshold for AC estimation (refer Figure 10)
freq_high = f_nyq/5;            % [MHz] Upper threshold for AC estimation (refer Figure 10)

% In an experimental context, lower and upper frequency threshold for  
% AC estimation is chosen as the intersection of -20dB bandwidths of
% echoes, Er and E3.

%% loading sos data 
sense = h5read(sosName, '/p');
sensor_data_p = sum(sensor);
% load(sosName);    % use this if rf-data are in *.mat format

sos_rf = receive_gain*conv(sensor_data_p,g_pulse);
sos_time = t_step*(0:length(sos_rf)-1)*1e6;
figure(4);hold on;
plot(sos_time,sos_rf/max(sos_rf),'DisplayName','Setup-2');
xlabel('Time (us)');ylabel('Voltage (V)');
title('Received A-line');
hold off;
hilbert_focus = abs(hilbert(sos_rf));
figure(5);hold on;plot(sos_time,20*log10(hilbert_focus/max(hilbert_focus)),'DisplayName','Setup-2');
xlabel('Time (us)');ylabel('Normalized-pressure (dB)');
title('Sensor pressure envelope');
hold off;

%% loading att_top data 
sense = h5read(attName1_top, '/p');
sensor_data_p = sum(sensor);
% load(attName1_top);    % use this if rf-data are in *.mat format

att_top_rf = receive_gain*conv(sensor_data_p,g_pulse);
att_top_time = t_step*(0:length(att_top_rf)-1)*1e6;
figure(4);hold on;
plot(att_top_time,att_top_rf/max(att_top_rf),'DisplayName','Setup-3');
xlabel('Time (us)');ylabel('Voltage (V)');
title('Received A-line');
hold off;
hilbert_focus = abs(hilbert(att_top_rf));
figure(5);hold on;plot(att_top_time,20*log10(hilbert_focus/max(hilbert_focus)),'DisplayName','Setup-3');
xlabel('Time (us)');ylabel('Normalized-pressure (dB)');
title('Sensor pressure envelope');
hold off;

%% loading ref data 
sense = h5read(refName, '/p');
sensor_data_p = sum(sensor);
% load(refName);    % use this if rf-data are in *.mat format

ref_rf = receive_gain*conv(sensor_data_p,g_pulse);
ref_time = t_step*(0:length(ref_rf)-1)*1e6;
figure(4);hold on;
plot(ref_time,ref_rf/max(ref_rf),'DisplayName','Setup-1');
xlabel('Time (us)');ylabel('Voltage (V)');
title('Received A-line');
hold off;
hilbert_focus = abs(hilbert(ref_rf));
figure(5);hold on;plot(ref_time,20*log10(hilbert_focus/max(hilbert_focus)),'DisplayName','Setup-1');
xlabel('Time (us)');ylabel('Normalized-pressure (dB)');
title('Sensor pressure envelope');
hold off;

%% Formatting echo data - snipping echo and computing spectrum

[y_ref,t_ref] = snip_array(ref_rf,ref_time,ref_tstart,ref_tend);        % snipt Er
[r_ref,f_ref] = spect_atten(y_ref,t_step,freq_low,freq_high);           % compute Er-spectrum

[y_sos, t_sos] = snip_array(sos_rf,sos_time,sos_tstart,sos_tend);       % snipt E1+E2+E3
[y_att, t_att] = snip_array(sos_rf,sos_time,att_tstart,att_tend);       % snip E3
[r_att,f_att] = spect_atten(y_att,t_step,freq_low,freq_high);           % compute E3-spectrum

[y_att_top, t_att_top] = snip_array(att_top_rf,att_top_time,att_top_tstart,att_top_tend);   % snipt Em
[r_att_top,f_att_top] = spect_atten(y_att_top,t_step,freq_low,freq_high);                   % compute Em-spectrum


%% Calculating phantom thickness & phantom material sound speed

R2 = 20*log10(r_ref) - max(20*log10(r_ref));
R0 = 20*log10(r_att) - max(20*log10(r_att));
R0_top = 20*log10(r_att_top) - max(20*log10(r_att_top));       

fsize = 12;
figure(1);hold on;
plot(f_ref,R2, f_att, R0,f_att_top,R0_top,'LineWidth',2);
legend('E_r','E_3','E_m','FontSize',fsize);
xlabel('Freq (MHz)','FontSize',fsize);
ylabel('Normalized power (dB)','FontSize',fsize);
title('Received A-line spectrum','FontSize',fsize);
ylim([cutoff 0]);
hold off;

f_low = max(f_ref(1), f_att(1));
f_high = min(f_ref(end), f_att(end));
f_low_membrane = max(f_ref(1), f_att_top(1));
f_high_membrane = min(f_ref(end), f_att_top(end));
[r_ref,f_ref] = spect_atten(y_ref,t_step,f_low,f_high);
[r_att,f_att] = spect_atten(y_att,t_step,f_low,f_high);
[r_att_top,f_att_top] = spect_atten(y_att_top,t_step,f_low,f_high);
[r_ref_membrane,f_ref_membrane] = spect_atten(y_ref,t_step,f_low_membrane,f_high_membrane);
[r_att_top_membrane,f_att_top_membrane] = spect_atten(y_att_top,t_step,f_low_membrane,f_high_membrane);

% applying cross-correlation analysis
[yc lags]= xcorr(y_sos, y_ref);
tc = lags*(t_sos(2)-t_sos(1)) + t_sos(1) - t_ref(1);

%removing unncessary portions from Setup RF-data; retaining E1,E2,E3 
ycNeg = yc(tc<=0);
tcNeg = tc(tc<=0);
psd = pwelch(ycNeg);
[~, freqIdx] = max(psd(10:end));
freqIdx = freqIdx + 9;
extFactor = 4*4; % single side
extLength = round(extFactor/(freqIdx/length(psd)/2));

% peak detection for the RF-data from Setup2 (E1+E2+E3)
peaks = zeros(3,1); locs = zeros(3,1);
[peaks(1), locs(1)] = max(ycNeg);
pk1 = peaks(1);
tpk1 = tcNeg(locs(1));

mask1 = ones(size(ycNeg)); 
mask1(max(1,locs(1)-extLength) : min(locs(1)+extLength, length(mask1))) = 0;
[peaks(2), locs(2)] = max(ycNeg.*mask1);
mask2 = mask1;
mask2(max(1,locs(2)-extLength) : min(locs(2)+extLength, length(mask2))) = 0;
[peaks(3), locs(3)] = max(ycNeg.*mask2);
[locs sortIdx] = sort(locs);
peaks = peaks(sortIdx);
tPeaks = tcNeg(locs);

figure; plot(tc, yc, 'b', tc(tc>0), yc(tc>0), 'r');
hold on; 
plot(tcNeg(locs), peaks, '*k');
title('Cross Correlation between Ref and Sos Signals')
xlabel('Time lag (\mus)')
ylabel('Cross correlation')

cw = 1480;                                  % sound of speed in water (m/s)
%speedinwater(temperature);                 % use for experimental data
cp = cw*(-tPeaks(3) + tPeaks(2)-tPeaks(1))./(tPeaks(2)-tPeaks(1));  %speed of sound in phantom
thickness = 1e-4*cp.*(tPeaks(2)-tPeaks(1))/2;       % in cm (Equation 10 in paper)
height = 1e-4*cw.*(tPeaks(3)-tPeaks(2))/2;          % in cm        

% water attenuation with temperature (Fisher'77)
water_atten = (55.9 - 2.37*temperature + 0.0477*temperature^2- ...
    0.000348*temperature^3)*(1e-15)*(1-3.84e-4)*1e10*(20*log10(exp(1)));

c_all = cp;
z_all = thickness;
psd_r = r_ref;
f = f_ref';
psd_a = r_att;
psd_a = psd_a*10^(-(gain_att-gain_ref)/20);
aw0 = water_atten*f.^2;

psd_a_ph_top = r_att_top;
psd_a_ph_top = psd_a_ph_top*10^(-(gain_att_top-gain_ref)/20);

f_membrane = f_ref_membrane';
psd_r_membrane = r_ref_membrane;
psd_a_ph_top_membrane = r_att_top_membrane;
psd_a_ph_top = psd_a_ph_top*10^(-(gain_att_top-gain_ref)/20);

%Method0 - without compensating for membrane
psd_ratio_0 = psd_a./psd_r; 
atten_0 = -20*log10(psd_ratio_0)/thickness/2 + aw0;     

%% Non linear curve fit - 
rwp_measured = double(psd_a_ph_top_membrane./psd_r_membrane*0.3748);    % Membrane pressure reflection coeff (R1) 
r1 = 1000*cw;                                                           % Impedance of water 
modelfun = @(b)(abs(cos(2*pi.*f_membrane*b(1))*(1-r1/b(3))+1i*sin(2*pi.*f_membrane*b(1))*(b(2)/b(3)-r1/b(2)))./abs(cos(2*pi.*f_membrane*b(1))*(1+r1/b(3))+1i*sin(2*pi.*f_membrane*b(1))*(r1/b(2)+b(2)/b(3))))-rwp_measured;

%beta array - parameters of the membrane
%beta(1) = L/C2, beta(2) = R2, beta(3) = R3
beta0_mean = [20*25/1500;3*r1;3*r1];        
beta0_lb = [1/5000;0.5*r1;0.5*r1];
beta0_ub = [1000/1000;10*r1;10*r1];
beta0_step = [2.5;2.5;2.5];      
iter_max = 10;
mse = ones(iter_max+1,iter_max+1,iter_max+1);
beta = [];

for iter1 = 1:iter_max+1
    for iter2 = 1:iter_max+1
        for iter3 = 1:iter_max+1
            beta0 = beta0_mean.*(beta0_step.^[iter1-iter_max/2;iter2-iter_max/2;iter3-iter_max/2]);
            if((abs((beta0(3)-beta0(2))/(beta0(3)+beta0(2)))>0.01) && (abs((r1-beta0(2))/(r1+beta0(2)))>0.01)&&(beta0(2)>r1)&&(beta0(3)>r1)&&(beta0(1)>0)&& (abs((beta0(2)^2-r1*beta0(3))/(beta0(2)^2+r1*beta0(3)))>0.01))
                [beta,mse(iter1,iter2,iter3)] = lsqnonlin(modelfun,beta0,beta0_lb,beta0_ub);
            else
                 mse(iter1,iter2,iter3) = 100;
            end

        end
    end
end

% choosing the beta values with lowest rms error for curve-fitting R1(f)
[min_mse,linear_ind] = min(mse(:));
[min_mse_loc(1),min_mse_loc(2),min_mse_loc(3)] = ind2sub(size(mse),linear_ind);
beta0 = beta0_mean.*(beta0_step.^[min_mse_loc(1)-iter_max/2;min_mse_loc(2)-iter_max/2;min_mse_loc(3)-iter_max/2]);
[beta,err] = lsqnonlin(modelfun,beta0,beta0_lb,beta0_ub);
reflect_membrane_1 = modelfun(beta)+rwp_measured;

% computing the one-way membrane pressure transmisssion coeffcient from 
% the estimated membrane parameters
transmit_membrane_1 = 4*beta(3)*r1./abs((r1+beta(3))*cos(2*pi.*f_membrane*beta(1))+1i*(beta(2)+r1*beta(3)/beta(2)).*sin(2*pi.*f_membrane*beta(1))).^2;
Transmit = 4*beta(3)*r1./abs((r1+beta(3))*cos(2*pi.*f*beta(1))+1i*(beta(2)+r1*beta(3)/beta(2)).*sin(2*pi.*f*beta(1))).^2;

% removing roundtrip membrane transmission loss from the measured insertion
% loss (Equation 11 in paper)
psd_ratio_new3 = psd_a./psd_r./(Transmit.^2); 
atten_new3 = -20*log10(psd_ratio_new3)./thickness/2+aw0;     
atten_new = atten_new3;                  

thickness = z_all;
thickness_mean = mean(z_all);
speed = c_all;
speed_mean = mean(speed);
atten_mean_0 = mean(atten_0,2);
atten_mean_new = mean(atten_new,2);

%% Membrane pressure reflection and transmission coeff (ground truth)

freq_base = linspace(f_nyq/20*1e6,f_nyq/5*1e6,1024);
l = saran_L;
alpha = log(10)*5*alpha_coeff_saran*(freq_base/1e6).^alpha_power_tmm;
k2 = 2*pi*freq_base/speed_saran-1i*alpha;

z1 = speed_water * density_water;
z2= speed_saran * density_saran;
z3 = speed_tmm * density_tmm;
zp = speed_plexi * density_plexi;

Rplexi = (zp-z1)/(zp+z1);
Rwp = abs((cos(k2*l)*(z3-z1)*z2+1i*sin(k2*l)*(z2*z2-z1*z3))...
    ./(cos(k2*l)*(z1+z3)*z2+1i*sin(k2*l)*(z2*z2+z1*z3)));
Rpw = abs((cos(k2*l)*(z1-z3)*z2+1i*sin(k2*l)*(z2*z2-z1*z3))...
    ./(cos(k2*l)*(z1+z3)*z2+1i*sin(k2*l)*(z2*z2+z1*z3)));

Twp = abs(2*z2*z3./(cos(k2*l)*(z1+z3)*z2+1i*sin(k2*l)*(z2*z2+z1*z3)));
Tpw = abs(2*z2*z1./(cos(k2*l)*(z1+z3)*z2+1i*sin(k2*l)*(z2*z2+z1*z3)));

%% Plotting results

fsize = 12;
figure(10);
hold on;
plot(f,alpha_coeff_tmm*f.^alpha_power_tmm,f,atten_new,f,atten_0,'LineWidth',2);
title ('Attenuation coefficient','FontSize',fsize);
legend({'Ground truth','Estimated','Phantom material+membrane'}, 'location', 'northwest','FontSize',fsize)
xlabel('freq (MHz)','FontSize',fsize);
ylabel('Attenuation coefficient (dB/cm)','FontSize',fsize);
hold off;

figure(11);
hold on;
plot(f,alpha_coeff_tmm*f.^(alpha_power_tmm-1),f,atten_new./f,f,atten_0./f,'LineWidth',2);
title ('Attenuation coefficient slope','FontSize',fsize);
legend({'Equation','Estimated','Phantom material+Membrane'}, 'location', 'northwest','FontSize',fsize)
xlabel('freq (MHz)','FontSize',fsize);
ylabel('Attenuation coefficient (dB/cm-MHz)','FontSize',fsize);
ylim([0.2,1.6])
hold off;

figure(13);hold on;
plot(freq_base/1e6,Rwp,f_ref_membrane,r_att_top_membrane./r_ref_membrane*0.375,'LineWidth',2)
title('Membrane reflection coefficient');
xlabel('Frequency (MHz)');ylabel('Reflection coefficient (V/V)');
legend({'Ground truth','Measured'}, 'location', 'northwest','FontSize',fsize)
ylim([0 1])
hold off;

figure(14);hold on;
plot(freq_base/1e6,((Twp.*Tpw).^2),f_ref_membrane,((transmit_membrane_1.^2)),'LineWidth',2)
title('Round-trip transmission loss - membrane');
xlabel('Frequency (MHz)');ylabel('Round-trip transmission coefficient (V/V)');
legend({'Ground truth','Estimated'}, 'location', 'northwest','FontSize',fsize)
hold off;


