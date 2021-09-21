clearvars;
clear all;close all;

% adding the k-wave toolbox path (version 1.3)
addpath(genpath('/home/kjn3/kwave_sims/kWave_1.3'));

% create the computational grid
Nx = 1*128*12*2-2*20;           % number of grid points in the x (axial) direction
Ny = 1*128*4*1-2*20;            % number of grid points in the y (lateral) direction
Nz = 1*128*4*1-2*20;            % number of grid points in the z (azimuth) direction
dx = 21.3e-6;                   % grid point spacing in the x direction [m]
dy = 21.3e-6;                   % grid point spacing in the y direction [m]
dz = 21.3e-6;                   % grid point spacing in the y direction [m]

kgrid = kWaveGrid(Nx, dx, Ny, dy, Nz, dz);

txdc_x = 1.5e-3;                        % [m] transducer position 
txdc_D = 9e-3;                      	% [m] diameter of bowl 
txdc_F = 36e-3;                      	% [m] geometric focus
F_actual = 36e-3;                       % [m] actual focus (distance between transducer and plexi)
saran_L = 64e-6;                      	% [m] saran/membrane thickness
ph_height = 1.5e-3;                   	% [m] distance between plexi and phantom-bot
ph_thickness = 25e-3;               	% [m] thickness of phantom
domain_length = Nx*dx;                  % [m] total length of the simulation domain
plexi_x = txdc_x + F_actual;           	% [m] Plexiglas position (depth)

speed_water = 1480;         		% [m/s]
speed_plexi = 2758;                 % [m/s]
density_water = 1000;       		% [kg/m3]
density_plexi = 1180;               % [kg/m3]
alpha_coeff_water = 0.0162;   		% [dB/cm/MHz^alpha_power_tmm]
alpha_coeff_plexi = 0.000;   		% [dB/cm/MHz^alpha_power_tmm]
alpha_power_tmm = 1.3;         		% [unit-less] Global power term

speed_saran = 2400;         		% [m/s]
speed_tmm = 1540;                   % [m/s]
density_saran = 1690;       		% [kg/m3]
density_tmm = 1040;         		% [kg/m3]
alpha_coeff_saran = 0.7883;    		% [dB/cm/MHz^alpha_power_tmm]
alpha_coeff_tmm = 0.4;      		% [dB/cm/MHz^alpha_power_tmm]

z_water = speed_water*density_water;    % [kg/m2/s]
z_plexi = speed_plexi*density_plexi;    % [kg/m2/s]
z_saran = speed_saran*density_saran;    % [kg/m2/s]
z_tmm = speed_tmm*density_tmm;          % [kg/m2/s]

t_plexi = 2*(F_actual)/speed_water;   	% [s]

% define the properties of the propagation medium - water, plexi
medium.alpha_power = alpha_power_tmm; 

medium.sound_speed = speed_water * ones(Nx, Ny, Nz);                   
medium.sound_speed(round(plexi_x/dx):Nx, :,:) = speed_plexi;           
medium.sound_speed_ref = speed_water;                                  

medium.density = density_water * ones(Nx, Ny, Nz);                    
medium.density(round(plexi_x/dx):Nx, :,:) = density_plexi;             

medium.alpha_coeff = alpha_coeff_water * ones(Nx, Ny, Nz);                               
medium.alpha_coeff(round(plexi_x/dx):Nx, :,:) = alpha_coeff_plexi;           

% create time array
t_end = t_plexi*1.3;                        % [s]
cfl = 0.25;                                 % cfl = PPP/PPW
[kgrid.t_array, dt] = makeTime(kgrid,max([speed_water,speed_plexi]), cfl, t_end);

% define a curved transducer element
source.p_mask = zeros(Nx, Ny, Nz);
bowl_pos = [round(txdc_x/dx), Ny/2+1, Nz/2+1];                      % [grid points]  
radius = round(txdc_F/dx);                                          % [grid points]
diameter = round(txdc_D/dx)+(mod(round(txdc_D/dx),2)==0);           % [grid points]
focus_pos = [ceil((txdc_x+txdc_F)/dx), (Ny/2+1), (Nz/2+1)];         % [grid points]
source.p_mask = makeBowl([Nx, Ny, Nz], bowl_pos, radius, diameter, focus_pos);

%% Source excitation definition 

%txdc properties
f0 = 10e6;              % [Hz]
bw = 0.5;               % [%age] bw*f0 = bandwidth
bwr = -3;               % [dB] 

% Generate transducer impulse response generate an Gaussian pulse @f0 with a
% -3dB bandwith of 50%
tt  = -5/f0 : kgrid.dt/100 : 5/f0;
g_pulse = gauspuls(tt,f0,bw,bwr);

% Time domain integration of source excitation-to counteract time-differential of
% mass source term (Eq. 2.9 - k-Wave user manual version 1.1)
ginteg = dt*cumtrapz(g_pulse);
ginteg_pulse = ginteg(1:100:length(tt));
impulse_response = ginteg_pulse;

% An input kronecker delta of -100V magnitude is the excitation signal 
source_mag = -100;
transmit_gain = 1e6/100;                % [Pa/V] transmit sensitivity of transducer
receive_gain = 100/1e8;                 % [V/Pa] receive sensitivity of transducer
t_arr = zeros(1,length(kgrid.t_array));
t_arr(1) = source_mag;
t_arr2 = conv(t_arr,impulse_response);
t_arr2(length(impulse_response)+1:length(t_arr2))=t_arr2(length(impulse_response));
t_arr = t_arr2(1:length(t_arr));
source.p = transmit_gain*t_arr;

% filter the source to remove any high frequencies not supported by the grid
source.p = filterTimeSeries(kgrid, medium, source.p);
source_p = source.p;

%% Code continued...

% create a display mask to display the transducer
display_mask = source.p_mask;

% assign the input options
input_args = {'PlotLayout', true,'DisplayMask', display_mask, 'DataCast', 'single', 'PMLInside', false, 'PlotPML', false,'PMLSize', 20};

% define sensor as same as source
sensor.mask = source.p_mask;

% set the record mode to capture the final wave-field and the statistics at
% each sensor point
sensor.record = {'p'};

% Code is compiled. Input hdf5 file to k-Wave CPU code is generated. 
% submit_ref.sh can be run on a job submission system with appropriate 
% CPU-core and memory requirements 
filename = './Txdc_10M_F4_3D_ref.h5';
kspaceFirstOrder3D(kgrid, medium, source, sensor, 'SaveToDisk', filename, input_args{:});

% Code parameters are saved for post-processing of the RF-data obtained from
% the simulation
t_step = kgrid.dt;
save Txdc_10M_3D_ref_param Nx Ny dx dy...
    f0 bw bwr t_step source_p txdc_x txdc_D txdc_F F_actual...
    saran_L ph_height ph_thickness...
    speed_water speed_plexi density_water density_plexi...
    alpha_coeff_water alpha_coeff_plexi alpha_power_tmm...
    speed_saran speed_tmm density_saran density_tmm...
    alpha_coeff_saran alpha_coeff_tmm

