%Lunar Surface Radiation Generator
%Aksoy
%2019

%clear; clc;

%% Setup

% Radiometer Frequency and Bandwidth
freq = 1e9; %in Hz
bw = 1e8; %in MHz

% Lunar Parameters
sunincangle = 180; %Sun Incidence Angle in Degrees, 180 => Night, 0=> Day
thickness = -5; %Regolith Thickness inmeters
layerthickness = 0.01; % Thickness of each regolith layer in meters
pch = 10; %Percentage of FeO and TiO2 amount in the regolith
h = 20; %Regolith Densification parameter
Q = 0; %Geothermal heatflux gradient
permittivitymodel = 2; %1 with frequency dependence, 2 without it, 3 Keihm(1984)

% Outputs
TBrockbase = 0; %Surface Brightness Temperature (EM emission) when the bedrock is rock
TBicebase = 0; %Surface Brightness Temperature (EM emission) when the bedrock is ice
Regcont = 0; %Contribution of Regolith to the Surface Brightness Temperature (EM emission)
Rockcont = 0; %Contribution of Rock base to the Surface Brightness Temperature (EM emission) when the bedrock is rock
Icecont = 0; %Contribution of Ice base to the Surface Brightness Temperature (EM emission) when the bedrock is ice
                
%% Simulation

% Regolith as layered media
z = (0 : -layerthickness : thickness(end)).'; %in meters

% Density of Regolith versus depth (g/cm3)
density = 1.92 - (1.92-1.30).*exp(z/(h*0.01));

% Physical Temperature Profile versus depth (K)
To = 250; % Temperature in the isothermal region
dTo = 162; % Fluctuation amplitude
w = 2*pi/(30*24*3600);
c = 0.55;
k = 8.25e-5 + 7.65e-5*exp(0.25*(2+z*100));
alpha = 2*k./(density*c);
heatflux = Q;
t = 3600*24*30/4 + (sunincangle/360)*3600*24*30;
T = To + dTo.*exp(z*100.*sqrt(w./alpha)).*sin(w*t+z*100.*sqrt(w./alpha)) + -1*z + (-z*(heatflux/10)/0.75); %Damped Sinusoidal close to the surface, linear increase in deep regolith

% Electrical Properties of Regolith, Rock and Ice
a1re = 0.0272/1e9;
a2re = 0.2967;
bre = 0.027;
cre = 3.058;
dre = 0.27;
theta = 300./T(end) - 1;
alpha = (0.00504+0.0062*theta).*exp(-22.1*theta);
beta = (0.502+0.131*theta).*1e-4./(1+theta) + 0.542e-6*((1+theta)./(theta+0.0073)).^2;
eps_r_re = 10.^(dre*density); % Real part of the regolith permittivity
eps_r_ic = 3.1884 + 9.1e-4*(T(end)-273); % Real part of the ice permittivity
eps_ro_ap15 = 5.87 + 1i*0.0086; % Permittivity of the Apollo 15 site bedrock

eps_i_re = zeros(length(density),length(freq)); % Imaginary part of the regolith permittivity
eps_i_ic = zeros(1,length(freq)); % Imginary part of the ice permittivity
pdre = zeros(length(density),length(freq)); % Penetration depth in regolith in meters
abscoeffre = zeros(length(density),length(freq)); % Absorption coefficient in regolith


for tt = 1 : length(density)
    for ff = 1 : length(freq)
        
        if permittivitymodel == 1
            eps_i_re(tt,ff) = eps_r_re(tt)*10.^((a1re*freq(ff)+a2re)*density(tt)+bre*pch-cre);
        elseif permittivitymodel == 2
            eps_i_re(tt,ff) = eps_r_re(tt)*10.^(0.312*density(tt)+0.038*pch-3.260);
        else
            eps_i_re(tt,ff) = 6.079e-5 * density(tt) * eps_r_re(tt) * (freq(ff)*1e-3)^0.25;
        end
        
        eps_i_ic(ff) = (alpha/(freq(ff)*1e-9)) + beta*freq(ff)*1e-9;
        
        omega = 2*pi*freq(ff);
        perm = 4*pi*1e-7;
        permit = 8.854*1e-12;
        
        abscoeffre(tt,ff) = -1*imag(2*pi*freq(ff)*sqrt(perm*permit*(eps_r_re(tt)-1i*eps_i_re(tt,ff))))*2;
        
    end
end

pdre(tt,ff) = (1/(2*abscoeffre(tt,ff))); %in m

%Fresnel Coefficients for each layer
nre = sqrt(repmat(eps_r_re,1,length(freq))+1i*eps_i_re);
n_air = 1;
n_rock = sqrt(eps_ro_ap15);
n_ice = sqrt(eps_r_ic+1i*eps_i_ic);
ns = [n_air*ones(1,length(freq)); nre];

% transmission Coefficients at layer interfaces
trans = zeros(length(density),length(freq));
trans_rock_to_reg = 1-abs((n_rock-ns(end,:))./(n_rock+ns(end,:))).^2;
trans_ice_to_reg = 1-abs((n_ice-ns(end,:))./(n_ice+ns(end,:))).^2;

weightre = zeros(length(density),length(freq));

% Weight functions
for tt = 1 : length(density)
    for ff = 1 : length(freq)
        
        trans(tt,ff) = 1-abs((ns(tt,ff)-ns(tt+1,ff))./(ns(tt,ff)+ns(tt+1,ff))).^2;
        
        % Effect of surface roughness
        if tt == 1
            roughness = 0;
            h = 4*(2*roughness*pi/(3e8/freq(ff)))^2;
            trans(tt,ff) = 1-(abs((ns(tt,ff)-ns(tt+1,ff))./(ns(tt,ff)+ns(tt+1,ff))).^2)*exp(-1*h*eps_r_re(tt));
        end
        
        weightre(tt,ff) = abscoeffre(tt,ff).*exp(-1*sum(abscoeffre(tt:-1:1,ff).*layerthickness))...
            *prod(trans(tt:-1:1,ff));
        
    end
end

%% Simulated Brightness Temperature

contribution_reg = (T.'*weightre*layerthickness);
contribution_rockbase = (T(end)*trans_rock_to_reg).*exp(-1*sum(abscoeffre)*layerthickness).*prod(trans,1);
contribution_icebase = (T(end)*trans_ice_to_reg).*exp(-1*sum(abscoeffre)*layerthickness).*prod(trans,1);

TBrockbase = contribution_reg + contribution_rockbase;
TBicebase = contribution_reg + contribution_icebase;
Regcont = contribution_reg;
Rockcont = contribution_rockbase;
Icecont = contribution_icebase;


%% Corresponding Radiometer Input Voltages (Assuming mean radiometer gain and offset are 1 K/V and 0 K respectively, System Temperature is 0 K, and 100% Antenna Efficiency)
duration = timetorun; %seconds
samplingrate = sample_freq; % variable
vi_rockbase = sqrt(exprnd(TBrockbase, 1, round(duration*samplingrate)));
vi_icebase = sqrt(exprnd(TBicebase, 1, round(duration*samplingrate)));
%% Added this line ~ Ian
v_in = vi_icebase/50;
