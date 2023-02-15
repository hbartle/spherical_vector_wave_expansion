%% Basic vector far-field spherical wave expansion
% Implementation of the vectorial far-field spherical wave expansion
% following FEKO nomenclature and definiton.
%
% Hannes Bartle
% EPFL Microwaves and Antennas Group
% 2023

clear all
close all
clc

f = 9.5e9;
c = physconst("lightspeed");
lambda = c/f;
k = 2*pi/lambda;
eps0 = 8.854*1e-12;
mu0 = 4*pi*1e-7;
Z0 =sqrt(mu0/eps0);

example_folder = "examples/";
%% Read FEKO SWE coefficients
file_path = "rect_patch.out";
%file_path = "RSW_SAR.out";
% file_path = "SWE_surface_wave_modes.out";
% file_path = "SWE_surface_wave_modes_j_1.out";
% file_path = "SWE_surface_wave_modes_j_1_6.out";
% file_path = "SWE_surface_wave_modes_j_6.out";
% file_path = "SWE_surface_wave_modes_j_1_6_21.out";
% file_path = "SWE_surface_wave_modes_j_21.out";
%file_path = "SWE_surface_wave_modes_j_1_6_21_38.out";
%file_path = "SWE_surface_wave_modes_multimode_random_mag_phase.out";

% FEKO provides the SWE coefficients in the ".out" file of the simulation.
% Some annoying manipulation is needed to read the correct lines from this
% file. 
lines = strip(readlines(example_folder+file_path));
modes_start = find(~cellfun(@isempty,strfind(lines,"FAR FIELD MODAL COEFFICIENTS"))) + 10;
modes_end = find(cellfun(@isempty,strfind(lines(modes_start:end)," ")));
modes_end = modes_start +modes_end(1)-2;
mode_data = str2double(split(lines(modes_start:modes_end,1)));    
% In case there's only one mode, this ugly bit of code is needed to make
% Matlab stop complaining. There's probably a better syntax for this.
if size(mode_data,2) > 1 
    j_feko_tmp = mode_data(:,1);
    s_feko_tmp = mode_data(:,2);
    m_feko_tmp = mode_data(:,3);
    n_feko_tmp = mode_data(:,4);
    q_feko_tmp = mode_data(:,5).*exp(1i*mode_data(:,6)*pi/180);
else
    j_feko_tmp = mode_data(1);
    s_feko_tmp = mode_data(2);
    m_feko_tmp = mode_data(3);
    n_feko_tmp = mode_data(4);
    q_feko_tmp = mode_data(5).*exp(1i*mode_data(6)*pi/180);
end


% FEFO removes modes with negligible power from the export so we pad the
% array with them.
for n_tmp = 1:max(n_feko_tmp)
    for m_tmp = -n_tmp:n_tmp
        for s_tmp = 1:2
            j_tmp = 2*(n_tmp*(n_tmp+1)+ m_tmp -1) +s_tmp;
            n_feko(j_tmp) = n_tmp;
            m_feko(j_tmp) = m_tmp;
            s_feko(j_tmp) = s_tmp;
            if sum(s_feko_tmp == s_tmp & m_feko_tmp == m_tmp & n_feko_tmp == n_tmp)
                q_feko(j_tmp) = q_feko_tmp(s_feko_tmp == s_tmp & m_feko_tmp == m_tmp & n_feko_tmp == n_tmp);
            else
                q_feko(j_tmp) = 0;
            end
        end
    end
end


%% Read FEKO rE far-field data export
file_path = "rect_patch_FarField1.ffe";
%file_path = "RSW_SAR_FarField1.ffe";
% file_path = "SWE_surface_wave_modes_FarField1.ffe";
% file_path = "SWE_surface_wave_modes_FarField1_j_1.ffe";
% file_path = "SWE_surface_wave_modes_FarField1_j_1_6.ffe";
% file_path = "SWE_surface_wave_modes_FarField1_j_6.ffe";
% file_path = "SWE_surface_wave_modes_FarField1_j_1_6_21.ffe";
% file_path = "SWE_surface_wave_modes_FarField1_j_21.ffe";
% file_path = "SWE_surface_wave_modes_FarField1_j_1_6_21_38.ffe";
%file_path = "SWE_surface_wave_modes_FarField1_multimode_random_mag_phase.ffe";

E_ff=readmatrix(example_folder+file_path,'FileType','text');
Eff_theta = E_ff(:,3)+1i*E_ff(:,4);
Eff_phi = E_ff(:,5)+1i*E_ff(:,6);
theta = E_ff(:,1)*pi/180;
phi = E_ff(:,2)*pi/180;


%% Read HFSS rE far-field data export
% file_path = "rE_RSW_SAR_singleprobe_HFSS.csv";
% E_ff=readmatrix(example_folder+file_path);
% Eff_theta = E_ff(:,4)+1i*E_ff(:,5);
% Eff_phi = E_ff(:,6)+1i*E_ff(:,7);
% theta = E_ff(:,3)*pi/180;
% phi = E_ff(:,2)*pi/180;



%% Remove special point at theta = pi. Check the m*Pmn/sin(theta)!
Eff_theta( theta==pi) =[];
Eff_phi( theta==pi) =[];
phi( theta==pi) =[];
theta( theta==pi) =[];

%% SWE

q_own = vectorSWE_FF(Eff_theta,Eff_phi,theta,phi,pi/180,pi/180,length(q_feko));

figure
subplot(2,1,1)
bar(abs([[q_feko,zeros(length(q_own)-length(q_feko),1)']',q_own']))
legend(["Feko","Own"])
grid on
xlabel("Mode index j")
ylabel("Magnitude [\surdW]")

subplot(2,1,2)
bar(angle([[q_feko,zeros(length(q_own)-length(q_feko),1)']',q_own'])*180/pi)
legend(["Feko","Own"])          
grid on
xlabel("Mode index j")
ylabel("Phase [°]")

%% Compare power contained in all modes
power_feko = 0.5*sum(abs(q_feko).^2)
power_own = 0.5*sum(abs(q_own).^2)



%% Reconstruct E-field
close all

[Eff_theta_reconstructed, Eff_phi_reconstructed] = vectorSWS_FF(q_own,theta,phi);

% Plot Far-field Cuts
figure
phi_cut = 3*pi/2;
plot(theta(phi==phi_cut)*180/pi,20*log10(abs(sqrt(Eff_theta(phi==phi_cut).^2 + Eff_phi(phi==phi_cut).^2))))
%plot(theta(phi==phi_cut)*180/pi,20*log10(abs(Eff_theta(phi==phi_cut))))
%plot(theta(phi==phi_cut)*180/pi,20*log10(abs(Eff_phi(phi==phi_cut))))
hold on
grid on
plot(theta(phi==phi_cut)*180/pi,20*log10(abs(sqrt(Eff_theta_reconstructed(phi==phi_cut).^2 + Eff_phi_reconstructed(phi==phi_cut).^2))))
%plot(theta(phi==phi_cut)*180/pi,20*log10(abs(Eff_theta_reconstructed(phi==phi_cut))))
%plot(theta(phi==phi_cut)*180/pi,20*log10(abs(Eff_phi_reconstructed(phi==phi_cut))))
legend(["FEKO","Own"])
xlabel('Theta [deg]')
ylabel('Magnitude [dBV]')
title("E_{abs}, \phi=" + phi_cut*180/pi)


figure
subplot(1,2,1)
%patternCustom(abs(Eff_theta),theta*180/pi,phi*180/pi)
%patternCustom(abs(Eff_phi),theta*180/pi,phi*180/pi)
patternCustom(20*log10(abs(sqrt(Eff_theta.^2 + Eff_phi.^2))),theta*180/pi,phi*180/pi)
title("Feko")
subplot(1,2,2)
%patternCustom(abs(Eff_theta_reconstructed),theta*180/pi,phi*180/pi)
%patternCustom(abs(Eff_phi_reconstructed),theta*180/pi,phi*180/pi)
patternCustom(20*log10(abs(sqrt(Eff_theta_reconstructed.^2 + Eff_phi_reconstructed.^2))),theta*180/pi,phi*180/pi)
title("Reconstructed")

% Quantify error
E_abs = abs(sqrt(Eff_theta.^2 + Eff_phi.^2));
E_abs_reconstructed = abs(sqrt(Eff_theta_reconstructed.^2 + Eff_phi_reconstructed.^2));

%E_abs = E_abs/max(E_abs);
%E_abs_reconstructed = E_abs_reconstructed/max(E_abs_reconstructed);

delta_rms = sqrt(1/length(E_abs)* sum(abs(E_abs - E_abs_reconstructed).^2 ))

sgtitle(["Pattern Comparison","\Delta_{rms}="+delta_rms])


