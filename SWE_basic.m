%
%   Basic Far-field spherical wave expansion
%
clear all
close all
clc

f = 9.5e9;
c = physconst("lightspeed");
lambda = c/f;
k = 2*pi/lambda;
Z0 = 120*pi;


% FEKO SWE coefficients
file_path = "rect_patch.out";
% file_path = "SWE_surface_wave_modes.out";
% file_path = "SWE_surface_wave_modes_j_1.out";
% file_path = "SWE_surface_wave_modes_j_1_6.out";
% file_path = "SWE_surface_wave_modes_j_6.out";
% file_path = "SWE_surface_wave_modes_j_1_6_21.out";
% file_path = "SWE_surface_wave_modes_j_21.out";
%file_path = "SWE_surface_wave_modes_j_1_6_21_38.out";
%file_path = "SWE_surface_wave_modes_multimode_random_mag_phase.out";


lines = strip(readlines(file_path));
modes_start = find(~cellfun(@isempty,strfind(lines,"FAR FIELD MODAL COEFFICIENTS"))) + 10;
modes_end = find(cellfun(@isempty,strfind(lines(modes_start:end)," ")));
modes_end = modes_start +modes_end(1)-2;
% SWE
mode_data = str2double(split(lines(modes_start:modes_end,1)));    
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

% Add zero magnitude modes that are ommitted by FEKO for easier comparison
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


 file_path = "rect_patch_FarField1.ffe";
% file_path = "SWE_surface_wave_modes_FarField1.ffe";
% file_path = "SWE_surface_wave_modes_FarField1_j_1.ffe";
% file_path = "SWE_surface_wave_modes_FarField1_j_1_6.ffe";
% file_path = "SWE_surface_wave_modes_FarField1_j_6.ffe";
% file_path = "SWE_surface_wave_modes_FarField1_j_1_6_21.ffe";
% file_path = "SWE_surface_wave_modes_FarField1_j_21.ffe";
% file_path = "SWE_surface_wave_modes_FarField1_j_1_6_21_38.ffe";
%file_path = "SWE_surface_wave_modes_FarField1_multimode_random_mag_phase.ffe";

E_ff=readmatrix(file_path,'FileType','text');

Eff_theta = E_ff(:,3)+1i*E_ff(:,4);
Eff_phi = E_ff(:,5)+1i*E_ff(:,6);
theta = E_ff(:,1)*pi/180;
phi = E_ff(:,2)*pi/180;


% Eff_theta(theta==0 | theta==pi) =[];
% Eff_phi(theta==0 | theta==pi) =[];
% phi(theta==0 | theta==pi) =[];
% theta(theta==0 | theta==pi) =[];

% Construct Far-field spherical pattern functions
N=max(n_feko);
for n = 1:N

    M = repmat(0:n,[length(theta),1]);
    P = legendre(n,cos(theta))';

    % Associated Legendre function divided by sin(theta)
    Ps = P./sin(theta);

    % Handle Special points
    Ps(:,M(1,:)==0) = 0;
    
    % Check these values!
    l = 1:floor(n/2);
    Ps(theta==0,abs(M(1,:))==1) = -M(theta==0,abs(M(1,:))==1) *...
                                   sum((-1).^l .* factorial(2*n-2*l).*(n-2*l)...
                                   ./(2^n*factorial(l).*factorial(n-l)...
                                   .*factorial(n-2*l)).*(cos(0)).^(n-2*l-1));
    Ps(theta==pi,abs(M(1,:))==1) = -M(theta==0,abs(M(1,:))==1) *...
                                    sum((-1).^l .* factorial(2*n-2*l).*(n-2*l)...
                                    ./(2^n*factorial(l).*factorial(n-l)...
                                    .*factorial(n-2*l)).*(cos(pi)).^(n-2*l-1));

    Ps(theta==0,abs(M(1,:)) >=2) = 0;
    Ps(theta==pi,abs(M(1,:))>=2) = 0;
    
    % Derivative of Associated Legendre Function with respect to theta
    dP = cos(theta).*abs(M).*P./sin(theta);
    % Handle Special points
    dP(:,M(1,:)==0) = 0;
    
    l= 1:floor(n/2);
    dP(theta==0,abs(M(1,:))==1) = -1 * sum((-1).^l .* factorial(2*n-2*l)...
                                  .*(n-2*l)./(2^n*factorial(l).*factorial(n-l)...
                                  .*factorial(n-2*l)).*(cos(0)).^(n-2*l-1));
    dP(theta==pi,abs(M(1,:))==1) = -1 * sum((-1).^l .* factorial(2*n-2*l)...
                                   .*(n-2*l)./(2^n*factorial(l).*factorial(n-l)...
                                   .*factorial(n-2*l)).*(cos(pi)).^(n-2*l-1));

    dP(theta==0,abs(M(1,:)) >=2) = 0;
    dP(theta==pi,abs(M(1,:))>=2) = 0;

    Pshifted = circshift(P,-1,2);
    Pshifted(:,end) = 0;

    for m = -n:n
        for s = 1:2
            j = 2*(n*(n+1)+ m -1) +s;
            
            if m == 0
                Cm = 1;
            else
                Cm = (-m./abs(m)).^m;
            end

            C = Cm*sqrt(1/(n*(n+1)))*...
                sqrt((2*n+1)/2 *factorial(n-abs(m))./factorial(n+abs(m)))...
                *exp(1i*m*phi);
            
            qj(j) =j;
            s_idx(j) = s;
            m_idx(j) = m;
            n_idx(j) = n;

            A = k*sqrt(Z0);
            A = 1;
            if s==1
                K_theta(:,j) = A*C *(1i)^(n+1).* 1i .*m.*Ps(:,M(1,:)==abs(m));
                K_phi(:,j) = -A*C*(1i)^(n+1).* (dP(:,M(1,:)==abs(m)) + Pshifted(:,M(1,:) ==abs(m))); 
            elseif s==2
                K_theta(:,j) = A*C * (1i)^(n).*(dP(:,M(1,:)==abs(m)) + Pshifted(:,M(1,:) ==abs(m)));
                K_phi(:,j) = A*C * (1i)^(n) .* 1i.* m.*Ps(:,M(1,:)==abs(m));
            end
        end
    end
end

% Get SWE coefficients

dtheta = pi/180;
dphi = pi/180;

q = repmat(Eff_theta,[1,length(qj)]).*conj(K_theta).*repmat(sin(theta),[1,length(qj)])...
    + repmat(Eff_phi,[1,length(qj)]).*conj(K_phi).*repmat(sin(theta),[1,length(qj)]);

for j = 1:length(qj)
    
    Q = reshape(q(:,j),[],length(unique(theta)));
    A = sqrt(2*pi)/(k*sqrt(Z0));
    A = 4*pi^2;
    q_own(j) = 1/A*trapz(pi/180,trapz(pi/180,Q));


    % Filter out low magnitude modes
%     if abs(q_own(j)) < 5e-3
%         q_own(j) = 0;
%     end
end

figure
subplot(2,1,1)
bar(abs([[q_feko,zeros(length(q_own)-length(q_feko),1)']',q_own']))
legend(["Feko","Own"])
grid on
xlabel("Mode index j")
title("Magnitude")

subplot(2,1,2)
bar(angle([[q_feko,zeros(length(q_own)-length(q_feko),1)']',q_own'])*180/pi)
legend(["Feko","Own"])          
grid on
xlabel("Mode index j")
title("Phase")


power_feko = 0.5*sum(abs(q_feko).^2)
power = 0.5*sum(abs(q_own).^2)



% % Check modes for orthogonality
check_mode = 1;
for j= 1:max(qj)

    K1_theta = reshape(K_theta(:,check_mode),[],length(unique(theta)));
    K1_phi = reshape(K_phi(:,check_mode),[],length(unique(theta)));
    K2_theta = reshape(K_theta(:,j),[],length(unique(theta)));
    K2_phi = reshape(K_phi(:,j),[],length(unique(theta)));
    
    sint = reshape(sin(theta),[],length(unique(theta)));

    orth_theta(j) = abs(trapz(pi/180,trapz(pi/180,K1_theta.*conj(K2_theta).*sint)));
    orth_phi(j) = abs(trapz(pi/180,trapz(pi/180,K1_phi.*conj(K2_phi).*sint)));
end
%[qj',s_idx',n_idx',m_idx',round(orth_theta'),round(orth_phi')]

%% Reconstruct E-field
close all



mode_cutoff = length(q_own);
%mode_cutoff = 20

A = k*sqrt(Z0/(2*pi));

% Somehow we have q with inverted phase. For reconstruction use complex
% conjugate
Eff_theta_reconstructed = A*K_theta(:,1:mode_cutoff)*conj(q_feko(1:mode_cutoff)');
Eff_phi_reconstructed = A*K_phi(:,1:mode_cutoff)*conj(q_feko(1:mode_cutoff)');
Eff_theta_reconstructed = A*K_theta(:,1:mode_cutoff)*conj(q_own(1:mode_cutoff)');
Eff_phi_reconstructed = A*K_phi(:,1:mode_cutoff)*conj(q_own(1:mode_cutoff)');

% Plot Far-field Cuts
figure
phi_cut = pi/2;
plot(theta(phi==phi_cut)*180/pi,10*log10(abs(sqrt(Eff_theta(phi==phi_cut).^2 + Eff_phi(phi==phi_cut).^2))))
%plot(theta(phi==phi_cut)*180/pi,10*log10(abs(Eff_theta(phi==phi_cut))))
%plot(theta(phi==phi_cut)*180/pi,10*log10(abs(Eff_phi(phi==phi_cut))))
hold on
grid on
plot(theta(phi==phi_cut)*180/pi,10*log10(abs(sqrt(Eff_theta_reconstructed(phi==phi_cut).^2 + Eff_phi_reconstructed(phi==phi_cut).^2))))
%plot(theta(phi==phi_cut)*180/pi,10*log10(abs(Eff_theta_reconstructed(phi==phi_cut))))
%plot(theta(phi==phi_cut)*180/pi,10*log10(abs(Eff_phi_reconstructed(phi==phi_cut))))
legend(["FEKO","Own"])
xlabel('Theta [deg]')

figure
subplot(1,2,1)
%patternCustom(abs(Eff_theta),theta*180/pi,phi*180/pi)
%patternCustom(abs(Eff_phi),theta*180/pi,phi*180/pi)
patternCustom(abs(sqrt(Eff_theta.^2 + Eff_phi.^2)),theta*180/pi,phi*180/pi)
title("Feko")
subplot(1,2,2)
%patternCustom(abs(Eff_theta_reconstructed),theta*180/pi,phi*180/pi)
%patternCustom(abs(Eff_phi_reconstructed),theta*180/pi,phi*180/pi)
patternCustom(abs(sqrt(Eff_theta_reconstructed.^2 + Eff_phi_reconstructed.^2)),theta*180/pi,phi*180/pi)
title("Reconstructed")


% Quantify error
E_abs = abs(sqrt(Eff_theta.^2 + Eff_phi.^2));
E_abs_reconstructed = abs(sqrt(Eff_theta_reconstructed.^2 + Eff_phi_reconstructed.^2));

E_abs = E_abs/max(E_abs);
E_abs_reconstructed = E_abs_reconstructed/max(E_abs_reconstructed);

delta_rms = sqrt(1/length(E_abs)* sum(abs(E_abs - E_abs_reconstructed).^2 ))
%% Generate FF from any mode
close all



q = zeros(1,6);
q(1) = 1*exp(1j*0);
q(5) = 0.8*exp(1j*0);
 q(17) = 0.6*exp(1j*-20*pi/180);
q(20) = 0.5*exp(1j*-20*pi/180);
q(30) = 0.4*exp(1j*-50*pi/180);
q(44) = 0.3*exp(1j*-50*pi/180);
q(55) = 0.1*exp(1j*-90*pi/180);
q(56) = 0.05*exp(1j*-90*pi/180);
q(80) = 0.1*exp(1j*-180*pi/180);

mode = 1:length(q);

Eff_theta_reconstructed = K_theta(:,mode)*q';
Eff_phi_reconstructed = K_phi(:,mode)*q';


figure
subplot(1,2,1)
patternCustom(abs(Eff_phi),theta*180/pi,phi*180/pi)
title("Feko")
subplot(1,2,2)
patternCustom(abs(Eff_phi_reconstructed),theta*180/pi,phi*180/pi)
title("Reconstructed")
sgtitle('E_{\phi}')

figure
subplot(1,2,1)
patternCustom(abs(Eff_theta),theta*180/pi,phi*180/pi)
title("Feko")
subplot(1,2,2)
patternCustom(abs(Eff_theta_reconstructed),theta*180/pi,phi*180/pi)
title("Reconstructed")
sgtitle('E_{\theta}')


figure
subplot(1,2,1)
patternCustom(abs(sqrt(Eff_theta.^2 + Eff_phi.^2)),theta*180/pi,phi*180/pi)
title("Feko")
subplot(1,2,2)
patternCustom(abs(sqrt(Eff_theta_reconstructed.^2 + Eff_phi_reconstructed.^2)),theta*180/pi,phi*180/pi)
title("Reconstructed")
sgtitle('E_{abs}')