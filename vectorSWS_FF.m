function [Eff_theta,Eff_phi] = vectorSWS_FF(q,theta,phi)
%VECTORSWS_FF Far-field vectorial spherical wave synthesis
%   Implementation of the far-field spherical pattern functions follows
%   closely the one implemented in FEKO. 
%
%   For reference, see:
%   https://2021.help.altair.com/2021.1/feko/topics/feko/user_guide/appendix/editfeko_cards/control/card_as_feko_r.htm
%
%   The handling of the special points at theta=0,pi follows this
%   derivation:
%   Sutinjo, Adrian T.. “Calculating Far-Field Radiation Based on FEKO Spherical Wave Coefficients.” (2015).
%
%
%   Input: 
%       q       :   Spherical wave expansion coefficients. Dimensions 1xJ.
%                   J is the number of modes in compressed notation.
%                   
%       theta   :   Theta vector where E-field shall be reconstructed. 
%                   Dimension 1xT*P where T is the number of Theta points 
%                   and P is the number of Phi points.
%       phi     :   Phi vector. Dimension 1xT*P.
%
%
%   Output:
%       Eff_theta   :   Far-field E-field Theta component calculated through SWS.
%                       Dimension 1xT*P
%       Eff_phi     :   Far-field E-field Phi component calculated through SWS.
%                       Dimension 1xT*P
%
%   Hannes Bartle
%   EPFL Microwaves and Antennas Group
%   2023

eps0 = 8.854*1e-12;
mu0 = 4*pi*1e-7;
Z0 =sqrt(mu0/eps0);

J = length(q);
N = -1 + sqrt(1+J/2);
for n = 1:N

    M = repmat(0:n,[length(theta),1]);
    P = legendre(n,cos(theta))';

    % Associated Legendre function divided by sin(theta)
    Ps = P./sin(theta);

    % Handle Special points
    Ps(:,M(1,:)==0) = 0;
    
    l = 0:floor(n/2);
    Ps(theta==0,abs(M(1,:))==1) = -M(theta==0,abs(M(1,:))==1) *...
                                   sum((-1).^l .* factorial(2*n-2*l).*(n-2*l)...
                                   ./(2^n*factorial(l).*factorial(n-l)...
                                   .*factorial(n-2*l)).*(1).^(n-2*l-1));
%     Ps(theta==pi,abs(M(1,:))==1) = -M(theta==0,abs(M(1,:))==1) *...
%                                     sum((-1).^l .* factorial(2*n-2*l).*(n-2*l)...
%                                     ./(2^n*factorial(l).*factorial(n-l)...
%                                     .*factorial(n-2*l)).*(-1).^(n-2*l-1));

    Ps(theta==0,abs(M(1,:)) >=2) = 0;
    Ps(theta==pi,abs(M(1,:))>=2) = 0;
    
    % Derivative of Associated Legendre Function with respect to theta
    dP = cos(theta).*abs(M).*P./sin(theta);
    % Handle Special points
    dP(:,M(1,:)==0) = 0;
    
    l= 0:floor(n/2);
    dP(theta==0,abs(M(1,:))==1) = -1 * sum((-1).^l .* factorial(2*n-2*l)...
                                  .*(n-2*l)./(2^n*factorial(l).*factorial(n-l)...
                                  .*factorial(n-2*l)).*(1).^(n-2*l-1));
%     dP(theta==pi,abs(M(1,:))==1) = -1 * sum((-1).^l .* factorial(2*n-2*l)...
%                                    .*(n-2*l)./(2^n*factorial(l).*factorial(n-l)...
%                                    .*factorial(n-2*l)).*(-1).^(n-2*l-1));


    dP(theta==0,abs(M(1,:)) >=2) = 0;
    dP(theta==pi,abs(M(1,:))>=2) = 0;

    Pshifted = circshift(P,-1,2);
    Pshifted(:,end) = 0;
    
    % Loop through S,M,N to make code more readable
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
            
            j_idx(j) =j;
            s_idx(j) = s;
            m_idx(j) = m;
            n_idx(j) = n;

            if s==1
                K_theta(:,j) = C *(1i)^(n+1).* 1i .*m.*Ps(:,M(1,:)==abs(m));
                K_phi(:,j) = -C*(1i)^(n+1).* (dP(:,M(1,:)==abs(m)) + Pshifted(:,M(1,:) ==abs(m))); 
            elseif s==2
                K_theta(:,j) = C * (1i)^(n).*(dP(:,M(1,:)==abs(m)) + Pshifted(:,M(1,:) ==abs(m)));
                K_phi(:,j) = C * (1i)^(n) .* 1i.* m.*Ps(:,M(1,:)==abs(m));
            end
        end
    end
end

% Amplitude scaling to match input Electric Far-field
A = sqrt(Z0/(2*pi)); 

% Somehow we have q with inverted phase. For reconstruction use complex
% conjugate. This might have something to do with how the complex E-field
% was imported.
mode_cutoff = 1:length(q);
Eff_theta = A*K_theta(:,mode_cutoff)*conj(q');
Eff_phi = A*K_phi(:,mode_cutoff)*conj(q');

end

