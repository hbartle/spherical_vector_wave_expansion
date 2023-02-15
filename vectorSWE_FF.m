function [q,j_idx,s_idx,m_idx,n_idx] = vectorSWE_FF(Eff_theta,Eff_phi,theta,phi,dtheta,dphi,J)
%VECTORSWE  Far-field vectorial spherical wave expansion
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
%   Input: 
%       Eff_theta   :   Far-field E-field Theta component.
%                       Dimension 1xT*P
%       Eff_phi     :   Far-field E-field Phi component.
%                       Dimension 1xT*P
%       theta   :   Theta vector where E-field shall be reconstructed. 
%                   Dimension 1xT*P where T is the number of Theta points 
%                   and P is the number of Phi points.
%       phi     :   Phi vector. Dimension 1xT*P.
%       dtheta  :   Theta sample spacing
%       dphi    :   Phi sample spacing
%       J       :   Maximum number of modes (compressed notation)
%                   J = 2*N*(N+2)
%
%
%   Output:
%       q       :   Spherical wave expansion coefficients. Dimensions 1xJ.
%                   J is the number of modes in compressed notation.
%       j_idx,
%       s_idx,
%       m_idx,  :   Mode indices. Dimension 1xJ
%       n_idx
%
%
%   Hannes Bartle
%   EPFL Microwaves and Antennas Group
%   2023

N = -1+sqrt(1+J/2);
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

%% Get SWE coefficients
%   Again for better readability, the SWE is implemented using a for loop
%   over all modes. This could be optimized for speed if a large number
%   of modes is expected.


q_tmp = repmat(Eff_theta,[1,length(j_idx)]).*conj(K_theta).*repmat(sin(theta),[1,length(j_idx)])...
    + repmat(Eff_phi,[1,length(j_idx)]).*conj(K_phi).*repmat(sin(theta),[1,length(j_idx)]);

for j = 1:length(j_idx)
    
    Q = reshape(q_tmp(:,j),[],length(unique(theta)));
    % Amplitude Scaling to (somewhat) match FEKO. This needs to be checked!
    A = 16*pi; 
    q(j) = 1/A*trapz(dtheta,trapz(dphi,Q));
end
end

