 function [var_eps,VAF_SL] = AOloopRWSlopes(G,H,Cphi0,sigmae,phisim)
% This MATLAB routine creates the Random Walk model adjusted for optimizing
% slopes and computes the variance of the residual wavefront var_eps, from system
% matrices G, H, covariance matrices Cphi0, sigmae*I, Cw, and the wavefront data phiSim
    % Determine sizes of matrices
    [m2,N] = size(phisim); 
    o = size(G,1);
    p = size(H,1);
    % Allocate matrixes
    eps = zeros(p,N);
    sSL = zeros(o,N);
    u = zeros(p,N);
    %Define shorthand
    Gamma = Cphi0*G'*inv(G*Cphi0*G' + sigmae^2*eye(o));
    % Set first values
    eps(:,1)     = phisim(:,1) -H*u(:,1); 
    s(:,1)       = G*eps(:,1) + sigmae^2*eye(o)*randn(o,1);
    eps_est(:,1) = Gamma*s(:,1);

    delta_u(:,1) = pinv(G*H)*G*Gamma*s(:,1);
    u(:,1)       = delta_u(:,1) + u(:,1);
    phi_est(:,1) = eps_est(:,1);
    
    % Create the loop based on the models in the report
    for k = 1:N-1
        phi_t(:,k+1)    = phisim(:,k);
        eps(:,k+1)      = phisim(:,k+1) -H*u(:,k);
        s(:,k+1)        = G*eps(:,k+1) +  sigmae^2*eye(o)*randn(o,1);
        eps_est(:,k+1)  = Gamma*s(:,k+1);
        phi_est(:,k+1)  = eps_est(:,k+1) + H*u(:,k);
        
        delta_uSL(:,k+1)  = pinv(G*H)*G*Gamma*s(:,k+1);
        u(:,k+1)          = delta_uSL(:,k+1) + u(:,k);
        
        eps_est_1(:,k+1) = eps_est(:,k+1) -H*delta_uSL(:,k+1);
        phi_est(:,k+1)   = eps_est_1(:,k+1)+H*u(:,k);
        eps_est_1_no_mean(:,k+1)  = eps_est_1(:,k+1) -mean(eps_est_1(:,k+1));
    end 
    VAF_SL= vaf(phisim,phi_est);
    var_eps = mean(var((eps-eps_est_1_no_mean)));
end
