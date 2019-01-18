function [var_eps,VAF_RW] = AOloopRW(G,H,Cphi0,sigmae,phisim)
% This MATLAB routine creates the Random Walk model
% and computes the variance of the residual wavefront var_eps, from system
% matrices G, H, covariance matrices Cphi0, sigmae*I, Cw, and the wavefront data phiSim

    % Compute sizes of given matrices
    [m2,N] = size(phisim);  
    o = size(G,1);
    p = size(H,1);
    % Allocate 
    eps = zeros(p,N);
    s = zeros(o,N);
    u = zeros(p,N);
    
    Gamma = Cphi0*G'*inv(G*Cphi0*G' + sigmae^2*eye(o));

    eps(:,1)     = phisim(:,1) -H*u(:,1); 
    s(:,1)       = G*eps(:,1) + sigmae^2*eye(o)*randn(o,1);
    eps_est(:,1) = Gamma*s(:,1);
    eps_est_no_mean(:,1) = eps_est(:,1) -mean(eps_est(:,1));

    delta_u(:,1) = pinv(H)*Gamma*s(:,1); 
    u(:,1)       = delta_u(:,1) + u(:,1);
    phi_est(:,1) = eps_est(:,1);

    for k = 1:N-1
        phi_t(:,k+1)    = phisim(:,k);
        eps(:,k+1)      = phisim(:,k+1) -H*u(:,k);
        s(:,k+1)        = G*eps(:,k+1) +  sigmae^2*eye(o)*randn(o,1);
        eps_est(:,k+1)  = Gamma*s(:,k+1);
        phi_est(:,k+1)  = eps_est(:,k+1) + H*u(:,k);
        eps_est_no_mean(:,k+1)  = eps_est(:,k+1) -mean(eps_est(:,k+1));
        
        
        delta_u(:,k+1)  = pinv(H)*Gamma*s(:,k+1);
        u(:,k+1)        = delta_u(:,k+1) + u(:,k);
        
        eps_est_1(:,k+1) = eps_est(:,k+1) -H*delta_u(:,k+1);
        phi_est(:,k+1)   = eps_est_1(:,k+1)+H*u(:,k);
    end 
    VAF_RW = vaf(phisim,phi_est);
    var_eps = var((eps-eps_est_1),0,2);
end
