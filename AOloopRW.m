function [eps_est,delta_u,s,VAF] = AOloopRW(G,H,Cphi0,sigmae,phisim)
    [m2,N] = size(phisim);  %Change phisim to phiSim!!
    o = size(G,1);
    p = size(H,1);
    
    eps = zeros(p,N);
    s = zeros(o,N);
    u = zeros(p,N);
    
    Gamma = Cphi0*G'*inv(G*Cphi0*G' + sigmae^2*eye(o));

% [U,S,Vt] = svd(H);
% Lf = length(nonzeros(round(diag(S),10))); % Find the number of nonzero elements of the singular values of H rounded up to 10 decimal points
% if Lf == size(H,1) % F is  full row-rank
%     FR = 1;
% else % F is not full row-rank
%     FR = 0;
%     U1 = U(:,1:Lf);
%     S1 = S(1:L,1:Lf);
%     V1 = Vt(1:Lf,:)';
% end

    eps(:,1)     = phisim(:,1) -H*u(:,1); % % %
    s(:,1)       = G*eps(:,1) + sigmae^2*eye(o)*randn(o,1);
    eps_est(:,1) = Gamma*s(:,1);
    eps_est(:,1) = eps_est(:,1) -mean(eps_est(:,1));

    delta_u(:,1) = pinv(H)*Gamma*s(:,1); % % %
    u(:,1)       = delta_u(:,1) + u(:,1);

    for k = 1:N-1
        eps(:,k+1)      = phisim(:,k+1) -H*u(:,k);
        s(:,k+1)        = G*eps(:,k+1) +  sigmae^2*eye(o)*randn(o,1);
        eps_est(:,k+1)  = Gamma*s(:,k+1);
        eps_est(:,k+1)  = eps_est(:,k+1) -mean(eps_est(:,k+1));

        delta_u(:,k+1)  = pinv(H)*Gamma*s(:,k+1);
        u(:,k+1)        = delta_u(:,k+1) + u(:,k);
    end
    VAF = vaf(eps,eps_est(:,1:end))
end
