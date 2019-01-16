function [eps_estSL,var_epsSL,delta_uSL,sSL] = AOloopRWSlopes(G,H,Cphi0,sigmae,phisim)
    [m2,N] = size(phisim);  %Change phisim to phiSim!!
    o = size(G,1);
    p = size(H,1);
    
    eps = zeros(p,N);
    sSL = zeros(o,N);
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
    sSL(:,1)       = G*eps(:,1) + sigmae^2*eye(o)*randn(o,1);
    eps_estSL(:,1) = Gamma*sSL(:,1);
    eps_estSL(:,1) = eps_estSL(:,1) -mean(eps_estSL(:,1));

    delta_uSL(:,1) = (pinv(G*H))*G*Gamma*sSL(:,1); % % %
    u(:,1)       = delta_uSL(:,1) + u(:,1);

    for k = 1:N-1
        eps(:,k+1)      = phisim(:,k+1) -H*u(:,k);
        sSL(:,k+1)        = G*eps(:,k+1) +  sigmae^2*eye(o)*randn(o,1);
        eps_estSL(:,k+1)  = Gamma*sSL(:,k+1);
        eps_estSL(:,k+1)  = eps_estSL(:,k+1) -mean(eps_estSL(:,k+1));

        delta_uSL(:,k+1)  =(pinv(G*H))*G*Gamma*sSL(:,k+1);
        u(:,k+1)        = delta_uSL(:,k+1) + u(:,k);
        
        %res_slopes(:,k+1) = s(:,k+1) - G*H*delta_u(:,k+1);   
    end
    sSL =sSL ;
  eps_estSL = eps_estSL; 
  var_epsSL = var(eps);
  %res_slopes = res_slopes;
end
