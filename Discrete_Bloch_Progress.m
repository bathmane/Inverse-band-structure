
% Build the rigidity matrix K_{ij} = a(e_j, e_i)
function [K] = adaptive_rigidity_matrix(adapt_N,adapt_P,q,V)
% The full version / See the SPARSE version below 
K = zeros(2*adapt_N+1,2*adapt_N+1); 
Vhat0 = V(adapt_P+1,1);   % Need the Fourier coef for k=0 
for k=-adapt_N:adapt_N
    for l=-adapt_N:adapt_N
        kp= k+adapt_N+1;
        lp=  l+adapt_N+1; 
        if (k==l)
           K(kp,lp) = (k+q)^2 + Vhat0; 
        else
           if (k-l<=adapt_N)&&(k-l>=-adapt_N)
               ind = k-l +adapt_P+1; 
               if(ind <= 2*adapt_P+1)&&(ind>=1)
                  K(kp,lp) = V(ind,1);
               else 
                   K(kp,lp) =0;
               end                    
           else
               K(kp,lp) = 0;
           end
        end
    end
end





% The Sparse version  
% K = sparse(2*N+1,2*N+1); 
% Vhat0 = V(1,1);  % Need the Fourier coef for k=0 
% Np= 2*N+1;
% Nlin = 1:1:Np; 
% diag= (Nlin+Nlin).^2 + Vhat0; 
% other = ??? ; 
% D = sparse(1:Np,1:Np,diag*ones(1,Np),Np,Np);
% E = sparse(2:Np,1:Np-1,ones(1,Np-1),Np,Np);
% K = E+D+E' ; 










end 