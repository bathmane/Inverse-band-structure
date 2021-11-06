% Build the discere Bloch operator K_{ij} = a(e_j, e_i)
function [K] = Discrete_Bloch(q_point,V)
global N;
global P;

K = zeros(2*N+1,2*N+1);
Vhat0 = V(P+1,1);  % Need the Fourier coef for k=0
for k=-N:N
    for l=-N:N
        kp= k+N+1;
        lp=  l+N+1;
        if (k==l)
            K(kp,lp) = (k+q_point)^2 + Vhat0;
        else
            if (k-l<=N)&&(k-l>=-N)
                ind = k-l +P+1;
                if(ind <= 2*P+1)&&(ind>=1)
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
end 