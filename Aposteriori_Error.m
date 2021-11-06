function [post_error]=Aposteriori_Error(reduced_bands,reduced_eigen_vectors,Vpot,size_Pot,reduced_N, full_N)
global Brillouin;  
global Q;
global M;



post_error = zeros(M,Q);
for kpoint =1:Q
    full_Bloch = Discrete_Bloch_Progress(full_N, size_Pot, Brillouin(kpoint), Vpot);
    
%         if(~ishermitian(full_Bloch))
%             error('full_Bloch is not HERMITIAN');
%         end
%         if(~all(eig(full_Bloch))>0)
%             error(' full_Bloch is not positive');
%        end
    for band =1:M
        reduced_value=reduced_bands(band,kpoint);
        reduced_vector=reduced_eigen_vectors(:,band,kpoint);
        
        kappa =0.5;
        theta = 0.1; 
        
        lam_N = reduced_value;
        delta = theta*(lam_N-kappa); 
        
        
        
        b = lam_N+delta;
        bmat = b*eye(2*full_N+1,2*full_N+1);
        
        %%%% Prior bounds
        lam_tild = 0;
        a = lam_tild+delta;
        amat = a*eye(2*full_N+1,2*full_N+1);
        
        %%%% projection of eigenvector
        UN = zeros(2*full_N+1,1);
        UN(full_N+1-reduced_N:full_N+1+reduced_N,1) = reduced_vector;
        
        %%%% solve the linear system
        Residual = full_Bloch*UN-lam_N*UN;
        error_val = Residual'*((inv(full_Bloch-bmat))*(full_Bloch-amat)*(inv(full_Bloch-bmat)))*Residual;
        post_error(band,kpoint) =error_val ;
    end
end

end 