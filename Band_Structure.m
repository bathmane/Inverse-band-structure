% Computes the 2N+1 first energy bands with Vhat 
function [All_bands, All_eigen_vectors] = Band_Structure(Vhat)
global Brillouin;   
global Q;         
global N;                 

% Remember that we use only N+1 Fourier Coefs. Hence 2N+1 bands 
All_bands=zeros(2*N+1,Q);
All_eigen_vectors= zeros(2*N+1,2*N+1,Q); 



for q =1:Q
    q_point     = Brillouin(q); 
    Kmat        = Discrete_Bloch(q_point, Vhat);
    [Vect,D]    = eig(Kmat);  
    All_bands(:,q) = diag(D); 
    All_eigen_vectors(:,:,q) =Vect;
end
end 
