% Computes the 2N+1 first energy bands with Vhat 
function [All_bands, All_eigen_vectors] = Band_Structure_Progress(adapt_N,adapt_P,Vhat)
global Brillouin; 
global Q;               

All_bands=zeros(2*adapt_N+1,Q);
All_eigen_vectors= zeros(2*adapt_N+1,2*adapt_N+1,Q); 


for q =1:Q
    q_point = Brillouin(q); 
    Kmat = Discrete_Bloch_Progress(adapt_N,adapt_P,q_point, Vhat);
    [Vect,D] = eig(Kmat); 
    All_bands(:,q) = diag(D); 
    All_eigen_vectors(:,:,q) =Vect;
    
   
end
end 
