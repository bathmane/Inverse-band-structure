
% Computes the cost between the computed Bandes and the Target energies 
function cosst = cost(Bandes,Target)
global Brillouin
global Gamma;   
global Q;       
global N;         
global P; 
global L;        
global M;       


temp = 0; 
for m =1:M
    for q=1:Q
            temp = temp +(Bandes(m,q) - Target(m,q))*(Bandes(m,q) - Target(m,q));    
    end 
end
cosst = (1/Q)*temp; 

% Ok 


end