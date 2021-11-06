% Computes the Gradient of J in the point W 
function grad = Gradient_real(W,Target,values,vectors,Method)
global Brillouin;  
global Gamma;   
global Q;        
global N;         
global P; 
global L;        
global M;       


switch Method 
    case 'FINITE_DIFFERENCES'        
        grad = zeros(2*P+1,1);
        h=10^(-10); 
        [get_bandes, get_eigen_vectors] = Band_Structure(W); 
        co = cost(get_bandes,Target);      
        
        for k=1:P
            Wh = W; 
            Wh(k,1)=W(k,1)+h; 
            Wh(2*P+1-k+1,1)=conj(Wh(k,1));           
            %fprintf('\n%.10f\n%.10f\n%.10f\n', Wh(1,1), Wh(2,1), Wh(3,1)); 
            [get_bandes_h, get_eigen_vectors] = Band_Structure(Wh); 
            c_h= cost(get_bandes_h,Target); 
            grad(k,1) = (c_h - co)/h;
            grad(2*P+1-k+1,1)= conj(grad(k,1));          
        end       
        
            Wp = W; 
            Wp(P+1,1)=W(P+1,1)+h; 
            %fprintf('\n%.10f\n%.10f\n%.10f\n', Wp(1,1), Wp(2,1), Wp(3,1)); 
            [get_bandes_h, get_eigen_vectors] = Band_Structure(Wp); 
            c_h= cost(get_bandes_h,Target); 
            grad(P+1,1) = (c_h - co)/h;

        
end 
end