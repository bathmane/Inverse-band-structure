% Computes the Gradient of J in the point W 
function grad = Gradient_Adaptive(adapt_P,adapt_N,W,Target,values,vectors,Method)
   
switch Method 
    case 'FINITE_DIFFERENCES'        
        grad = zeros(2*adapt_P+1,1);
        h=10^(-8); 
        [get_bandes, get_eigen_vectors] = Band_Structure_Progress(adapt_N,adapt_P,W); 
        co = cost(get_bandes,Target);      
        
        for k=1:adapt_P
            Wh = W; 
            Wh(k,1)=W(k,1)+h; 
            Wh(2*adapt_P+1-k+1,1)=conj(Wh(k,1));           
            %fprintf('\n%.10f\n%.10f\n%.10f\n', Wh(1,1), Wh(2,1), Wh(3,1)); 
            [get_bandes_h, get_eigen_vectors] = Band_Structure_Progress(adapt_N,adapt_P,Wh); 
            c_h= cost(get_bandes_h,Target); 
            grad(k,1) = (c_h - co)/h;
            grad(2*adapt_P+1-k+1,1)= conj(grad(k,1));          
        end       
        
            Wp = W; 
            Wp(adapt_P+1,1)=W(adapt_P+1,1)+h; 
            %fprintf('\n%.10f\n%.10f\n%.10f\n', Wp(1,1), Wp(2,1), Wp(3,1)); 
            [get_bandes_h, get_eigen_vectors] = Band_Structure_Progress(adapt_N,adapt_P,Wp); 
            c_h= cost(get_bandes_h,Target); 
            grad(adapt_P+1,1) = (c_h - co)/h;

        
end 
end