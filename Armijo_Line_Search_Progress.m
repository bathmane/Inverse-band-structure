function res = Armijo_Line_Search_Progress(adapt_N, adapt_P,cost_curr,dir,grad_curr, W_curr, Target)
   
    ainit = 1000;    
    sig = 0.1; 
    max_iter= 10000;  
    loop_a=ainit; 
    it =0; 
    dir2 = zeros(2*adapt_P+1,1); 
    dir2(1:adapt_P+1,1)= dir'; 
    for k=1:adapt_P
        dir2(2*adapt_P+1-k+1,1)= conj(dir(k)); 
    end
    dir = dir2; 
    %% Computing the cost at the point xk = Beta+loop_a*s_curr
    Xk = W_curr + loop_a*dir; 
    [eigen_values, eigen_vectors] = Band_Structure_Progress(adapt_N,adapt_P,Xk); 
    loopf = cost(eigen_values,Target); 
    
    %while (loopf > cost_curr + lam*loop_a*grad_curr*s_curr)&&(it<max_iter) % Relaxed
    while (loopf > cost_curr)&&(it <max_iter)
        loop_a = sig*loop_a;    
        Xk = W_curr + loop_a*dir; 
        [eigen_values, eigen_vectors] = Band_Structure_Progress(adapt_N,adapt_P,Xk); 
        loopf = cost(eigen_values,Target);
        it = it+1; 
    end 
    res = loop_a; 
end 
