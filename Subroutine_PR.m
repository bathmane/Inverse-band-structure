function [Wopt,ffval,nor_dW]=Subroutine_PR(adapt_N,adapt_P,W0,prec_loss, prec_grad,maxit,gradtype,Target)
W = W0;
it=1;
P = adapt_P; 
N = adapt_N;
[eigen_values, eigen_vectors] = Band_Structure_Progress(adapt_N,adapt_P,W);
current_cost = cost(eigen_values,Target);
dW = Gradient_Adaptive(adapt_P,adapt_N,W,Target,eigen_values,eigen_vectors,gradtype);
dW2= dW(1:P+1,1);
dir2 = -dW2;
tau2=Armijo_Line_Search_Progress(adapt_N,adapt_P,current_cost,dir2,dW, W,Target);
s2 = tau2*dir2;
s22(1:P+1,1)= s2';
for k=1:P+1
    s22(2*P+1-k+1,1)= conj(s2(k));
end
W = W + s22;
ss_prev2 = dir2;
dir_prev2= dir2;
nor_dW= norm(dW);
while (nor_dW>prec_grad)&&(current_cost>prec_loss)&&(it < maxit) 
    if (mod(it,20)==0)
        fprintf('\n\t[sub PR] : Iter %d, cost =%d, nor_dW= %d, step [%f]',it,current_cost, nor_dW,tau2);
    end   
    old_cost = current_cost;
    it = it+1;
    [eigen_values, eigen_vectors] = Band_Structure_Progress(adapt_N,adapt_P,W);
    dW= Gradient_Adaptive(adapt_P,adapt_N,W,Target,eigen_values,eigen_vectors,gradtype);
    dW2= dW(1:P+1,1);
    nor_dW = norm(dW);
    dir_current2 = - dW2;
    bb2= (dir_current2'*(dir_current2-dir_prev2))/(dir_prev2'*dir_prev2);
    ss_current2 = dir_current2 + bb2*ss_prev2;
    tau2= Armijo_Line_Search_Progress(adapt_N,adapt_P,current_cost,ss_current2,dW,W,Target);
    s2 = tau2*ss_current2;
    s22(1:P+1,1)= s2';
    for k=1:P+1
        s22(2*P+1-k+1,1)= conj(s2(k));
    end
    W = W + s22;
    ss_prev2 = ss_current2;
    dir_prev2 = dir_current2;
    
    [eigen_values, eigen_vectors] = Band_Structure_Progress(adapt_N,adapt_P,W);
    current_cost= cost(eigen_values,Target);
    if (old_cost < current_cost)
        fprintf('\n ALARM  ! Cost is increasing ');
    end
end
Wopt = W;
IT = it;
ffval = current_cost;
end