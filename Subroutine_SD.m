function [Wopt,ffval,nor_dW]=Subroutine_SD(adapt_N,adapt_P,W0,prec_loss, prec_grad,maxit,gradtype,Target)
P=adapt_P; 
N= adapt_N; 
W = W0;
it = 1;
[eigen_values, eigen_vectors] = Band_Structure_Progress(adapt_N,adapt_P,W);
[eigen_values, eigen_vectors] = Band_Structure_Progress(adapt_N,adapt_P,W);
current_cost = cost(eigen_values,Target);
dW = Gradient_Adaptive(adapt_P,adapt_N,W,Target,eigen_values,eigen_vectors,gradtype);
nor_dW = norm(dW);
while (it < maxit)&&(nor_dW>prec_grad)%&&(current_cost>prec_loss)
    it = it+1;     
    [eigen_values, eigen_vectors] = Band_Structure_Progress(adapt_N,adapt_P,W);
    dW= Gradient_Adaptive(adapt_P,adapt_N,W,Target,eigen_values,eigen_vectors,gradtype);
    dW2= dW(1:P+1,1);
    tau= Armijo_Line_Search_Progress(adapt_N,adapt_P,current_cost,-dW2,dW,W,Target) ;
    W = W-tau*dW;
    nor_dW = norm(dW);
    [eigen_values, eigen_vectors] = Band_Structure_Progress(adapt_N,adapt_P,W);
    current_cost = cost(eigen_values,Target);
end
Wopt= W;
ffval= current_cost;
end