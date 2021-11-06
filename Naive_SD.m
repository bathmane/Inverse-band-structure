function [Wopt, IT, ffval]=Naive_SD(W0,prec_loss, prec_grad,maxit,gradtype,Target)
global P;
global N;
global M; 
global wkspace;



W = W0;
it = 1;
[eigen_values, eigen_vectors] = Band_Structure(W);
init_cost = cost(eigen_values,Target);
if init_cost < prec_loss
    error('\n [OPTIMAL FOUND] : Your initial guess is already optimal ! \n');
end
[eigen_values, eigen_vectors] = Band_Structure(W);
current_cost = cost(eigen_values,Target);
dW = Gradient_real(W,Target,eigen_values,eigen_vectors,gradtype);
nor_dW = norm(dW);
history_pot(it,1:2*P+1) = W';
history_cost(it,1)= current_cost;
history_grad(it,1)= nor_dW;
tau = Inf;
while (it < maxit)&&(nor_dW>prec_grad)%&&(current_cost>prec_loss)
     if (mod(it,20)==0)
        fprintf(1,'\n [SD] : Iteration %d: cost=%d, nor_dW=%d, step=[%d]',it,current_cost, nor_dW,tau);
     end
    it = it+1;
    [eigen_values, eigen_vectors] = Band_Structure(W);
    dW= Gradient_real(W,Target,eigen_values,eigen_vectors,gradtype);
    dW2= dW(1:P+1,1);
    tau= Armijo_Line_Search(current_cost,-dW2,dW,W,Target);
    if (tau < eps)
        break; 
    end 
    W = W-tau*dW;
    nor_dW = norm(dW);
    [eigen_values, eigen_vectors] = Band_Structure(W);
    current_cost = cost(eigen_values,Target);
    history_pot(it,1:2*P+1) = W';
    history_cost(it,1)= current_cost;
    history_grad(it,1)= nor_dW;
end


Wopt= W;
IT= it;
ffval= current_cost;
fprintf(1,'\n\n[Steepest Descent]: %d iteration to converge with final cost %d and grad %d\n',it,current_cost,nor_dW);
History_pot_mat =  strcat(wkspace,'/SD/',wkspace,'_History_Pot.mat');
save(History_pot_mat, 'history_pot');
History_cost_mat =  strcat(wkspace,'/SD/',wkspace,'_History_Cost.mat');
save(History_cost_mat, 'history_cost');
History_grad_mat =  strcat(wkspace,'/SD/',wkspace,'_History_Grad.mat');
save(History_grad_mat, 'history_grad');

end
