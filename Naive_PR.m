function [Wopt, IT, ffval]=Naive_PR(W0,prec_loss, prec_grad,maxit,gradtype,Target)
global N;
global M;
global P;
global wkspace;


W = W0;
it=1;
[eigen_values, eigen_vectors] = Band_Structure(W);
current_cost = cost(eigen_values,Target);
history_pot(it,1:2*P+1) = W';
history_cost(it,1)= current_cost;
if current_cost < prec_loss
    error('\n [OPTIMAL FOUND] : Your initial guess is already optimal \n');
end
dW = Gradient_real(W,Target,eigen_values,eigen_vectors,gradtype);
dW2= dW(1:P+1,1);
dir2 = -dW2;
tau2=Armijo_Line_Search(current_cost,dir2,dW, W, Target);
s2 = tau2*dir2;
s22(1:P+1,1)= s2';
for k=1:P+1
    s22(2*P+1-k+1,1)= conj(s2(k));
end
W = W + s22;
ss_prev2 = dir2;
dir_prev2= dir2;
nor_dW= norm(dW);
history_grad(it,1)= nor_dW;
while (it < maxit) && (nor_dW>prec_grad) %&&(current_cost>prec_loss)
    if (mod(it,20)==0)
        fprintf('\n[PR] : Iter %d, cost =%d, nor_dW= %d, step [%d]',it,current_cost, nor_dW,tau2);
    end
    it = it+1;
    [eigen_values, eigen_vectors] = Band_Structure(W);
    dW= Gradient_real(W,Target,eigen_values,eigen_vectors,gradtype);
    dW2= dW(1:P+1,1);
    nor_dW = norm(dW);
    dir_current2 = - dW2;
    bb2= (dir_current2'*(dir_current2-dir_prev2))/(dir_prev2'*dir_prev2);
    ss_current2 = dir_current2 + bb2*ss_prev2;
    tau2= Armijo_Line_Search(current_cost,ss_current2,dW,W,Target);
    if tau2 < eps 
      break; 
    end
    s2 = tau2*ss_current2;
    s22(1:P+1,1)= s2';
    for k=1:P+1
        s22(2*P+1-k+1,1)= conj(s2(k));
    end
    W = W + s22;
    ss_prev2 = ss_current2;
    dir_prev2 = dir_current2;   
    [eigen_values, eigen_vectors] = Band_Structure(W);
    current_cost= cost(eigen_values,Target);
    history_pot(it,1:2*P+1) = W';
    history_cost(it,1)= current_cost;
    history_grad(it,1)= nor_dW;
end
Wopt = W;
IT = it;
ffval = current_cost;
fprintf(1,'\n\n[POLAK RIBIERE] : %d iteration to converge with final cost %d and grad %d\n',it,current_cost, nor_dW);
History_pot_mat =  strcat(wkspace,'/PR/',wkspace,'_History_Pot.mat');
save(History_pot_mat, 'history_pot');
History_cost_mat =  strcat(wkspace,'/PR/',wkspace,'_History_Cost.mat');
save(History_cost_mat, 'history_cost');
History_grad_mat =  strcat(wkspace,'/PR/',wkspace,'_History_Grad.mat');
save(History_grad_mat, 'history_grad');
end