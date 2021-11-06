function [Wopt, IT, ffval]=Naive_BFGS(W0,prec_loss, prec_grad,maxit,gradtype,Target)
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
[eigen_values, eigen_vectors] = Band_Structure(W);
current_cost = cost(eigen_values,Target);
dW = Gradient_real(W,Target,eigen_values,eigen_vectors,gradtype);
dW2= dW(1:P+1,1);
nor_dW = norm(dW);
history_grad(it,1)= nor_dW;
H2 = eye(P+1);
tau2=Inf;
while  (it < maxit) && (nor_dW>prec_grad)%&&(current_cost>prec_loss)
    if (mod(it,20)==0)
        fprintf('\n[BFGS] : Iter %d, cost =%d, nor_dW= %d, step [%f]',it,current_cost, nor_dW,tau2);
    end
    old_cost = current_cost;
    it = it+1;
    dir2 = - pinv(H2)*dW2;
    tau2 = Armijo_Line_Search(current_cost,dir2,dW2, W, Target);
    s2 = tau2*dir2;
    s22(1:P+1,1)= s2';
    for k=1:P+1
        s22(2*P+1-k+1,1)= conj(s2(k));
    end
    W = W + s22;
    [eigen_values, eigen_vectors] = Band_Structure(W);
    current_cost= cost(eigen_values,Target);
    dW_old = dW;
    dW = Gradient_real(W,Target,eigen_values,eigen_vectors,gradtype) ;
    dW2= dW(1:P+1,1);
    dW_old2 = dW_old(1:P+1,1);
    y2  = dW2-dW_old2;
    nor_dW = norm(dW);
    % update hessian with bfgs rule
    H2 = H2+ (y2*y2')/(y2'*s2) - ((H2*s2)*(H2*s2)')/(s2'*H2*s2);
    if (old_cost < current_cost)
        fprintf('\n ALARM  ! Cost is increasing ');
    end
    history_pot(it,1:2*P+1) = W';
    history_cost(it,1)= current_cost;
    history_grad(it,1)= nor_dW;
end

Wopt = W;
IT = it;
ffval = current_cost;
fprintf(1,'\n\n[BFGS] : %d iteration to converge with final cost %.10f and grad %.10f\n',it,current_cost, nor_dW);
History_pot_mat =  strcat(wkspace,'/BFGS/',wkspace,'_History_Pot.mat');
save(History_pot_mat, 'history_pot');
History_cost_mat =  strcat(wkspace,'/BFGS/',wkspace,'_History_Cost.mat');
save(History_cost_mat, 'history_cost');
History_grad_mat =  strcat(wkspace,'/BFGS/',wkspace,'_History_Grad.mat');
save(History_grad_mat, 'history_grad');
end