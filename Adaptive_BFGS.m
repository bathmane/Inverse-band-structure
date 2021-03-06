function [Wopt, IT, ffval,Progr]=Adaptive_BFGS(P0,N0,W0,prec_loss, prec_grad,eta,maxit,gradtype,Target)
global Q;
global M;
global wkspace;
global full_N;

% Initialization 
adapt_P =P0;
adapt_N= N0;
W = W0;
it = 0;
tau = Inf;
stop_criterion = true;

% Pre-Descent 
[eigen_values, eigen_vectors] = Band_Structure_Progress(adapt_N,adapt_P,W);
init_cost = cost(eigen_values,Target);
current_cost = init_cost;
dW = Gradient_Adaptive(adapt_P,adapt_N,W,Target,eigen_values,eigen_vectors,gradtype);
dW2= dW(1:adapt_P+1,1);
dir2 = -dW2;
tau2=Armijo_Line_Search_Progress(adapt_N, adapt_P,current_cost,dir2,dW, W, Target);
s2 = tau2*dir2;
s22(1:adapt_P+1,1)= s2';
for k=1:adapt_P+1
    s22(2*adapt_P+1-k+1,1)= conj(s2(k));
end
W = W + s22;
nor_dW= sqrt(sum(dW.^2));

% Optim Loop
it = 1;
history_pot{it} = W;
history_cost(it,1)= current_cost;
history_grad(it,1) = nor_dW;
Progr(it,1)= adapt_N;
Progr(it,2)= adapt_P;
while (stop_criterion)    
    it = it+1;
    history_pot{it} = W;
    history_cost(it,1)= min(current_cost, history_cost(it-1,1));
    history_grad(it,1) = nor_dW;
    Progr(it,1)= adapt_N;
    Progr(it,2)= adapt_P;
    
    fprintf(1,'\n---[Adaptive BFGS] : Iter %d: cost=%d, nor_dW=%d, tau=%d, P=%d, N=%d',it,current_cost,nor_dW,tau,adapt_P,adapt_N);
    % Perform Inner Descent 
    [W,ffval,nor_dW]=Subroutine_BFGS(adapt_N,adapt_P,W,prec_loss,prec_grad,maxit,gradtype,Target);
    [bands_W, vectors_W] = Band_Structure_Progress(adapt_N,adapt_P,W);
    current_cost = cost(bands_W, Target);    
    
    % Error contribution E1
    Post_Err_2 = Aposteriori_Error_2(bands_W,vectors_W,W,adapt_P,adapt_N,full_N);
    E1 = 0;
    for m=1:M
        for q=1:Q
            %D_mq  = abs(exact_bands(m,q) - bands_W(m,q));
            D_mq = Post_Err_2(m,q);
            E1 = E1 + (2*abs( Target(m,q) - bands_W(m,q)) + D_mq )*D_mq ;
        end
    end
    E1 =E1/Q;
    if(E1 > eta)
        adapt_N = adapt_N +1;
    end
    
    
    % Error contribution 2
    P2 = 2*adapt_P;
    W2P  = zeros(2*P2+1,1);
    W2P(P2+1,1) = W(adapt_P+1, 1);
    for ind = 1:adapt_P
        W2P(P2+1+ind,1) = W(adapt_P+1+ind,1);
        W2P(P2+1-ind,1) = W(adapt_P+1-ind,1);
    end
    [P2_bands,P2_vectors] = Band_Structure_Progress(adapt_N,P2,W2P);
    dW2P= Gradient_Adaptive(P2,adapt_N,W2P,Target, P2_bands, P2_vectors,gradtype);
    E_gr  = norm(dW2P);
    [C,I] = max( dW2P(  P2+1+adapt_P+1: 2*P2+1));
    Ptild = I + adapt_P;
    Wtild = zeros(2*Ptild+1, 1);
    Wtild(Ptild+1,1) = W(adapt_P+1, 1);
    for ind = 1:adapt_P
        Wtild(Ptild+1+ind,1) = W(adapt_P+1+ind,1);
        Wtild(Ptild+1-ind,1) = W(adapt_P+1-ind,1);
    end
    [Ptild_bands,Ptild_vectors] = Band_Structure_Progress(adapt_N,Ptild,Wtild);
    dWtild= Gradient_Adaptive(Ptild,adapt_N,Wtild,Target, Ptild_bands, Ptild_vectors,gradtype);
    Diff = dWtild;
    Diff(Ptild+1,1) =0;
    for ind = 1:adapt_P
        Diff(Ptild+1+ind,1)=0;
        Diff(Ptild+1-ind,1)=0;
    end
    E2 = norm(Diff);
    if (E2 > eta)
        adapt_P=Ptild;
        W = Wtild;
    end
    
    % Stopping Criterion 
    stop_criterion = (E1 > eta)||(E2 > eta)||(nor_dW > prec_grad);   
end
% Save Outputs 
Wopt = W;
IT = it;
ffval = current_cost;
fprintf(1,'\n\n[Adaptive BFGS] : %d iteration to converge with final cost %f\n',it,current_cost);
History_pot_mat =  strcat(wkspace,'/BFGS/',wkspace,'_History_Pot_progress.mat');
save(History_pot_mat, 'history_pot');
History_cost_mat =  strcat(wkspace,'/BFGS/',wkspace,'_History_Cost_progress.mat');
save(History_cost_mat, 'history_cost');
History_grad_mat =  strcat(wkspace,'/BFGS/',wkspace,'_History_Grad_progress.mat');
save(History_grad_mat, 'history_grad');
end