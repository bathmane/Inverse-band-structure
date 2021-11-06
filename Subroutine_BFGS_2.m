function [Wopt,fval]=Subroutine_BFGS_2(adapt_N,adapt_P,W0,prec_loss, prec_grad,maxit,gradtype,Target)
 P= adapt_P; 
 N = adapt_N; 


 
 
 
 
     function [f,g] = Objective(W,adapt_N,adapt_P,Target)
             [bands_W, vects_W] = Band_Structure_Progress(adapt_N,adapt_P,W); 
             f=cost(bands_W,Target);
             g= Gradient_Adaptive(adapt_P,adapt_N,W,Target,bands_W,vects_W,gradtype);
     end
    options = optimoptions(@fminunc,'Display','iter-detailed','Algorithm','quasi-newton', 'HessUpdate', 'bfgs',...
        'MaxIterations', maxit,'SpecifyObjectiveGradient',true, 'MaxFunctionEvaluations', 10^6,'StepTolerance',10^(-25),...
        'OptimalityTolerance', 10^-8);
    
    [Wopt, fval, exitflag] = fminunc(@(W) Objective(W,adapt_N,adapt_P,Target),W0,options); 
%     while (fval > prec_loss) 
%     [Wopt, fval, exitflag2] = fminunc(@(W) Objective(W,Target),Wopt,options)
%     end
 
end 
 
 
 
 
 
 
 
 
 
 
 
 
 
% %  
% %  
% % W = W0;
% % it=1;
% % [eigen_values, eigen_vectors] = Band_Structure_Progress(adapt_N,adapt_P,W);
% % current_cost = cost(eigen_values,Target);
% % dW = Gradient_Adaptive(adapt_P,adapt_N,W,Target,eigen_values,eigen_vectors,gradtype);
% % nor_dW = norm(dW);
% % H2 = eye(2*P+1);
% % while  (it < maxit) && (nor_dW>prec_grad)%&&(current_cost>prec_loss)
% %     %fprintf('\n[BFGS] : Iter %d, cost =%.5f, nor_dW= %.5f, step [%f]',it,current_cost, nor_dW,tau2);
% %     old_cost = current_cost;
% %     it = it+1;
% %     dir = - pinv(H2)*dW;
% %     tau2 = Armijo_Line_Search_Progress(adapt_N,adapt_P,current_cost,dir2,dW2, W, Target);
% %     s2 = tau2*dir2;
% %     s22(1:P+1,1)= s2';
% %     for k=1:P+1
% %         s22(2*P+1-k+1,1)= conj(s2(k));
% %     end
% %     W = W + s22; 
% %     [eigen_values, eigen_vectors] = Band_Structure_Progress(adapt_N,adapt_P,W);
% %     current_cost= cost(eigen_values,Target);
% %     dW_old = dW;
% %     dW = Gradient_Adaptive(adapt_P,adapt_N,W,Target,eigen_values,eigen_vectors,gradtype) ;
% %     dW2= dW(1:P+1,1);
% %     dW_old2 = dW_old(1:P+1,1);
% %     y2  = dW2-dW_old2;
% %     nor_dW = norm(dW);
% %     % update hessian with bfgs rule
% %     H2 = H2+ (y2*y2')/(y2'*s2) - ((H2*s2)*(H2*s2)')/(s2'*H2*s2);
% % end
% % 
% % Wopt = W;
% % ffval = current_cost;
% % end