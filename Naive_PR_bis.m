function [Wopt]=Naive_PR_bis(W0,prec_loss, prec_grad,maxit,gradtype,Target)
global N;       
global M;     
global P; 
global wkspace; 




% %%%% With gradient 
    function [f,g] = Objective(W, Target)
             [bands_W, vects_W] = Band_Structure(W); 
             f=cost(bands_W,Target);
             g= Gradient_real(W,Target,bands_W,vects_W,gradtype);
    end
    options = optimoptions(@fminunc,'Display','iter-detailed','Algorithm','quasi-newton', 'HessUpdate', 'dfp',...
        'MaxIterations', maxit,'SpecifyObjectiveGradient',true, 'MaxFunctionEvaluations', 10^6,'StepTolerance',10^(-25),...
        'OptimalityTolerance', 10^-8);
    
    [Wopt, fval, exitflag] = fminunc(@(W) Objective(W,Target),W0,options)
    while (fval > prec_loss) 
    [Wopt, fval, exitflag2] = fminunc(@(W) Objective(W,Target),Wopt,options)
    end
end