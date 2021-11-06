clear all; 
clc; 
close all; 
fil =1;
fprintf(fil,'\n\t **********************************************************************');
fprintf(fil,'\n\t **                  Code:     Inverse Hill                        ****');
fprintf(fil,'\n\t **                Author:     Athmane Bakhta                      ****');
fprintf(fil,'\n\t **               Contact:     https://cermics.enpc.fr/~bakhtaa/   ****');
fprintf(fil,'\n\t **           Last update:     Dec  2017                           ****');
fprintf(fil,'\n\t **                                                                ****');
fprintf(fil,'\n\t **********************************************************************\n\n');


%%  global variables 
global Brillouin;   % Discrete 1st Brillouin zone
global Gamma;       % Gamma: unit cell of the lattice 
global Q;           % Number of k-points in Brillouin zone 
global L;           % Number of x points in Gamma
global M;           % Number of target bands in the cost function
global N;           % Numbe of basis functions 
global P;           % Number of Fourier modes in the potential
global full_N;      % Used in the progressive algorithm 
global wkspace; 
global fig;
global texFile; 
global out_files; 
fig=1; 

%% Initializations
N=40;                   
P=1;                    
M=1;                    
Q=20; 
L=60; 
full_N=1000; 

%% Create workspaces 
Nb_mod ='Inverse_'; 
dirp = 25; 
wkspace=sprintf('%s_P%d',Nb_mod,dirp);
mkdir(wkspace); 
wkspace2=sprintf('%s/BFGS', wkspace);
mkdir(wkspace2); 
wkspace2=sprintf('%s/PR', wkspace);
mkdir(wkspace2); 
wkspace2=sprintf('%s/SD', wkspace);
mkdir(wkspace2); 
texFile = fopen(sprintf('%s/Report_Run.txt', wkspace),'a') ;
out_files = [1, texFile]; 


%% Discretization of the reduced first Brillouin zone 
left=0;  
right=0.5;
Brillouin = left:(right-left)/(Q-1):right; 
  
%% Discretization of unit cell 
left=-pi;  
right=pi;
Gamma = left:(right-left)/(L-1):right ; 





Target = zeros(M,Q);
%%%-----------------------------------
%%%         Direct problem 
%%%----------------------------------- 
% %% The periodic potential Vper(x) 
Vper=zeros(2*P+1,1);
for k=1:P
    Vper(k,1)= rand(); 
    Vper(2*P+1-k+1,1)= conj(Vper(k,1));
end
Vbuild_Vper = zeros(1,L);
for x=1:L
    xc=Gamma(x);
    for k=-P:P
        kp = k+P+1;
        Vbuild_Vper(1,x) =Vbuild_Vper(1,x)+Vper(kp,1)*(exp(-1i*k*xc))/sqrt(2*pi);
    end
end
% 
% %% Shift the spectrum
% [get_bands_V, get_eigen_vectors_V] = Band_Structure(Vper); 
% min_eigen= abs(min(get_bands_V(1,:)))+1; 
% Vper(P+1) = min_eigen; 
% Vbuild_Vper = Vbuild_Vper+Vper(P+1); 
% 
% 
% %% Plot the potential 
% fig=fig+1;
% FIG=figure(fig);
% set(FIG, 'name', 'Periodic Potential Vper');
% hold on ;  
% plot(Gamma,Vbuild_Vper,'k', 'linewidth', 1.1);
% xlabel(gca,'x \in \Gamma ','FontSize',14);
% ylabel(gca,'V_{per}','FontSize',14);
% path_jpg =  sprintf('%s/Pot_P%d.jpg',wkspace,P); 
% saveas(FIG,path_jpg); 
% path_fig =  sprintf('%s/Pot_P%d.fig',wkspace,P); 
% saveas(FIG,path_fig);   
% 
% 
[bands_V, eigen_vectors_V] = Band_Structure(Vper); 
fig=fig+1;
FIG = figure(fig);
set (FIG , 'name', 'Target Energy bands' ) ;
hold on ; 
for r=1:M
   subplot(M,1,M-r+1);
   plot(Brillouin,bands_V(r,:), 'k', 'linewidth', 1);   
   Str = sprintf('band %d', r); 
   ylabel(Str, 'FontSize',11); 
   if (r==1)
       xlabel(' q \in  \Gamma^{*}', 'FontSize',11);   
   end 
   hold on ; 
end
Target = bands_V; 

P = dirp; 

%%%--------------------------------
%%% Create piecewise target 
%%%--------------------------------
% % b1 = bands_V(1,1); 
% % a1 = (bands_V(1,floor(Q/3)) - b1)/Brillouin(floor(Q/3));
% % seg_1 =  a1.*Brillouin(1:floor(Q/3)) + b1; 
% % Target(1,1:floor(Q/3))  = seg_1; 
% % a2 = 0.05 ;
% % b2 =  bands_V(1,floor(Q/3)) - a2 * Brillouin(floor(Q/3)); 
% % seg_2 =  a2.*Brillouin(floor(Q/3) :Q) + b2;
% % Target(1,floor(Q/3):Q)  = seg_2; 


Target = 0.*Target; 
Target(1,:)= 1; 




%%%-----------------------------------
%%%         Invserse  problem 
%%%-----------------------------------

% Initial guess 
W0 = zeros(2*P+1,1);
Vbuild_W0 = zeros(1,L); 
for x=1:L
    xc=Gamma(x); 
    for k=-P:P
        kp = k+P+1; 
        Vbuild_W0(1,x) =Vbuild_W0(1,x)+W0(kp,1)*(exp(-1i*k*xc))/sqrt(2*pi);     
    end
end


% optimization parameters 
gradtype = 'EXACT'; 
gradtype = 'FINITE_DIFFERENCES';
prec_loss = 10^-5; 
prec_grad = 10^-5; 
maxit     = 5*10^3;


% Test gradient 
% [values, vectors]= Band_Structure(W0); 
% cost(Target, values)
% grad = Gradient_real(W0,Target,values,vectors,gradtype)
% 


%% BFGS
tic; 
[BFGS_Wopt, BFGS_IT, f_BFGS]=Naive_BFGS(W0,prec_loss,prec_grad,maxit,gradtype,Target);
BFGS_time = toc; 


%% Polak Ribiere 

tic; 
[PR_Wopt, PR_IT, f_PR]=Naive_PR(W0,prec_loss, prec_grad,maxit,gradtype,Target);
% PR_Wopt=BFGS_Wopt;
% PR_IT=BFGS_IT;
% f_PR=f_BFGS;
PR_time = toc; 


%% Steepest Descent
tic; 
[SD_Wopt, SD_IT, f_SD] = Naive_SD(W0,prec_loss, prec_grad,maxit,gradtype,Target); 
% SD_Wopt=BFGS_Wopt;
% SD_IT=BFGS_IT;
% f_SD=f_BFGS;
SD_time = toc; 



%% Plot results of Naive search
[opt_bands_BFGS,get_eigen_vectors_V]=Band_Structure(BFGS_Wopt); 
[opt_bands_PR,get_eigen_vectors_V]=Band_Structure(PR_Wopt); 
[opt_bands_SD,get_eigen_vectors_V]=Band_Structure(SD_Wopt); 

Vbuild_bfgs = zeros(1,L); 
Vbuild_pr = zeros(1,L);
Vbuild_sd = zeros(1,L); 
Vbuild_orig = zeros(1,L); 
for x=1:L
    xc=Gamma(x); 
    for k=-P:P
         kp = k+P+1; 
%         Vbuild_orig(1,x) = Vbuild_orig(1,x)+Vper(kp,1)*(exp(-1i*k*xc))/sqrt(2*pi);
         Vbuild_bfgs(1,x) =Vbuild_bfgs(1,x)+BFGS_Wopt(kp,1)*(exp(-1i*k*xc))/sqrt(2*pi);
         Vbuild_pr(1,x) =Vbuild_pr(1,x)+PR_Wopt(kp,1)*(exp(-1i*k*xc))/sqrt(2*pi);  
         Vbuild_sd(1,x) =Vbuild_sd(1,x)+SD_Wopt(kp,1)*(exp(-1i*k*xc))/sqrt(2*pi);    
    end
end 




fig_pot=100;
FIG=figure(fig_pot);
hold on ;  
plot(Gamma, Vbuild_orig, '-kp',...
        Gamma,Vbuild_W0,'--k',...
        Gamma,Vbuild_bfgs,'-bs',...
        Gamma,Vbuild_pr, '-ro',...    
        Gamma,Vbuild_sd,'-m^');
ylabel(gca,'Potentials','FontSize',14);
xlim([-3.15 3.15]); 

fig_bands=200;
FIG = figure(fig_bands);
set (FIG , 'name', 'Target-optimal Energy bands' ) ; 
for r=1:M
   subplot(M,1,M-r+1);
   plot(Brillouin,Target(r,:),'-k',...
        Brillouin,opt_bands_BFGS(r,:),'-bs',...
        Brillouin,opt_bands_PR(r,:),'-ro',...
        Brillouin,opt_bands_SD(r,:),'-m^');
   Str = sprintf('band %d', r); 
   ylabel(Str, 'FontSize',11); 
   if (r==1)
       xlabel(' q \in  \Gamma^{*}', 'FontSize',11);   
   end 
   hold on ; 
end



%% 2) Progressive search algorithms 
% Initial values for P and N 
P0=2; 
N0=2; 
W0= 0.001.*ones(2*P0+1,1); 

% optimization parameters 
gradtype='EXACT';  
gradtype= 'FINITE_DIFFERENCES';
prec_loss= 10^-5; 
prec_grad= 10^-5;
eta = 10^-5; 
maxit= 5*10^3;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 2.1 Adaptive BFGS  
tic; 
[progress_BFGS_Wopt, progress_BFGS_IT, f_progress_BFGS, NP_evolution_BFGS] = Adaptive_BFGS(P0,N0,W0,prec_loss, prec_grad,eta,maxit,gradtype,Target); 
Adaptive_BFGS_time=toc;


%% 2.1 Adaptive PR
tic; 
[progress_PR_Wopt, progress_PR_IT, f_progress_PR, NP_evolution_PR] = Adaptive_PR(P0,N0,W0,prec_loss, prec_grad,eta,maxit,gradtype,Target); 
% progress_PR_Wopt =progress_BFGS_Wopt ;
% progress_PR_IT=progress_BFGS_IT;
% f_progress_PR= f_progress_BFGS;
% NP_evolution_PR= NP_evolution_BFGS;
Adaptive_PR_time=toc;

 
%% 2.1 Adaptive SD

tic; 
[progress_SD_Wopt, progress_SD_IT, f_progress_SD, NP_evolution_SD] = Adaptive_SD(P0,N0,W0,prec_loss, prec_grad,eta,maxit,gradtype,Target); 
% progress_SD_Wopt =progress_BFGS_Wopt ;
% progress_SD_IT=progress_BFGS_IT;
% f_progress_SD= f_progress_BFGS;
% NP_evolution_SD= NP_evolution_BFGS;
Adaptive_SD_time=toc; 





%% Evolution of the parameters N and P for BFGS
fN = NP_evolution_BFGS(progress_BFGS_IT,1); 
fP = NP_evolution_BFGS(progress_BFGS_IT,2);
[opt_bands_progress_BFGS, vectors] = Band_Structure_Progress(fN,fP,progress_BFGS_Wopt); 
Vbuild_Progress_BFGS = zeros(1,L); 
for x=1:L
    xc=Gamma(x); 
    for k=-fP:fP
        kp = k+fP+1; 
        Vbuild_Progress_BFGS(1,x) =Vbuild_Progress_BFGS(1,x)+progress_BFGS_Wopt(kp,1)*(exp(-1i*k*xc))/sqrt(2*pi);     
    end
end
FIG=figure(fig_pot);
hold on ;  
plot(Gamma, Vbuild_Progress_BFGS,'-.bs');
ylabel(gca,'Potentials','FontSize',14);
xlim([-3.15 3.15]); 





%% Evolution of the parameters N and P for PR
fN = NP_evolution_PR(progress_PR_IT,1); 
fP = NP_evolution_PR(progress_PR_IT,2);
maxNP= 2*fP+1;
maxNN= 2*fN+1; 
[opt_bands_progress_PR, vectors] = Band_Structure_Progress(fN,fP,progress_PR_Wopt); 
Vbuild_Progress_PR = zeros(1,L); 
for x=1:L
    xc=Gamma(x); 
    for k=-fP:fP
        kp = k+fP+1; 
        Vbuild_Progress_PR(1,x) =Vbuild_Progress_PR(1,x)+progress_PR_Wopt(kp,1)*(exp(-1i*k*xc))/sqrt(2*pi);     
    end
end
FIG=figure(fig_pot);
hold on ;  
plot(Gamma, Vbuild_Progress_PR,'-.ro');
ylabel(gca,'Potentials','FontSize',14);
xlim([-3.15 3.15]); 

%% Evolution of the parameters N and P for SD 
fN = NP_evolution_SD(progress_SD_IT,1); 
fP = NP_evolution_SD(progress_SD_IT,2); 
[opt_bands_progress_SD, vectors] = Band_Structure_Progress(fN,fP,progress_SD_Wopt); 
Vbuild_Progress_SD = zeros(1,L); 
for x=1:L
    xc=Gamma(x); 
    for k=-fP:fP
        kp = k+fP+1; 
        Vbuild_Progress_SD(1,x) =Vbuild_Progress_SD(1,x)+progress_SD_Wopt(kp,1)*(exp(-1i*k*xc))/sqrt(2*pi);     
    end
end
FIG=figure(fig_pot);
hold on ;  
plot(Gamma, Vbuild_Progress_SD,'-.m^');
ylabel(gca,'Potentials','FontSize',14);
xlim([-3.15 3.15]);
legend('Target', 'W0','BFGS', 'PR', 'SD','Progressive BFGS','Progressive PR','Progressive SD');






%% Plot the bands of the different algos 
FIG = figure(fig_bands);
set (FIG , 'name', 'Target-optimal Energy bands' ) ;
hold on ; 
for r=1:M
   subplot(M,1,M-r+1);
   plot(Brillouin,opt_bands_progress_BFGS(r,:),'-.bs',...
        Brillouin,opt_bands_progress_PR(r,:),'-.ro',...
         Brillouin,opt_bands_progress_SD(r,:),'-.m^'); 
   ylabel(Str, 'FontSize',11); 
   if (r==1)
       xlabel(' q \in  \Gamma^{*}', 'FontSize',11);   
   end 
   hold on ; 
end
legend('BFGS', 'PR', 'SD','Progressive BFGS','Progressive PR','Progressive SD');






%% Plot the history of Cost and grad   
History_cost_mat =  strcat(wkspace,'/BFGS/',wkspace,'_History_Cost_progress.mat'); 
Cost_BFGS_Progress= importdata(History_cost_mat);
History_cost_mat =  strcat(wkspace,'/PR/',wkspace,'_History_Cost_progress.mat'); 
Cost_PR_Progress= importdata(History_cost_mat);
History_cost_mat =  strcat(wkspace,'/SD/',wkspace,'_History_Cost_progress.mat'); 
Cost_SD_Progress= importdata(History_cost_mat);

History_cost_mat =  strcat(wkspace,'/BFGS/',wkspace,'_History_Cost.mat'); 
Cost_BFGS= importdata(History_cost_mat);
History_cost_mat =  strcat(wkspace,'/PR/',wkspace,'_History_Cost.mat'); 
Cost_PR= importdata(History_cost_mat);
History_cost_mat =  strcat(wkspace,'/SD/',wkspace,'_History_Cost.mat'); 
Cost_SD= importdata(History_cost_mat);

fig=fig+1;
FIG=figure(fig);
plot(Cost_BFGS,'-bs'); hold on;
plot(Cost_PR,'-ro');  hold on;
plot(Cost_SD,'-m^'); hold on;
plot(Cost_BFGS_Progress,'-.bs'); hold on;
plot(Cost_PR_Progress,'-.ro'); hold on;
plot(Cost_SD_Progress,'-.m^');  hold on;
legend('BFGS', 'PR', 'SD','Progressive BFGS','Progressive PR','Progressive SD');

path_jpg =  sprintf('%s/ConvergenceN%dP%d.jpg',wkspace,N,P);
saveas(FIG,path_jpg);
path_fig =  sprintf('%s/ConvergenceN%dP%d.fig',wkspace,N,P);
saveas(FIG,path_fig);


for f= 1:length(out_files)
    fil = out_files(f)
fprintf(fil,'\n ----- CPU times:-------\n'); 
fprintf(fil,'\n Naive SD = %.5f', SD_time); 
fprintf(fil,'\n Naive PR = %.5f', PR_time);
fprintf(fil,'\n Naive BFGS = %.5f', BFGS_time);
fprintf(fil,'\n Adaptive SD = %.5f', Adaptive_SD_time); 
fprintf(fil,'\n Adaptive PR = %.5f', Adaptive_PR_time);
fprintf(fil,'\n Adaptive BFGS = %.5f', Adaptive_BFGS_time);


fprintf(fil,'\n ----- Relative CPU times:-------\n'); 
fprintf(fil,'\n Naive SD = %.5f', SD_time/SD_time); 
fprintf(fil,'\n Naive PR = %.5f', PR_time/SD_time);
fprintf(fil,'\n Naive BFGS = %.5f', BFGS_time/SD_time);
fprintf(fil,'\n Adaptive SD = %.5f', Adaptive_SD_time/SD_time); 
fprintf(fil,'\n Adaptive PR = %.5f', Adaptive_PR_time/SD_time);
fprintf(fil,'\n Adaptive BFGS = %.5f', Adaptive_BFGS_time/SD_time);
end 













%% Plot the evolution of N with iterations 
fig=fig+1; 
FIG = figure(fig);
set (FIG , 'name', 'Evolution of s with iterations');
plot(2.*NP_evolution_BFGS(:,1)+1, '-bs'); hold on ; 
plot(2.*NP_evolution_PR(:,1)+1, '-r^'); hold on ; 
plot(2.*NP_evolution_SD(:,1)+1, 'color',[0.5 0.5 0]); hold on ; 
legend('progressive BFGS',  'progressive PR' , 'progressive SD'); 
xlabel(gca,'iterations','FontSize',14);
ylabel(gca,'s','FontSize',14);
path_jpg =  sprintf('%s/Evolution_N_iterations.jpg',wkspace); 
saveas(FIG,path_jpg); 
path_fig =   sprintf('%s/Evolution_N_iterations.fig',wkspace); 
saveas(FIG,path_fig);



%% Plot the evolution of P with iterations
fig=fig+1;
FIG = figure(fig);
set (FIG , 'name', 'Evolution of p with iterations');
plot(2.*NP_evolution_BFGS(:,2)+1, '-bs'); hold on ; 
plot(2.*NP_evolution_PR(:,2)+1, '-r^'); hold on ;
plot(2.*NP_evolution_SD(:,2)+1, 'color',[0.5 0.5 0]); hold on ;
legend('progressive BFGS',  'progressive PR' , 'progressive SD');
xlabel(gca,'iterations','FontSize',14);
ylabel(gca,'p','FontSize',14);
path_jpg =  sprintf('%s/Evolution_P_iterations.jpg',wkspace); 
saveas(FIG,path_jpg); 
path_fig =   sprintf('%s/Evolution_P_iterations.fig',wkspace); 
saveas(FIG,path_fig);



















fclose(texFile); 