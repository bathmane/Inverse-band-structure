clear; 
clc; 
close all; 
fil =1;
fprintf(fil,'\n\t **********************************************************************');
fprintf(fil,'\n\t **                  Code:     REVISED  code for  Flat Band        ****');
fprintf(fil,'\n\t **                Author:     Athmane Bakhta                      ****');
fprintf(fil,'\n\t **               Contact:                                         ****');
fprintf(fil,'\n\t **           Last update:     Dec  2018                           ****');
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
N=30;                   
P=20;                    
M=1;                    
Q=15; 
L=60; 
full_N=1000; 

%% Create workspaces 
dirp = P; 
Nb_mod ='Inverse_Flat_'; 
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


%% Target bands 
Target = zeros(M,Q);
Target(1,:)= 1; 


%%%=======================================================
%%%         N A I V E     O P T I M  
%%%=======================================================
W0 = zeros(2*P+1,1);
Vbuild_W0 = zeros(1,L); 
for x=1:L
    xc=Gamma(x); 
    for k=-P:P
        kp = k+P+1; 
        Vbuild_W0(1,x) =Vbuild_W0(1,x)+W0(kp,1)*(exp(-1i*k*xc))/sqrt(2*pi);     
    end
end
gradtype = 'FINITE_DIFFERENCES';
prec_loss = 10^-10; 
prec_grad = 10^-8; 
maxit     = 10000;


%% Polak Ribiere 
tic; 
[PR_Wopt, PR_IT, f_PR]=Naive_PR(W0,prec_loss, prec_grad,maxit,gradtype,Target);
PR_time = toc; 


%% Steepest Descent
tic; 
[SD_Wopt, SD_IT, f_SD] = Naive_SD(W0,prec_loss, prec_grad,maxit,gradtype,Target); 
SD_time = toc; 


%% Results of Naive search
[opt_bands_PR,get_eigen_vectors_V]=Band_Structure(PR_Wopt); 
[opt_bands_SD,get_eigen_vectors_V]=Band_Structure(SD_Wopt); 
% Build the potentials 
Vbuild_pr = zeros(1,L);
Vbuild_sd = zeros(1,L); 
Vbuild_orig = zeros(1,L); 
for x=1:L
    xc=Gamma(x); 
    for k=-P:P
         kp = k+P+1; 
         Vbuild_pr(1,x) =Vbuild_pr(1,x)+PR_Wopt(kp,1)*(exp(-1i*k*xc))/sqrt(2*pi);  
         Vbuild_sd(1,x) =Vbuild_sd(1,x)+SD_Wopt(kp,1)*(exp(-1i*k*xc))/sqrt(2*pi);    
    end
end 


fig_pot=100;
FIG=figure(fig_pot);
plot(Gamma,Vbuild_W0,'--k',...
        Gamma,Vbuild_pr, '-ro',...    
        Gamma,Vbuild_sd,'-m^'); 
xlim([-3.15 3.15]);
ylabel(gca,'Potentials','FontSize',14);
xlim([-3.15 3.15]);
legend('W0','PR', 'SD');


%%%=======================================================
%%%         A D A P T I V E     O P T I M  
%%%=======================================================
P0=3; 
N0=8; 
W0= 0.001.*ones(2*P0+1,1); 
gradtype= 'FINITE_DIFFERENCES';
prec_loss= 10^-10; 
prec_grad= 10^-8;
eta   = 10^-8; 
maxit_inter = 5000;



%% Adaptive PR
tic; 
[progress_PR_Wopt, progress_PR_IT, f_progress_PR, NP_evolution_PR] = Adaptive_PR(P0,N0,W0,prec_loss, prec_grad,eta,maxit_inter,gradtype,Target); 
Adaptive_PR_time=toc;

%% Adaptive BFGS  
tic; 
[progress_BFGS_Wopt, progress_BFGS_IT, f_progress_BFGS, NP_evolution_BFGS] = Adaptive_BFGS(P0,N0,W0,prec_loss, prec_grad,eta,maxit_inter,gradtype,Target); 
% progress_BFGS_Wopt = progress_PR_Wopt; 
% progress_BFGS_IT = progress_PR_IT; 
% f_progress_BFGS = f_progress_PR; 
Adaptive_BFGS_time=toc;

%% Adaptive SD
tic; 
%[progress_SD_Wopt, progress_SD_IT, f_progress_SD, NP_evolution_SD] = Adaptive_SD(P0,N0,W0,prec_loss, prec_grad,eta,maxit_inter,gradtype,Target); 
progress_SD_Wopt = progress_PR_Wopt; 
progress_SD_IT = progress_PR_IT; 
f_progress_SD = f_progress_PR; 
Adaptive_SD_time=toc; 


return 

%%%=======================================================
%%%         P L O T    R E S U L T S  
%%%=======================================================

%BFGS
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
%PR
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
%SD 
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


fig_pot=100;
FIG=figure(fig_pot);
plot(Gamma,Vbuild_W0,'--k',...
        Gamma,Vbuild_bfgs,'-bs',...
        Gamma,Vbuild_pr, '-ro',...    
        Gamma,Vbuild_sd,'-m^', ...
        Gamma, Vbuild_Progress_SD,'-.m^',...
        Gamma, Vbuild_Progress_BFGS,'-.bs',... 
        Gamma, Vbuild_Progress_PR,'-.ro'); hold on; 
xlim([-3.15 3.15]);
ylabel(gca,'Potentials','FontSize',14);
xlim([-3.15 3.15]);
legend('W0','BFGS', 'PR', 'SD','Progressive BFGS','Progressive PR','Progressive SD');
path_jpg =  sprintf('%s/PotentialsN%dP%d.jpg',wkspace,N,P);
saveas(FIG,path_jpg);
path_fig =  sprintf('%s/PotentialsN%dP%d.fig',wkspace,N,P);
saveas(FIG,path_fig);

%% Plot the bands of the different algos 
fig_bands=200;
FIG = figure(fig_bands);
set (FIG , 'name', 'Target-optimal Energy bands' ) ;
hold on ;
for r=1:M
    plot(Brillouin,Target(r,:),'-k',...
        Brillouin,opt_bands_BFGS(r,:),'-bs',...
        Brillouin,opt_bands_PR(r,:),'-ro',...
        Brillouin,opt_bands_SD(r,:),'-m^',...
        Brillouin,opt_bands_progress_BFGS(r,:),'-.bs',...
        Brillouin,opt_bands_progress_PR(r,:),'-.ro',...
        Brillouin,opt_bands_progress_SD(r,:),'-.m^');
    ylabel('Bands', 'FontSize',11);
    xlabel(' q \in  \Gamma^{*}', 'FontSize',11);
    hold on;
end
legend('Target','BFGS', 'PR', 'SD','Progressive BFGS','Progressive PR','Progressive SD');
path_jpg =  sprintf('%s/BandssN%dP%d.jpg',wkspace,N,P);
saveas(FIG,path_jpg);
path_fig =  sprintf('%s/BandssN%dP%d.fig',wkspace,N,P);
saveas(FIG,path_fig);


%% Plot the Cost function  
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
loglog(Cost_BFGS,'-bs'); hold on;
loglog(Cost_PR,'-ro');  hold on;
loglog(Cost_SD,'-m^'); hold on;
loglog(Cost_BFGS_Progress,'-.bs'); hold on;
loglog(Cost_PR_Progress,'-.ro'); hold on;
loglog(Cost_SD_Progress,'-.m^');  hold on;
xlim([0 maxit]); 
legend('BFGS', 'PR', 'SD','Progressive BFGS','Progressive PR','Progressive SD');
path_jpg =  sprintf('%s/ConvergenceN%dP%d.jpg',wkspace,N,P);
saveas(FIG,path_jpg);
path_fig =  sprintf('%s/ConvergenceN%dP%d.fig',wkspace,N,P);
saveas(FIG,path_fig);

% Print cpu times 
for f= 1:length(out_files)
    fil = out_files(f);
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
semilogx(NP_evolution_BFGS(:,1), '-bs'); hold on ; 
semilogx(NP_evolution_PR(:,1), '-r^'); hold on ; 
semilogx(NP_evolution_SD(:,1), 'color',[0.5 0.5 0]); hold on ; 
xlim([0 maxit]);
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
semilogx(NP_evolution_BFGS(:,2), '-bs'); hold on ; 
semilogx(NP_evolution_PR(:,2), '-r^'); hold on ;
semilogx(NP_evolution_SD(:,2), 'color',[0.5 0.5 0]); hold on ;
xlim([0 maxit]);
legend('progressive BFGS',  'progressive PR' , 'progressive SD');
xlabel(gca,'iterations','FontSize',14);
ylabel(gca,'p','FontSize',14);
path_jpg =  sprintf('%s/Evolution_P_iterations.jpg',wkspace); 
saveas(FIG,path_jpg); 
path_fig =   sprintf('%s/Evolution_P_iterations.fig',wkspace); 
saveas(FIG,path_fig);
fclose(texFile); 