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

 

%%  global vars 
global Brillouin;   % Discrete 1st Brillouin zone
global Gamma;       % Gamma: unit cell of the lattice 
global Q;           % Number of k-points in Brillouin zone
global N;           % Numbe of basis functions
global P;           % Number of Fourier modes in the potential 
global L;           % Number of x points in Gamma
global M;           % Number of target bands in the cost function 
global fig;
global full_N;      % Used in the progressive algorithm 


fig=1; 




%% Initializations
N=5;           
P=2;         
Q=10; 
L=60;  
full_N=1000; 
M= 3; 
 



%% Discretization of the first Brillouin zone 
left= 0;  
right=0.5;
Brillouin = left:(right-left)/(Q-1):right; 
  
%% Discretization of Gamma 
left=-pi;  
right=pi;
Gamma = left:(right-left)/(L-1):right ; 


%% The periodic potential Vper(x)  
Vper = zeros(2*P+1,1); 
for k=1:P
    Vper(k,1)= rand();
    Vper(2*P+1-k+1,1) = Vper(k,1);
end
Vbuild_Vper = zeros(1,L); 
for x=1:L
    xc=Gamma(x); 
    for k=-P:P
        kp = k+P+1; 
        Vbuild_Vper(1,x) =Vbuild_Vper(1,x)+Vper(kp,1)*(exp(-1i*k*xc))/sqrt(2*pi);     
    end
end

 

fig=fig+1;
FIG=figure(fig);
set(FIG, 'name', 'Periodic Potential Vper');
hold on ;  
plot(Gamma,Vbuild_Vper,'k', 'linewidth', 1.1);
xlabel(gca,'x \in \Gamma ','FontSize',14);
ylabel(gca,'V_{per}','FontSize',14);
 
%% To be sure that The spectrum is positive
Vper


Vper(P+1,1) = - real(min(Vbuild_Vper(1,:)))+1;  
Vbuild_Vper = Vbuild_Vper + Vper(P+1,1); 

Vper

fig=fig+1;
FIG=figure(fig);
set(FIG, 'name', 'Periodic Potential with nonnegative spectrum');
hold on ;  
plot(Gamma,Vbuild_Vper,'k', 'linewidth', 1.1);
xlabel(gca,'x \in \Gamma ','FontSize',14);
ylabel(gca,'V_{per}','FontSize',14);


%% Compute the band structure of Vper
[get_bands_V, get_eigen_vectors_V] = Band_Structure(Vper); 

fig=fig+1;
FIG = figure(fig);
set (FIG , 'name', 'Target Energy bands' ) ;
hold on ; 
for r=1:M
   plot(Brillouin,get_bands_V(r,:));   
   ylabel('Bands', 'FontSize',11); 
   hold on ; 
end
xlabel(' q \in  \Gamma^{*}', 'FontSize',11);   



%% Exact errors 
Exact_Err = zeros(M,Q); 

tic; 
[exact_bands, exact_vectors] = Band_Structure_Progress(full_N,P,Vper);
E1 = 0;
for m=1:M
    for q=1:Q
        Exact_Err(m,q) = abs(exact_bands(m,q) - get_bands_V(m,q));
    end
end
exact_time = toc; 


% 
% fig=fig+1;
% FIG = figure(fig);
% set (FIG , 'name', 'EXACT BANDS' ) ;
% hold on ; 
% for r=1:M
%    plot(Brillouin,exact_bands(r,:));   
%    ylabel('Bands', 'FontSize',11); 
%    hold on ; 
% end
% xlabel(' q \in  \Gamma^{*}', 'FontSize',11);   
% 
% 
% 
% fig=fig+1;
% FIG = figure(fig);
% set (FIG , 'name', 'Exact Errors' ) ;
% hold on ; 
% for r=1:M
%    plot(Brillouin,Exact_Err(r,:), '-ko');   
%    hold on ; 
% end


%% Aposteriori errors 
tic; 
Post_Err = Aposteriori_Error(get_bands_V,get_eigen_vectors_V,Vper,P,N,full_N); 
estim_time = toc; 
 
tic; 
Post_Err_2 = Aposteriori_Error_2(get_bands_V,get_eigen_vectors_V,Vper,P,N,full_N); 
estim_time_2 = toc; 

% fig=fig+1;
% FIG = figure(fig);
% set (FIG , 'name', 'Apost Errors' ) ;
% hold on ; 
% for r=1:M
%    plot(Brillouin,Post_Err(r,:), '-rs');   
%    hold on ; 
% end
% xlabel(' q \in  \Gamma^{*}', 'FontSize',11);
% ylabel('estimators', 'FontSize',11); 

fprintf('\n\n ------------ CPU times comparison ---------\n'); 
disp('Exact error time '); 
disp(exact_time);
disp('Estim error with matrix inversion '); 
disp(estim_time); 
disp('Estim bis  error with BiCGSTAB'); 
disp(estim_time_2); 
%% Comparison 
fig=fig+1;
FIG = figure(fig);
set (FIG , 'name', 'Error for band 1' ) ;
hold on ; 
plot(Brillouin,Exact_Err(1,:), 'k');  
plot(Brillouin,Post_Err(1,:), 'r'); 
plot(Brillouin,Post_Err_2(1,:), 'm'); 
legend('Exact error', ' Estimator (matrix inversion)', 'Estimator (BiCGSTAB)'); 



fig=fig+1;
FIG = figure(fig);
set (FIG , 'name', 'Error for band 2 ' ) ;
hold on ; 
plot(Brillouin,Exact_Err(2,:), 'k');  
plot(Brillouin,Post_Err(2,:), 'r'); 
plot(Brillouin,Post_Err_2(2,:), 'm'); 
legend('Exact error', ' Estimator (matrix inversion)', 'Estimator (BiCGSTAB)'); 

fig=fig+1;
FIG = figure(fig);
set (FIG , 'name', 'Error for band 3' ) ;
hold on ; 
plot(Brillouin,Exact_Err(3,:), 'k');  
plot(Brillouin,Post_Err(3,:), 'r'); 
plot(Brillouin,Post_Err_2(3,:), 'm'); 
legend('Exact error', ' Estimator (matrix inversion)', 'Estimator (BiCGSTAB)'); 




fprintf('\n-- -- -- -- -- *** END *** -- -- -- -- --'); 