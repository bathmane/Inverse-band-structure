clear all; 
clc; 
close all; 

fil =1;
fprintf(fil,'\n\t **********************************************************************');
fprintf(fil,'\n\t **                  Code:     Inverse Hill                       ****');
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
global wkspace; 
global fig;
global full_N;      % Used in the adaptive algorithm 


fig=1; 




%% Initializations
N=10;           
P=2;          
Q=51; 
L=60; 
maxN=40;
maxP=20; 
full_N=30; 

%% Create workspaces 
wkspace='Band_structure';
mkdir(wkspace); 



%% Discretization of the first Brillouin zone 
left=-0.5;  
right=0.5;
Brillouin = left:(right-left)/(Q-1):right; 
  
%% Discretization of Gamma 
left=-pi;  
right=pi;
Gamma = left:(right-left)/(L-1):right ; 


%% The periodic potential Vper(x)  
Vper = zeros(2*P+1,1); 
for k=1:P
    Vper(k,1)= rand()+5*rand()*1i;
    Vper(2*P+1-k+1,1) = conj( Vper(k,1)  );
end
%Vper = 0.*Vper; 

Vbuild_Vper = zeros(1,L); 
for x=1:L
    xc=Gamma(x); 
    for k=-P:P
        kp = k+P+1; 
        Vbuild_Vper(1,x) =Vbuild_Vper(1,x)+Vper(kp,1)*(exp(-1i*k*xc))/sqrt(2*pi);     
    end
end
 
%%% To be sure that The spectrum is positive 
Vper(P+1,1) = floor(abs(min(real(Vbuild_Vper(1,:))) +1)) ; 
Vbuild_Vper = Vbuild_Vper + Vper(P+1,1);  

fig=fig+1;
FIG=figure(fig);
set(FIG, 'name', 'Periodic Potential Vper');
hold on ;  
plot(Gamma,Vbuild_Vper,'k', 'linewidth', 1.1);
xlabel(gca,'x \in \Gamma ','FontSize',14);
ylabel(gca,'V_{per}','FontSize',14);


path_jpg =  sprintf('%s/Pot_P%d.jpg',wkspace,P); 
saveas(FIG,path_jpg); 
path_fig =  sprintf('%s/Pot_P%d.fig',wkspace,P); 
saveas(FIG,path_fig);   
 


%% Compute the band structure of Vper
[get_bands_V, get_eigen_vectors_V] = Band_Structure(Vper); 
 
M= 4; 
fig=fig+1;
FIG = figure(fig);
set (FIG , 'name', 'Target Energy bands' ) ;
hold on ; 
for r=1:M
   subplot(M,1,M-r+1);
   plot(Brillouin,get_bands_V(r,:), 'k', 'linewidth', 1.5);   
   Str = sprintf('band %d', r); 
   ylabel(Str, 'FontSize',11); 
   if (r==1)
       xlabel(' q \in  \Gamma^{*}', 'FontSize',11);   
   end 
   hold on ; 
end

path_jpg=sprintf('%s/Bands_P%d.jpg',wkspace,P); 
saveas(FIG,path_jpg); 
path_fig=sprintf('%s/Bands_P%d.fig',wkspace,P); 
saveas(FIG,path_fig);   
 





fprintf('\n-- -- -- -- -- *** END *** -- -- -- -- --'); 