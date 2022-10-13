clear; clc; close('all');
% inputs
D = 1 ; % metre 
ray = D/2 ; 
Ce = 10 ; % mol/m3 
k = 1.8*10^(-9) ; % s^(-1) 
Deff = 10^(-10) ; % m2/s 
s = 10^(-8) ; % mol/m3/s 
Nr = 100; 
tsim = 10^10; % secondes
Nt = 100; 
annees = tsim/(86400*365);

% steps 
ri = 0 ; rf = D/2  ; dr = (rf - ri)/Nr; 
ti = 0;  tf = tsim ; dt = (tf-ti)/Nt ;

h = dr; 
R = ri:dr:rf; nr = numel(R) ; 
T = ti:dt:tf; nt = numel(T) ;

[r,t] = meshgrid(R,T); r = r'; t = t'; % x and t matrices


% Concentration Matrix
%c = ones(nr,nt)*(h^2*(1+k*dt));
c = zeros(nr,nt);
% c(1,1) = 0;
% c(nr,1) = Ce; 
 

m = zeros(nr,nr);
m(1,1) = -3 ;
m(1,2) = 4 ; 
m(1,3) = -1 ;
m(nr,nr) = 1; 

for i=2:Nr  
 alpha = (-Deff*dt*h/(2*R(i))-Deff*dt) / (h^2*(1-k*dt)); 
 beta = (h^2+2*Deff*dt)/ (h^2*(1-k*dt)); 
 gamma = (Deff*dt*h/(2*R(i))-Deff*dt)/ (h^2*(1-k*dt)); 
    m(i,i) = beta ;
    m(i,i+1) = alpha ; 
    m(i,i-1) = gamma ;   
end

matcoeffs= [ alpha; beta ; gamma];


for t=1:nt-1
    K = [0 ; c(2:nr-1,t) ; Ce] ;  % Cold
    X = m\K ; 
    c(:,t+1) = X;
end     


% solution analytique 
C_analytique = zeros(nr,1); 
for i=1:nr
C_analytique (i) =  0.25 * s/Deff * ray^2 * (R(i)^2/ray^2 -1) + Ce; 
end 


%% Results

% Animation
    for j = 1:nt
        plot(R,c(:,j)); xlabel('Distance [m]'); ylabel('Concentration [mol]');
        title(['Concentration Versus Distance at Time = ' num2str(round(T(j),3)) ' s']) ; 
        grid('on'); drawnow; %pause(1);
    end

    figure(2)
    plot(R,C_analytique)
    hold on 
    plot(R,c(:,nt),'o')
    grid on 
    title(['Concentration Versus Distance ' num2str(round(T(j),3)) ' s']) ;
    legend('solution analytique', 'solution numérique')
    
    
    