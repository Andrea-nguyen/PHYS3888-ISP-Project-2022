function [T,N,C,M,t] = modChemotherapy(T0,M0,treatmentType,treatment,tspan,xlimAll,logyT,logyAll,KN,KC,kT)
%% Modified Chemotherapy MODEL to incorporate treament length  
% This function numerically integrates the (planar) equations corresponding
% to the ODE:
%
% 
% T' = a*T*(1-b*T)-c1*N*T-KT*M*T;
% N' = alpha1 + g*T*N/(s+T) - mu*N - p*T*N - KN*M*N;
% C' = alpha2 - beta*C - KC*M*C;
% M' = -gamma*M + VM - kT*M*T;
%
%
% Where T,N,C,M >= 0 are variables denoting
%   T(t): Tumor cell population
%   N(t): Effector immune cell
%   C(t): Circulating lymphocyte population
%   M(t): Concentration of chemotherapeutic drugs in tissue
%
%The file then plots five different figures.
%
% Figure 1: a plot of the time series of T
% Figure 2: a plot of the time series of M
% Figure 3: a plot of the time series of N
% Figure 4: a plot of the time series of C
% Figure 6: a plot of the time series of all with normalised values
% Figure 5: a plot of the solution trajectory of the identical solution in (T,N)
%(phase) space.
%
% Parameters:
%   T0 - initial value of tumor cells
%   M0 - initial value of therapy concentration
%   treatmentType - type of treatment
%       'c' - constant treatment
%       'h' - spike injection modelled by heaviside function
%       'g' - injection modelled by gamma distribution
%   treatment = [VM, VMstart, VMend]
%       VM - height of pulse
%       VMstart - when the treatment starts
%       VMend - when the treatment ends
%   or treatment = [VM, a, b, s]
%       VM - 0 or 1, either np pulse or otherwise
%       a,b - parameters of gamma function
%       s - scaling factor so that peak is close to 1
%   logyT - logical to plot T(t) in semilogY
%   logyAll - logical to plot all time series in semilogY for figure 6
%   kT - drug and tumour combinational rate
%
% Returns: T(t), N(t), C(t), M(t), t


%% Check inputs and set defaults:
if nargin < 1 || isempty(T0)
    T0 = 1e7;
end
if nargin < 2 || isempty(M0)
    M0 = 0.5;
end
if nargin < 3 || isempty(treatmentType)
    treatmentType = 'g';
end
if nargin < 4 || isempty(treatment)
    treatment = [1, 2, 0.2, 14];
end
if nargin < 5 ||  isempty(tspan)
    tspan = [0 2000]; % output as matrix for colourmap instead
end
if nargin < 6 || isempty(xlimAll)
    xlimAll = tspan;
end
if nargin < 7 || isempty(logyT)
    logyT=false;
end
if nargin < 8 || isempty(logyAll)
    logyAll=false;
end
if nargin < 9 || isempty(KN)
    KN = 6e-1;
end
if nargin < 10 || isempty(KC)
    KC = 6e-1;
end
if nargin < 11 || isempty(kT)
    kT = 0;
end



%% PREAMBLE

%clear all; %Deletes all the data from the current workspace: good coding practice to ensure that
% %outdated values do not unknowingly affect your current run.
% 
close all; %Closes all figure windows.

format long; %Displays data to 15 decimal places

set(0,'defaultTextInterpreter','latex');
set(0,'defaultAxesFontSize',15);



%% INITIALISATION

% Parameter values
a = 4.31e-3;
b = 1.02e-14;
c1 = 3.41e-10;

alpha1 = 1.2e4;
alpha2 = 7.5e8;
beta = 1.2e-2;
%eta = 0; % changable, 0 <= eta <= 1
gamma = 9e-1;
mu = 4.12e-2;

g = 1.5e-2;
s = 2.02e1;
p = 2e-11;

%KC = 6e-1;
%KN = 6e-1;
KT = 8e-1;

%kT = 0; %3.2e-9;

%VM dependent on time, define later in the ode, changable, 0 <= VM <= 1

% Integration time interval
%tspan = [0 2000];

% Initial condition (point in the plane) z0 = (T(0),N(0),C(0),M(0)) = (T0,N0,C0,M0)

z0 = [T0 3e5 6.25e10 M0];





%% NUMERICAL INTEGRATION 

%The numerical integration step. The function 'ode45' takes as input
%arguments the ode function 'odefcn' (specified at the end of the file, where you actually 
%define the ODEs), the time interval
%tspan to be integrated, the initial condition z0, and the error tolerances
% specified in opts.

% The output are a set of arrays of equal length (but different
% widths): the discrete values for the time of integration t and the
% corresponding discrete solution values stored in z.

% Specify error tolerance for integration step
opts = odeset('RelTol',1e-12,'AbsTol',1e-12);

%Integration step
[t,z] = ode15s(@(t,z) odefcn(t,z,a,b,c1,KT,alpha1,g,s,mu,p,KN,alpha2,beta,KC,gamma,treatmentType,treatment,kT), tspan, z0, opts);


T = z(:,1); % the first column stores the x1 components along the discretised solution curve
N = z(:,2); % the second column stores the x2 components along the discretised solution curve
C = z(:,3); % the third column stores the x3 components along the discretised solution curve
M = z(:,4); % the forth column stores the x4 components along the discretised solution curve








%% FIGURE PLOTTING


% %Figure 1: the time series.
% figure(1); hold on; %opens up a figure window. The 'hold on' command allows us to plot multiple
% %objects within one window.
% 
% plot(t,T,'b','LineWidth',1); %Plots the time series for the x1-component in blue
% plot(t,M,'r','LineWidth',1); %Plots the time series for the x2-component in red
% 
% xlabel('$t$ (days)');  %labels the x-axis
% ylabel('Cell populations');
% legend('$T(t)$','$M(t)$','Interpreter',"latex"); %adds a legend
% grid on; %adds a grid to guide the eye
% 
% movegui('northwest'); %places the figure on the top left.



% %Figure 1: the time series T.
% 
% figure(1); hold on; 
% 
% plot(t,T,'color',[0 0.4470 0.7410],'LineWidth',1); %Plots the time series for the x1-component in blue
% ylabel('Tumour cell population $T(t)$ [cell]');
% xlabel('$t$ (days)');  %labels the x-axis
% legend('$T(t)$','Interpreter',"latex"); %adds a legend
% 
% %grid on; %adds a grid to guide the eye
% if logyT
% 	set(gca, 'YScale', 'log')
% end
% movegui('northeast'); %places the figure on the top right.
% 
% 
% 
% % Figure 2: the time series N.
% 
% figure(2); hold on; 
% 
% plot(t,N,'color',[0.4940 0.1840 0.5560],'LineWidth',1); %Plots the time series for the x2-component in black
% ylabel('Effector immune cell population $N(t)$ [cell]');
% xlabel('$t$ (days)');  %labels the x-axis
% legend('$N(t)$','Interpreter',"latex"); %adds a legend
% 
% %grid on; %adds a grid to guide the eye 
% movegui('west'); %places the figure on the top left.
%  
% 
% 
% % %Figure 3: the time series C.
% 
% figure(3); hold on; 
% 
% plot(t,C,'color',[0.8500 0.3250 0.0980],'LineWidth',1); %Plots the time series for the x3-component in green
% ylabel('Circulating lymphocyte population $C(t)$ [cell]');
% xlabel('$t$ (days)');  %labels the x-axis
% legend('$C(t)$','Interpreter',"latex"); %adds a legend
% %grid on; %adds a grid to guide the eye
% 
% movegui('east'); 
% 
% 
% 
%Figure 4: the time series M.

figure(4); hold on; 

plot(t,M,'color',[0.4660 0.6740 0.1880],'LineWidth',1); %Plots the time series for the x4-component in red
ylabel('Concentration of chemotherapeutic drugs $M(t)$')
xlabel('$t$ (days)');  %labels the x-axis
legend('$M(t)$','Interpreter',"latex"); %adds a legend
%grid on; %adds a grid to guide the eye
xlim([0 25]);
%xline(9.125)
movegui('northwest'); %places the figure on the left.



%Figure 5: the phase space plot.

figure(5); hold on; 

plot(T,M,'k','LineWidth',1); %Plots the solution trajectory in black in phase space.
 
plot(T(1),M(1),'ko','MarkerSize',8,'MarkerFaceColor','g'); %Plots the initial condition as a green dot.
plot(T(end),M(end),'ko','MarkerSize',8,'MarkerFaceColor','m'); %Plots the final point in the trajectory as a magenta dot.

xlabel('$T(t)$');  
ylabel('$M(t)$');


%grid on; 
movegui('north'); %places the figure on the right.
 



%Figure 6: the time series normalised by initial condition.

Tnorm = T ./ z0(1);
Nnorm = N ./ z0(2);
Cnorm = C ./ z0(3);
Mnorm = M;

figure(6); hold on;

plot(t,Tnorm,'LineWidth',2); %Plots the time series for the x1-component in blue
plot(t,Nnorm,'--','LineWidth',2); %Plots the time series for the x2-component in red
plot(t,Cnorm,'-.','LineWidth',2); %Plots the time series for the x3-component in red
plot(t,Mnorm,':','LineWidth',2); %Plots the time series for the x4-component in red

xlabel('$t$ (days)');  %labels the x-axis
ylabel('States (Normalised to initial conditions)');
%xlim([0 100]);
legend('Tumour $T(t)$','Effector immune cell $N(t)$', 'Circulating lymphocytes $C(t)$', 'Drug concentration $M(t)$','Interpreter',"latex"); %adds a legend
%grid on; %adds a grid to guide the eye
if logyAll
	set(gca, 'YScale', 'log')
else
    ylim([0 1.5]);
    yticks(0:0.5:1.5);
end
%xline(9.125)
xlim(xlimAll);
movegui('northwest'); %places the figure on the top left.




%Figure 7: the time series normalised by initial condition.

figure(7); hold on;

plot(t,T,'LineWidth',2); %Plots the time series for the x1-component in blue
plot(t,N,'--','LineWidth',2); %Plots the time series for the x2-component in red
%plot(t,Cnorm,'-.','LineWidth',1); %Plots the time series for the x3-component in red
%plot(t,Mnorm,':','LineWidth',2); %Plots the time series for the x4-component in red

xlabel('$t$ (days)');  %labels the x-axis
ylabel('Cells');
%xlim([0 100]);
%grid on; %adds a grid to guide the eye
%set(gca, 'YScale', 'log')
if logyT
	set(gca, 'YScale', 'log')
end
title('Initial values $(T_0, N_0, C_0, M_0) = (10^7, 3\times 10^5, 6.25 \times 10^{10}, 0)$')
legend('Tumour $T(t)$','Effector immune cell $N(t)$','Interpreter',"latex"); %adds a legend
%xline(9.125)
%legend('Tumour $T(t)$','Effector immune cell $N(t)$','$x=9.125$','Interpreter',"latex"); %adds a legend
%xlim(xlimAll);
movegui('northwest'); %places the figure on the top left.





%% DEFINITION OF THE ODEs


%A warning shows up for 't' to let us know that the time variable is unused
%in the equations (i.e. they are *autonomous*). MATLAB also allows us to
%integrate nonautonomous equations (for eg. to study problems with time-dependent
%forcing.)

function dzdt = odefcn(t,z,a,b,c1,KT,alpha1,g,s,mu,p,KN,alpha2,beta,KC,gamma,treatmentType,treatment,kT)
dzdt = zeros(4,1);

%These are the equations:

% Drug for a period of time
if treatmentType == 'c'
    VMmod = treatment(1);
elseif treatmentType == 'h'
    VMmod = treatment(1)*(heaviside(t-treatment(2))-heaviside(t-treatment(3)));
elseif treatmentType == 'g'
    VMmod = medModel(treatment(1), t, treatment(2),treatment(3),treatment(4));
end
% T' = a*T*(1-b*T)-c1*N*T-KT*M*T;
% N' = alpha1 + g*T*N/(s+T) - mu*N - p*T*N - KN*M*N;
% C' = alpha2 - beta*C - KC*M*C;
% M' = -gamma*M + VM - kT*M*T;

dzdt(1) = a*z(1)*(1-b*z(1)) - c1*z(2)*z(1) - KT*z(4)*z(1);
dzdt(2) = alpha1 + g*z(1)/(s+z(1))*z(2) - mu*z(2) - p*z(1)*z(2) - KN*z(4)*z(2);
dzdt(3) = alpha2 - beta*z(3) - KC*z(4)*z(3);
dzdt(4) = -gamma*z(4) + VMmod - kT*z(1)*z(4);

end

end