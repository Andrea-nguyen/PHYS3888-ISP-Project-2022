function [V,N,t] = ImmuneViralInfection(tspan,logyV,logyN,V0,N0,gamma,K,cN,r,m)
%% Immune response to viral infection MODEL
% This m-file numerically integrates the (planar) equations corresponding
% to the ODE:
%
% f = V^m / (V^m + KN^m)
%
% V' = gamma*V*(1 - V/K) - cN*V*N - c*V;
% N' = sN + r*N*f - deltaN*N;
%
% Where  V(t), N(t) >= 0 are variables denoting Viral concentration and T-cell
% population at time t correspondingly.
%
% The file then plots 5 different figures.
%
% Figure 1: a plot of the time series V.
% Figure 2: a plot of the time series N.
% Figure 3: a plot of the time series of the given solution.
% Figure 4: a plot of the time series of the given solution on semilogX scale.
% Figure 5: a plot of the solution trajectory in (V,N) phase space.
%
% Parameters:
%   tspan - integration time
%   logyV - logical to plot V(t) in semilogY
%   logyN - logical to plot N(t) in semilogY
%   V0 - Initial value for viral load V(0)
%   N0 - Initial value for T-cell T(0)
%   gamma - grow rate of virus
%   K - maximum concentration of viral load
%   cN - killing rate of virus by T-cells
%   r - Proliferation rate of T-cells
%   m - width of sigmoidal function f
%
% Returns: V(t), N(t), t


%% Check inputs and set defaults:
if nargin < 1 || isempty(tspan)
    tspan = [0 150];
end
if nargin < 2 || isempty(logyV)
    logyV = false;
end
if nargin < 3 || isempty(logyN)
    logyN = false;
end
if nargin < 4 || isempty(V0)
    V0 = 0.31;
end
if nargin < 5 || isempty(N0)
    N0 = 1e6;
end
if nargin < 6 || isempty(gamma)
    gamma = 3.98;
end
if nargin < 7 || isempty(K)
    K = 7.11e8;
end
if nargin < 8 || isempty(cN)
    cN = 1.58e-8;
end
if nargin < 9 || isempty(r)
    r = 0.794;
end
if nargin < 10 || isempty(m)
    m = 2;
end


%% PREAMBLE

close all; %Closes all figure windows.

format long; %Displays data to 15 decimal places

set(0,'defaultTextInterpreter','latex');
set(0,'defaultAxesFontSize',15);



%% INITIALISATION

% Parameter values

% gamma = 3.98;
% K = 7.11e8;
% cN = 1.58e-8;
c = 2.4;

deltaN = 0.1;
sN = N0 * deltaN;
% r = 0.794;
KN = 1e3;
%m = 2;

% Integration time interval
%tspan = [0 150];

% Initial condition (point in the plane) z0 = (V(0),N(0)) = (V0,N0)

z0 = [V0, N0];





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
opts = odeset('RelTol',1e-3,'AbsTol',1e-6);

%Integration step
[t,z] = ode15s(@(t,z) odefcn(t,z,gamma,K,cN,c,deltaN,sN,r,KN,m), tspan, z0,opts);


V = z(:,1); % the first column stores the x components along the discretised solution curve
N = z(:,2); % the second column stores the y components along the discretised solution curve







%% FIGURE PLOTTING


%Figure 1: the time series V(t).

figure(1); hold on; 

plot(t,V,'b','LineWidth',1); %Plots the time series for the x-component in blue

xlabel('$t$ (days)');  %labels the x-axis
ylabel('Viral concentration (copies/mL)')
legend('$V(t)$','interpreter','latex'); %adds a legend
%grid on; %adds a grid to guide the eye
if logyV
    set(gca, 'YScale', 'log')
end
movegui('northwest'); %places the figure on the top left.



%Figure 2: the time series N(t).

figure(2); hold on; 

plot(t,N,'r','LineWidth',1); %Plots the time series for the y-component in red
ylabel('T-cell population (cells)')
xlabel('$t$ (days)');  %labels the x-axis
legend('$N(t)$','interpreter','latex'); %adds a legend
%grid on; %adds a grid to guide the eye
if logyN
    set(gca, 'YScale', 'log')
end
movegui('northeast'); %places the figure on the top left.


%Figure 3: the time series.
figure(3); hold on; %opens up a figure window. The 'hold on' command allows us to plot multiple
%objects within one window.

plot(t,V,'b','LineWidth',1); %Plots the time series for the x-component in blue
plot(t,N,'r','LineWidth',1); %Plots the time series for the y-component in red
ylabel('State')
xlabel('$t$ (days)');  %labels the x-axis
legend('$V(t)$','$N(t)$','interpreter','latex'); %adds a legend

%grid on; %adds a grid to guide the eye

movegui('west'); %places the figure on the top left.


%Figure 4: the time series on semilog.

figure(4); hold on; 

plot(t,V,'b','LineWidth',1); %Plots the time series for the x-component in blue
plot(t,N,'r','LineWidth',1); %Plots the time series for the y-component in red
ylabel('State')
xlabel('$t$ (days)');  %labels the x-axis
legend('$V(t)$','$N(t)$','interpreter','latex'); %adds a legend
set(gca, 'YScale', 'log')
%grid on; %adds a grid to guide the eye

movegui('east'); %places the figure on the top left.


%Figure 5: the phase space plot.

figure(5); hold on; 

plot(V,N,'k','LineWidth',1); %Plots the solution trajectory in black in phase space.

plot(V(1),N(1),'ko','MarkerSize',8,'MarkerFaceColor','g'); %Plots the initial condition as a green dot.
plot(V(end),N(end),'ko','MarkerSize',8,'MarkerFaceColor','m'); %Plots the final point in the trajectory as a magenta dot.


xlabel('Viral concentration (copies/mL)');  
ylabel('T-cell population (cells)');  
grid on; 
set(gca, 'XScale', 'log')

movegui('south'); 


%Figure 5: the phase space plot.
% 
% figure(5); hold on; 
% 
% plot(V,N,'k','LineWidth',1); %Plots the solution trajectory in black in phase space.
% 
% plot(V(1),N(1),'ko','MarkerSize',8,'MarkerFaceColor','g'); %Plots the initial condition as a green dot.
% plot(V(end),N(end),'ko','MarkerSize',8,'MarkerFaceColor','m'); %Plots the final point in the trajectory as a magenta dot.
% 
% 
% xlabel('Viral concentration (copies/mL)');  
% ylabel('T-cell population (cells)');  
% grid on; 
% %set(gca, 'XScale', 'log')
% 
% movegui('north'); 







%% DEFINITION OF THE ODEs


%A warning shows up for 't' to let us know that the time variable is unused
%in the equations (i.e. they are *autonomous*). MATLAB also allows us to
%integrate nonautonomous equations (for eg. to study problems with time-dependent
%forcing.)

function dzdt = odefcn(t,z,gamma,K,cN,c,deltaN,sN,r,KN,m)
dzdt = zeros(2,1);

%These are the equations:
% f = V^m/(V^m+KN^m)
%
% V' = gamma*V*(1-V/K)-cN*V*N-c*V;
% N' = sN + r*N*f - deltaN*N;

f = z(1)^m/(z(1)^m+KN^m);
     
dzdt(1) = gamma*z(1)*(1-z(1)/K)-cN*z(1)*z(2)-c*z(1);
dzdt(2) = -deltaN*z(2) + r*z(2)*f + sN;
end
end