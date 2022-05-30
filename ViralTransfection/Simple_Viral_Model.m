p0 = [0.99, 0.0001];
t0 = -5;
tfinal = 30;

[t,p] = ode45(@viral,[t0 tfinal], p0);

% to graph viral load:
semilogy(t,p(:,2))
xlabel('days')
ylabel('viral load (RNA copies/mL)')

% to graph fraction of uninfected cells
%semilogy(t,p(:,1))
%xlabel('days')
%ylabel('fraction of uninfected cells)')

function dpdt = viral(t,p)
%p(1) = f(t)
%p(2) = V(t)

beta = 7.8e-6; % rate constant for virus infection
gamma = 3.91; % maximum rate constant for viral replication
delta = 0.53; % death rate of infected cells 


dpdt = [-beta*p(1)*p(2);
    gamma*p(1)*p(2) - delta*p(2)];
end 