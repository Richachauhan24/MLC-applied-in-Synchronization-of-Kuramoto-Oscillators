%% Code for the weighted network of 6 nodes of Kuramoto oscillators with Kuramoto equation solutions showing simulations of the network.

clear all
clc

%% Network
N = 6; %% no. of nodes
s = [1 2 2 3 4 5 6 6 6];
t = [6 4 6 5 1 2 2 3 4];
weights = [10 5 5 10 9 9 7 2 2];
names = {'1','2','3','4','5','6'};
G = digraph(s,t,weights,names); %% directed strongly connected graph
A = adjacency(G,'weighted');
A = full(A); %% weighted adjacency matrix


%% Network
%N = 12; %% no. of nodes
%A=[ones(N/2) ones(N/2);ones(N/2) ones(N/2)] - eye(N); %% weighted adjacency matrix
%G = digraph(A); %% directed strongly connected graph


%% integration parameters
dt = 0.01;
maximum_step = 300/dt; %1000/dt



%% parameters
omega = ones(N,1)*1/100; % natural frequency parameters

%% initialization
theta = zeros (N, maximum_step); % phase initial state
theta(:,1) = 2* pi* rand (N,1); 
dthetadt = zeros (N,1);

for t = 1:maximum_step
    
    %% Main Integration
    coup = zeros (N,1);
   
    for i = 1:N
        coup (i) = sum(A(:,i).* sin(theta(:,t) - theta(i,t)));
        %coup (i) = 0.001*sum(A(:,i).* sin(theta(:,t) - theta(i,t))); % make weak coupling for desynchronization
        dthetadt (i) = omega(i) + 1/N*coup (i);
        theta (i,t+1) = mod(theta(i,t) + dt*(dthetadt(i)), 2* pi) ;
    end
end

%% order parameter
r = abs( mean ( exp(j*(theta(1:N,:))) ,1));



%% plots
figure (1), plot (r) 
xlabel('Time $t$','interpreter','latex','fontsize',20)
ylabel('order parameter $r$','interpreter','latex','fontsize',20)


figure(2), for i = 1:N
           plot(sin(theta(i,:))); hold on;
           xlabel('Time $t$','interpreter','latex','fontsize',20')
           ylabel('$\sin(\theta_{i})$','interpreter','latex','fontsize',20)
           end



