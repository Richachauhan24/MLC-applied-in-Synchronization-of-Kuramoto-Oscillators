function J=unidim_DS_evaluator(ind,mlc_parameters,i,fig)
%% Obtaining parameters from MLC object.
Tf=mlc_parameters.problem_variables.Tf;
objective=mlc_parameters.problem_variables.objective;
gamma=mlc_parameters.problem_variables.gamma;
Tevmax=mlc_parameters.problem_variables.Tevmax;

N=mlc_parameters.problem_variables.N;
A=mlc_parameters.problem_variables.A;
omega=mlc_parameters.problem_variables.omega;
kappa=mlc_parameters.problem_variables.kappa;
B=mlc_parameters.problem_variables.B;
 
%% Interpret individual.
m=readmylisp_to_formal_MLC(ind);
m=strrep(m,'S0','y');
K=@(y)(1/N*sum(B.*sin(y-y'))'); % K=@(y)(1/N*sum(A.*sin(y-y'))'); when matrix A is not weighted then matrix B is not required.
eval(['K=@(y)(' m ');']);
f=@(t,y) (omega + kappa/N*sum(A.*sin(y-y'))' + K(y) + testt(toc,Tevmax));
   

%% Evaluation
try                       % Encapsulation in try/catch.
tic 
[T,Y]=ode45(f,[0 Tf],rand(N,1)*2*pi);  % Integration.
if T(end)==Tf             % Check if Tf is reached.
    b=Y*0+K(Y);           % Computes b.
    Jt=1/Tf*cumtrapz(T,(abs( mean ( exp(j*(Y)),2))-objective).^2 + gamma*sum(b.^2,2)); % Computes J. Mean r value should be one, objective = 1
    J=Jt(end);
else
    J=mlc_parameters.badvalue;  % Return high value if Tf is not reached.
end
catch err
   %J=mlc_parameters.badvalue % Return high value if ode45 fails.
   J = 1e+33
end
    
if nargin>3   % If a fourth argument is provided, plot the result
    subplot(3,1,1)
    plot(T,Y,'-k','linewidth',1.2)
    ylabel('$\theta_{i}$','interpreter','latex','fontsize',20)
    subplot(3,1,2)
    plot(T,b,'-k','linewidth',1.2)
    ylabel('$u$','interpreter','latex','fontsize',20)
    subplot(3,1,3)
    plot(T,Jt,'-k','linewidth',1.2)
    ylabel('$(r-1)^2+\gamma u^2$','interpreter','latex','fontsize',20)
    xlabel('Time $t$','interpreter','latex','fontsize',20)


    figure (2), plot (abs( mean ( exp(j*(Y)),2))) 
    xlabel('Time $t$','interpreter','latex','fontsize',20)
    ylabel('order parameter $r$','interpreter','latex','fontsize',20)

    figure(3), for i = 1:N
               plot(sin(Y)); hold on;
               end
               xlabel('Time $t$','interpreter','latex','fontsize',20)
               ylabel('$\sin(\theta_{i})$','interpreter','latex','fontsize',20)
end
