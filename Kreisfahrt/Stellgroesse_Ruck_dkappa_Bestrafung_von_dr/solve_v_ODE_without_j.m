t0 = 0; tf = 1; N = 100; fax = 1; fay = 1; fr = 1; kapparef = 0.01; sf = 1000; drf = 0; psirf = 0;
x0 = [0 5 0 0 0 kapparef].'; l0 = [0 0 0 0 0 0].'; %l0 = 0.1*randn(4,1);
p.fax = fax; p.fay = fay; p.fr = fr; p.kapparef = kapparef; p.sf = sf; p.drf = drf; p.psirf = psirf;
p.x0 = x0; p.l0 = l0; p.t0 = t0; p.tf = tf; p.N = N;  

%%
p.a0 = 0.766109705407045*1.02; % v0 = 5
% p.a0 = 1.403740933445438; % v0 = 0.1
[t_1,ode_sol_1] = ode45(@(t,X) v_Sol_ODE_1(t,X,p),[0 30],[p.x0(2) p.a0]);
[t_2,ode_sol_2] = ode45(@(t,X) v_Sol_ODE_2(t,X,p),[0 80],p.x0(2));

v_RL = (2/(3*p.fay*p.kapparef^2))^(1/4);
l1_RL = -2*p.fay*p.kapparef^2*v_RL^3;
c = p.a0^2 - p.fay/p.fax*p.kapparef^2*p.x0(2)^4 - 2*l1_RL/p.fax*p.x0(2);
ax_sol = sqrt(p.fay/p.fax*p.kapparef^2*real(ode_sol_2).^4 + 2*l1_RL/p.fax*real(ode_sol_2) + c);

figure(1)
subplot(2,1,1)
plot(t_2,real(ode_sol_2),'-.','Linewidth',2)
ylabel('$v\,[m/s]$','Interpreter','Latex')
xlabel('$t\,[s]$','Interpreter','Latex')
grid on
hold on
subplot(2,1,2)
plot(t_2,ax_sol,'-.','Linewidth',2)
ylabel('$a_x\,[m/s^2]$','Interpreter','Latex')
xlabel('$t\,[s]$','Interpreter','Latex')
grid on
hold on

%%
function dXdt = v_Sol_ODE_1(t,X,p)
v = X(1);
dv = X(2);

v_RL = (2/(3*p.fay*p.kapparef^2))^(1/4);
l1_RL = -2*p.fay*p.kapparef^2*v_RL^3;

d2v = 2*p.fay/p.fax*p.kapparef^2*v^3 + l1_RL/p.fax;

dXdt = [dv; d2v];
end

%%
function dXdt = v_Sol_ODE_2(t,X,p)
v = X;

v_RL = (2/(3*p.fay*p.kapparef^2))^(1/4);
l1_RL = -2*p.fay*p.kapparef^2*v_RL^3;
c = p.a0^2 - p.fay/p.fax*p.kapparef^2*p.x0(2)^4 - 2*l1_RL/p.fax*p.x0(2);

dv = sqrt(p.fay/p.fax*p.kapparef^2*v^4 + 2*l1_RL/p.fax*v + c);

dXdt = dv;
end