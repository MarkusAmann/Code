syms sr v a dr psir kappa l1 l2 l3 l4 l5 l6 j dkappa...
    sr_RL_sym v_RL_sym a_RL_sym dr_RL_sym psir_RL_sym kappa_RL_sym l1_RL_sym l2_RL_sym l3_RL_sym l4_RL_sym l5_RL_sym l6_RL_sym ...
    fax_sym fjx_sym fay_sym fjy_sym fr_sym kappar_sym

t0 = 0; tf = 1; N = 100; fjy = 1; fjx = 1; fax = 1; fay = 1; fr = 1; kapparef = 0.01; sf = 1000; drf = 0; psirf = 0;
x0 = [0 5 0 0 0 kapparef].'; l0 = [0 0 0 0 0 0].'; %l0 = 0.1*randn(4,1);
p.fjx = fjx; p.fjy = fjy; p.fax = fax; p.fay = fay; p.fr = fr; p.kapparef = kapparef; p.sf = sf; p.drf = drf; p.psirf = psirf;
p.x0 = x0; p.l0 = l0; p.t0 = t0; p.tf = tf; p.N = N;  

sr_RL = 10;
v_RL = ((3*p.fr)/(2*p.fay*p.kapparef^4)-sqrt(((3*p.fr)/(2*p.fay*p.kapparef^4))^2-((2*p.fr)/(p.fay^2*p.kapparef^6))))^(1/4);
a_RL = 0;
dr_RL = p.kapparef^3*v_RL^4*p.fay/p.fr;
psir_RL = 0;
kappa_RL = p.kapparef;
l1_RL = -2*p.fay*p.kapparef^2*v_RL^3;
l2_RL = 0;
l3_RL = 0;
l4_RL = 0;
l5_RL = -p.fay*p.kapparef*v_RL^3;
l6_RL = 0;
z_RL_num = [sr_RL; v_RL; a_RL; dr_RL; psir_RL; kappa_RL];

z_RL_sym = [sr_RL_sym; v_RL_sym; a_RL_sym; dr_RL_sym; psir_RL_sym; kappa_RL_sym];
z = [sr v a dr psir kappa].';

% j = -l3/fjx_sym; 
% dkappa = -(l6 + 2*fjy_sym*kappa*a*v^3)/(fjy_sym*v^4);

f = [v*cos(psir)/(1-dr*kappar_sym);...
    a;...
    j;...
    v*sin(psir);...
    kappa*v - kappar_sym*v*cos(psir)/(1-dr*kappar_sym);...
    dkappa];

A = jacobian(f,z);
A_RL_sym = subs(A,z,z_RL_sym);
A_RL_num = double(subs(subs(A,z,z_RL_num),[fjx_sym fjy_sym fax_sym fay_sym fr_sym kappar_sym],[p.fjx p.fjy p.fax p.fay p.fr p.kapparef]));

ew = eig(A_RL_num);
ew_real = real(ew);
ew_imag = imag(ew);
figure
plot(ew_real,ew_imag,'*','MarkerSize',8)
grid on
hold on
set(gcf,'renderer','Painters')
xlabel('Real','Interpreter','Latex')
ylabel('Imag','Interpreter','Latex')