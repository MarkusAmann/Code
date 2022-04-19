clear all
close all
% clc

% FUNKTIONIERT NICHT!!!
%% hier parameter vorgeben
fa = 1;            % Gewichtung Längsbeschleunigung
fj = 2;            % Gewichtung Ruck
fu = 1;            % Gewichtung zeitl. Ableitung Ruck

tf = 40;
sf = 200;   % Länge der Gerade  

s0 = 0; v0 = 10; vf = 0; a0 = 0; af = 0; j0 = 0; jf = 0;

[t,s,v,a,l1,l2,p] = double_int(s0,v0,a0,j0,sf,vf,af,jf,fa,fj,fu,tf);
figure
subplot(3,1,1)
plot(t,s,'--','LineWidth',2)
ylabel('s [m]')
grid on
hold on
subplot(3,1,2)
plot(t,v,'--','LineWidth',2)
ylabel('v [m/s]')
grid on
hold on
subplot(3,1,3)
plot(t,a,'--','LineWidth',2)
ylabel('a [m/s^2]')
xlabel('t [s]')
grid on
hold on

figure('Name','adj')
subplot(2,1,1)
plot(t,l1,'--','LineWidth',2)
ylabel('l_{1,opt}')
grid on
hold on
subplot(2,1,2)
plot(t,l2,'--','LineWidth',2)
ylabel('l_{2,opt}')
grid on
hold on

J_fun = 1/2*p.f*a.^2;
J = trapz(t,J_fun) + p.tf;
fprintf('l1: %f\ntf: %f\nJ: %f\n',p.l1,p.tf,J)

function [t,s,v,a,l1,l2,p] = double_int(s0,v0,a0,j0,sf,vf,af,jf,fa,fj,fu,tf)
% ----------------------------------------------------
% p: Parameterstruktur
p.fa=fa; p.fj=fj; p.fu=fu; % Parameter
p.s0=0; p.v0=v0; p.a0=a0; p.j0=j0; % Anfangsbedingungen
p.tf=tf; p.sf=sf; p.vf=vf; p.af=af; p.jf=jf; % Endbedingungen (Ubergabe aus Funktionsaufruf)
opt = optimset('Display','iter','MaxFunEvals',1e5,'Algorithm','levenberg-marquardt'); % Optionen
X0 = [0,0,0,0,0,0,0,0]; % Startwert 
Xopt = fsolve(@eqns,X0,opt,p); % Numerische Losung mit fsolve
end
% ----------------------------------------------------
function res = eqns(X,p) % Gleichungen in Residuenform
c1=X(1); c2=X(2); c3=X(3); cs=X(4); C1=X(5); C2=X(6); C3=X(7); C4=X(8); 
res = [sfct(X,0,p) - p.s0;...
    vfct(X,0,p) - p.v0;...
    afct(X,0,p) - p.a0;...
    jfct(X,0,p) - p.j0;...
    sfct(X,p.tf,p) - p.sf;...
%     v(tf) == vf;...
%     a(tf) == af;...
%     j(tf) == jf;...
    l2fct(X,p.tf,p) - 0;...
    l3fct(X,p.tf,p) - 0;...
    djfct(X,p.tf,p) - 0
    ];
end
% ----------------------------------------------------
function s = sfct(X,t,p) % Funktionen fur s, v und a
c1=X(1); c2=X(2); c3=X(3); cs=X(4); C1=X(5); C2=X(6); C3=X(7); C4=X(8); 
s = cs + (p.fa*((c1*t^3)/6 - (c2*t^2)/2 + c3*t) + c1*p.fj*t)/p.fa^2 - (2^(1/2)*C1*exp(-t*(p.fj/(2*p.fu) + (p.fj^2 - 4*p.fa*p.fu)^(1/2)/(2*p.fu))^(1/2)))/(p.fj/p.fu + (p.fj^2 - 4*p.fa*p.fu)^(1/2)/p.fu)^(1/2) + (2^(1/2)*C2*exp(t*(p.fj/(2*p.fu) - (p.fj^2 - 4*p.fa*p.fu)^(1/2)/(2*p.fu))^(1/2)))/(p.fj/p.fu - (p.fj^2 - 4*p.fa*p.fu)^(1/2)/p.fu)^(1/2) + (2^(1/2)*C4*exp(t*(p.fj/(2*p.fu) + (p.fj^2 - 4*p.fa*p.fu)^(1/2)/(2*p.fu))^(1/2)))/(p.fj/p.fu + (p.fj^2 - 4*p.fa*p.fu)^(1/2)/p.fu)^(1/2) - (2^(1/2)*C3*exp(-t*(p.fj/(2*p.fu) - (p.fj^2 - 4*p.fa*p.fu)^(1/2)/(2*p.fu))^(1/2)))/(p.fj/p.fu - (p.fj^2 - 4*p.fa*p.fu)^(1/2)/p.fu)^(1/2); 
end
function v = vfct(X,t,p)
c1=X(1); c2=X(2); c3=X(3); cs=X(4); C1=X(5); C2=X(6); C3=X(7); C4=X(8); 
v = C2*exp(t*((p.fj - (p.fj^2 - 4*p.fa*p.fu)^(1/2))/(2*p.fu))^(1/2)) + C3*exp(-t*((p.fj - (p.fj^2 - 4*p.fa*p.fu)^(1/2))/(2*p.fu))^(1/2)) + (c3*p.fa + c1*p.fj)/p.fa^2 + C1*exp(-t*((p.fj + (p.fj^2 - 4*p.fa*p.fu)^(1/2))/(2*p.fu))^(1/2)) + C4*exp(t*((p.fj + (p.fj^2 - 4*p.fa*p.fu)^(1/2))/(2*p.fu))^(1/2)) + (c1*t^2)/(2*p.fa) - (c2*t)/p.fa;
end
function a = afct(X,t,p)
c1=X(1); c2=X(2); c3=X(3); cs=X(4); C1=X(5); C2=X(6); C3=X(7); C4=X(8); 
a = C2*exp(t*((p.fj - (p.fj^2 - 4*p.fa*p.fu)^(1/2))/(2*p.fu))^(1/2))*((p.fj - (p.fj^2 - 4*p.fa*p.fu)^(1/2))/(2*p.fu))^(1/2) - c2/p.fa - C3*exp(-t*((p.fj - (p.fj^2 - 4*p.fa*p.fu)^(1/2))/(2*p.fu))^(1/2))*((p.fj - (p.fj^2 - 4*p.fa*p.fu)^(1/2))/(2*p.fu))^(1/2) - C1*exp(-t*((p.fj + (p.fj^2 - 4*p.fa*p.fu)^(1/2))/(2*p.fu))^(1/2))*((p.fj + (p.fj^2 - 4*p.fa*p.fu)^(1/2))/(2*p.fu))^(1/2) + C4*exp(t*((p.fj + (p.fj^2 - 4*p.fa*p.fu)^(1/2))/(2*p.fu))^(1/2))*((p.fj + (p.fj^2 - 4*p.fa*p.fu)^(1/2))/(2*p.fu))^(1/2) + (c1*t)/p.fa;
end
function j = jfct(X,t,p)
c1=X(1); c2=X(2); c3=X(3); cs=X(4); C1=X(5); C2=X(6); C3=X(7); C4=X(8); 
j = c1/p.fa + (C1*exp(-t*((p.fj + (p.fj^2 - 4*p.fa*p.fu)^(1/2))/(2*p.fu))^(1/2))*(p.fj + (p.fj^2 - 4*p.fa*p.fu)^(1/2)))/(2*p.fu) + (C4*exp(t*((p.fj + (p.fj^2 - 4*p.fa*p.fu)^(1/2))/(2*p.fu))^(1/2))*(p.fj + (p.fj^2 - 4*p.fa*p.fu)^(1/2)))/(2*p.fu) + (C2*exp(t*((p.fj - (p.fj^2 - 4*p.fa*p.fu)^(1/2))/(2*p.fu))^(1/2))*(p.fj - (p.fj^2 - 4*p.fa*p.fu)^(1/2)))/(2*p.fu) + (C3*exp(-t*((p.fj - (p.fj^2 - 4*p.fa*p.fu)^(1/2))/(2*p.fu))^(1/2))*(p.fj - (p.fj^2 - 4*p.fa*p.fu)^(1/2)))/(2*p.fu);
end
function dj = djfct(X,t,p)
c1=X(1); c2=X(2); c3=X(3); cs=X(4); C1=X(5); C2=X(6); C3=X(7); C4=X(8); 
dj = (C4*exp(t*((p.fj + (p.fj^2 - 4*p.fa*p.fu)^(1/2))/(2*p.fu))^(1/2))*(p.fj + (p.fj^2 - 4*p.fa*p.fu)^(1/2))*((p.fj + (p.fj^2 - 4*p.fa*p.fu)^(1/2))/(2*p.fu))^(1/2))/(2*p.fu) - (C1*exp(-t*((p.fj + (p.fj^2 - 4*p.fa*p.fu)^(1/2))/(2*p.fu))^(1/2))*(p.fj + (p.fj^2 - 4*p.fa*p.fu)^(1/2))*((p.fj + (p.fj^2 - 4*p.fa*p.fu)^(1/2))/(2*p.fu))^(1/2))/(2*p.fu) + (C2*exp(t*((p.fj - (p.fj^2 - 4*p.fa*p.fu)^(1/2))/(2*p.fu))^(1/2))*(p.fj - (p.fj^2 - 4*p.fa*p.fu)^(1/2))*((p.fj - (p.fj^2 - 4*p.fa*p.fu)^(1/2))/(2*p.fu))^(1/2))/(2*p.fu) - (C3*exp(-t*((p.fj - (p.fj^2 - 4*p.fa*p.fu)^(1/2))/(2*p.fu))^(1/2))*(p.fj - (p.fj^2 - 4*p.fa*p.fu)^(1/2))*((p.fj - (p.fj^2 - 4*p.fa*p.fu)^(1/2))/(2*p.fu))^(1/2))/(2*p.fu);
end
function l1 = l1fct(X,t,p)
c1=X(1); c2=X(2); c3=X(3); cs=X(4); C1=X(5); C2=X(6); C3=X(7); C4=X(8); 
l1 = c1;
end
function l2 = l2fct(X,t,p)
c1=X(1); c2=X(2); c3=X(3); cs=X(4); C1=X(5); C2=X(6); C3=X(7); C4=X(8); 
l2 = -c1.*t + c2;
end
function l3 = l3fct(X,t,p)
c1=X(1); c2=X(2); c3=X(3); cs=X(4); C1=X(5); C2=X(6); C3=X(7); C4=X(8); 
l3 = -p.fa*vfct(X,t,p) + l1fct(X,t,p)/2*t^2 - c2*t + c3;
end
