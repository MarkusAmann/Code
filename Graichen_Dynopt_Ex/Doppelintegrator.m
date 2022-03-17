clear all
close all
clc

h = 100; f = 1;
[t,s,v,a,p] = double_int(h,f);
figure
subplot(3,1,1)
plot(t,s)
ylabel('s [m]')
grid on
hold on
subplot(3,1,2)
plot(t,v)
ylabel('v [m/s]')
grid on
hold on
subplot(3,1,3)
plot(t,a)
ylabel('a [m/s^2]')
xlabel('t [s]')
grid on
hold on

p.l1*vfct(p.l1,p.a0,p.tf,p) + 1
p.l1
p.tf

roots = fsolve(@(X) funs(X,p),[-0.1, 11]);

function res = funs(X,p)
l1 = X(1); tf = X(2);
a0 = -l1*tf/p.f;
s_tf = 1/(6*p.f)*l1*tf^3 + 1/2*a0*tf^2 + p.v0*tf + p.s0;
v_tf = l1*tf^2/(2*p.f) + a0*tf + p.v0;
res = [s_tf - p.sf;...
    v_tf*l1 + 1];
end 
% syms  tf
% l1 = -3/2*p.f*(p.sf - p.s0 - p.v0*tf)/tf^3;
% a0 = -l1*tf/p.f;
% v_tf = l1*tf^2/(2*p.f) - l1*tf/p.f*tf + p.v0;
% tf_ana = solve(v_tf*l1 + 1 == 0, tf);

function [t,s,v,a,p] = double_int(h,f)
% ----------------------------------------------------
% h: gewunschte Enddistanz
% f: Gewichtungsfaktor fuer Bestrafung von a
% (t,s,v,a): Trajektorien des Partikels
% p: Parameterstruktur
p.f = f; % Parameter
p.s0=0; p.v0=5; % Anfangsbedingungen
p.sf=h; % Endbedingungen (Ubergabe aus Funktionsaufruf)
opt = optimset('Display','iter','MaxFunEvals',1e3); % Optionen
X0 = [-0.1,1,2]; % Startwert für [l1,a0,tf]
Xopt = fsolve(@eqns,X0,opt,p); % Numerische Losung mit fsolve
p.l1=Xopt(1); p.a0=Xopt(2); p.tf=Xopt(3); % Losung
t = linspace(0,p.tf,100); % Trajektorien
s = sfct(p.l1,p.a0,t,p);
v = vfct(p.l1,p.a0,t,p);
a = afct(p.l1,p.a0,t,p);
p.l2 = -a*p.f;
end
% ----------------------------------------------------
function res = eqns(X,p) % Gleichungen in Residuenform
l1=X(1); a0=X(2); tf=X(3);
res = [sfct(l1,a0,tf,p) - p.sf;
l1*vfct(l1,a0,tf,p) + 1;
-p.f*afct(l1,a0,tf,p)];
end
% ----------------------------------------------------
function s = sfct(l1,a0,t,p) % Funktionen fur s, v und a
s = 1/(6*p.f)*l1*t.^3 + 1/2*a0*t.^2 + p.v0*t + p.s0; 
end
function v = vfct(l1,a0,t,p)
v = 1/(2*p.f)*l1*t.^2 + a0*t + p.v0;
end
function a = afct(l1,a0,t,p)
a = l1/p.f*t + a0;
end
