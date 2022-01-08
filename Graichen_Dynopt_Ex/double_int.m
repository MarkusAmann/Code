function [t,s,v,a,p] = double_int(h,f)
% ----------------------------------------------------
% h: gewunschte Enddistanz
% f: Gewichtungsfaktor fuer Bestrafung von a
% (t,s,v,a): Trajektorien des Partikels
% p: Parameterstruktur
p.f = f; % Parameter
p.s0=0; p.v0=5; % Anfangsbedingungen
p.sf=h; % Endbedingungen (Ubergabe aus Funktionsaufruf)
opt = optimset('Display','iter'); % Optionen
X0 = [-0.1,1,2]; % Startwert für [l1,a0,tf]
Xopt = fsolve(@eqns,X0,opt,p); % Numerische Losung mit fsolve
p.l1=Xopt(1); p.a0=Xopt(2); p.tf=Xopt(3); % Losung
t = linspace(0,p.tf,100); % Trajektorien
s = sfct(p.l1,p.a0,t,p);
v = vfct(p.l1,p.a0,t,p);
a = afct(p.l1,p.a0,t,p);
p.l2 = -a*p.f;

% ----------------------------------------------------
function res = eqns(X,p) % Gleichungen in Residuenform
l1=X(1); a0=X(2); tf=X(3);
res = [sfct(l1,a0,tf,p) - p.sf;
l1*vfct(l1,a0,tf,p) + 1;
-p.f*afct(l1,a0,tf,p)];
% ----------------------------------------------------
function s = sfct(l1,a0,t,p) % Funktionen fur s, v und a
s = 1/(6*p.f)*l1*t.^3 + 1/2*a0*t.^2 + p.v0*t + p.s0; 
function v = vfct(l1,a0,t,p)
v = 1/(2*p.f)*l1*t.^2 + a0*t + p.v0;
function a = afct(l1,a0,t,p)
a = l1/p.f*t + a0;

