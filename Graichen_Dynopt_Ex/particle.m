function [t,x,y,p] = particle(xf,yf)
% ----------------------------------------------------
% (xf,yf): gewunschter Endpunkt
% (t,x,y): Trajektorien des Partikels
% p: Parameterstruktur
p.a = 1; % Parameter
p.x0=0; p.v0=0; p.y0=0; p.w0=1; % Anfangsbedingungen
p.xf=xf; p.yf=yf; % Endbedingungen (Ubergabe aus Funktionsaufruf)
opt = optimset('Display','iter'); % Optionen
X0 = [-1,0,1]; % Startwert
Xopt = fsolve(@eqns,X0,opt,p); % Numerische Losung mit fsolve
p.lx=Xopt(1); p.ly=Xopt(2); p.tf=Xopt(3); % Losung
t = linspace(0,p.tf,100); % Trajektorien
x = xfct(p.lx,p.ly,t,p);
y = yfct(p.lx,p.ly,t,p);
% ----------------------------------------------------
function res = eqns(X,p) % Gleichungen in Residuenform
lx=X(1); ly=X(2); tf=X(3);
res = [ xfct(lx,ly,tf,p) - p.xf;
yfct(lx,ly,tf,p) - p.yf;
lx*vfct(lx,ly,tf,p) + ly*wfct(lx,ly,tf,p) + 1 ];
% ----------------------------------------------------
function x = xfct(lx,ly,t,p) % Funktionen fur x und v
cosu = 1/sqrt(1+(ly/lx)^2);
x = p.x0 + p.v0*t + p.a/2*cosu*t.^2; % 't.^2' steht fur komponentenweise Auswertung
function v = vfct(lx,ly,t,p)
cosu = 1/sqrt(1+(ly/lx)^2);
v = p.v0 + p.a*cosu*t;
% ----------------------------------------------------
function y = yfct(lx,ly,t,p) % Funktionen fur y und w
sinu = ly/(lx*sqrt(1+(ly/lx)^2));
y = p.y0 + p.w0*t + p.a/2*sinu*t.^2; % 't.^2' steht fur komponentenweise Auswertung
function w = wfct(lx,ly,t,p)
sinu = ly/(lx*sqrt(1+(ly/lx)^2));
w = p.w0 + p.a*sinu*t;

