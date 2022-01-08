function [t,s,v,a,p] = kreisfahrt(sf,kapparef,fx,fy,varargin)
% ----------------------------------------------------
% sf: gewunschte Enddistanz
% fx: Gewichtungsfaktor fuer Bestrafung von ax
% fy: Gewichtungsfaktor fuer Bestrafung von ay
% varargin: vf gew unschte Endgeschw., falls nicht frei gelassen
% (t,s,v,a): Trajektorien des Partikels
% p: Parameterstruktur
if nargin == 5
    vf = varargin;
    p.vf = vf;
end
p.fx = fx; p.fy = fy; % Parameter
p.kapparef = kapparef; % konst. Krümmung
p.s0=0; p.v0=0; % Anfangsbedingungen
p.sf=sf; % Endbedingungen (Ubergabe aus Funktionsaufruf)
opt = optimset('Display','iter'); % Optionen
X0 = [-0.1,0.5,10]; % Startwert für [l1,a0,tf]
Xopt = fsolve(@eqns,X0,opt,p); % Numerische Losung mit fsolve
p.l1=Xopt(1); p.a0=Xopt(2); p.tf=Xopt(3); % Losung
t = linspace(0,p.tf,100); % Trajektorien
v_traj = v_eval(t,p);
s = sfct(p.l1,p.a0,t,p);
v = vfct(p.l1,p.a0,t,p);
a = afct(p.l1,p.a0,t,p);
p.l2 = -2*a*p.f;

% ----------------------------------------------------
function res = eqns(X,p) % Gleichungen in Residuenform
l1=X(1); a0=X(2); tf=X(3);
if isfield(p,'vf')
    res = [sfct(l1,a0,tf,p) - p.sf;
        vfct(l1,a0,tf,p) - p.vf;
        p.fy*vfct(l1,a0,tf,p)^4*p.kapparef^2 + l1*vfct(l1,a0,tf,p) + 1];
else
    res = [sfct(l1,a0,tf,p) - p.sf;
        afct(l1,a0,tf,p);
        p.fy*vfct(l1,a0,tf,p)^4*p.kapparef^2 + l1*vfct(l1,a0,tf,p) + 1];
end
% ----------------------------------------------------
function s = sfct(l1,a0,t,p) % Funktionen fur s, v und a
s = p.fy/(60*p.fx)*vfct(l1,a0,t,p).^6*p.kapparef^2 + l1/(12*p.fx).*t^3 + a0/2*t^2 + p.v0*t + p.s0; 
function vv = vfct(l1,a0,t,p)
persistent v_init
if isempty(v_init)
    v_init = 0;
end
% v = p.fy/(10*p.fx)*vfct(l1,a0,t,p).^5*p.kapparef^2 + l1/(4*p.fx).*t^2 + a0*t + p.v0;
fcn = @(v) p.fy/(10*p.fx)*v^5*p.kapparef^2 + l1/(4*p.fx).*t^2 + a0*t + p.v0 - v;
vv = fzero(fcn,v_init);
v_init = vv;
function a = afct(l1,a0,t,p)
a = p.fy/(2*p.fx)*vfct(l1,a0,t,p).^4*p.kapparef^2 + l1/(2*p.fx)*t + a0;
% -----------------------------------------------------
function v_traj = v_eval(t,p)
l1=p.l1; a0=p.a0; 
v_traj=[];
v_init = 0;
for i=1:length(t)
    ti = t(i);
    fcn = @(v) p.fy/(10*p.fx)*v^5*p.kapparef^2 + l1/(4*p.fx).*ti^2 + a0*ti + p.v0 - v;
    v_init = fzero(fcn,v_init);
    v_traj = [v_traj; v_init];
end
