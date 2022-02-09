clear all
close all
clc
[vec,sol,p] = schiff_rwa(-15,15)
u_deg = vec.u*180/pi;
u_opt_deg = -1/p.r*vec.adj'*[0; p.c(2); 0]*180/pi;
figure 
plot(vec.t,u_deg,vec.t,u_opt_deg)

function [vec,sol,p] = schiff_rwa(umin,umax)
% ----------------------------------------------------------------
p = schiff_param; % Schiffparameter
p.umin = umin*pi/180; % Stellgroenbeschrankungen
p.umax = umax*pi/180; % (Ubergabe in Grad)
opt = bvpset('Stats','on','RelTol',1e-4); % Optionen fur bvp4c
sol0 = bvpinit(linspace(0,p.tf,30),[p.x0;zeros(3,1)]); % Startlosung
sol = bvp4c(@schiff_kangln, @schiff_rb, sol0, opt, p); % Aufruf von bvp4c
vec.t = linspace(0,p.tf,100); % optimale Losung
vec.x = interp1(sol.x,sol.y(1:3,:)',vec.t)';
vec.adj = interp1(sol.x,sol.y(4:6,:)',vec.t)';
for i=1:length(vec.t)
vec.u(:,i) = schiff_uopt(vec.x(:,i),vec.adj(:,i),p); % Steuerung
end
end

% ----------------------------------------------------------------
function res = schiff_rb(X0,Xf,p) % Randbedingungen
x0 = X0(1:3); xf = Xf(1:3); adjf = Xf(4:6);
res = [ x0 - p.x0;
adjf - p.S*(xf-p.xf) ];
end

function p = schiff_param
% ----------------------------------------------------------------
p.c = [-0.26, 0.2, -1.87, 0.6]; % Schiffparameter
p.v = 3.5; % Geschwindigkeit
p.tf = 15; % Endzeit (Optimierungshorizont)
p.x0 = [0 ; 0; 0]; % Anfangsbedingungen
p.xf = [45*pi/180; 0; 0]; % (gewunschte) Endbedingungen
p.r = 0.1; % Gewichtungen im Integralanteil
p.Q = diag([1,1,1]);
p.S = diag([1,1,1]); % Gewichtung des Endzustands
end

function f = schiff_kangln(t,X,p)
% ----------------------------------------------------------------
x = X(1:3); adj = X(4:6);
u = schiff_uopt(x,adj,p);
f = [ fsys(x,u,p); % kanonische Gleichungen
-p.Q*(x-p.xf)-dfdx(x,adj,u,p)'*adj ];
end

function f = fsys(x,u,p)
% ----------------------------------------------------------------
f = [ x(2); % Schiffdynamik
p.c(1)*x(2)+p.c(2)*u;
p.c(3)*p.v*abs(x(3))*x(3)+p.c(4)*x(2) ];
end

function J = dfdx(x,adj,u,p)
% ----------------------------------------------------------------
J = [ 0, 1, 0; % Jacobimatrix
0, p.c(1), 0;
0, p.c(4), 2*p.c(3)*p.v*x(3)*sign(x(3)) ];
end

function u = schiff_uopt(x,adj,p)
% ----------------------------------------------------------------
u0 = -1/p.r*adj'*[0; p.c(2); 0]; % unbeschrankte Stellgroe
if u0> p.umin & u0<p.umax, u = u0;
elseif u0<=p.umin, u = p.umin;
else u = p.umax;
end
end