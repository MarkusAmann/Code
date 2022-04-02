clear all
close all
% clc

h = 100; v0 = 10; f = 1; s1 = 20;
[t,s,v,a,l1,l2,p] = double_int(h,f,v0,s1);
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
J = trapz(t,J_fun) + p.tf + 1*p.t1;
fprintf('l1_1: %f\nl1_2: %f\nnu_tilde: %f\nt1: %f\ntf: %f\nJ: %f\n',p.l1_1,p.l1_2,p.nu_tilde,p.t1,p.tf,J)

function [t,s,v,a,l1,l2,p] = double_int(h,f,v0,s1)
% ----------------------------------------------------
% h: gewunschte Enddistanz
% f: Gewichtungsfaktor fuer Bestrafung von a
% (t,s,v,a): Trajektorien des Partikels
% p: Parameterstruktur
p.f = f; % Parameter
p.s0=0; p.v0=v0; p.s1 = s1;% Anfangsbedingungen
p.sf=h; % Endbedingungen (Ubergabe aus Funktionsaufruf)
opt = optimset('Algorithm','levenberg-marquardt','Display','iter','MaxFunEvals',1e4,'MaxIter',1e4); % Optionen
X0 = [-0.07,-0.07,0.65,10,50,10,2,0]; % Startwert für [l1_1,l1_2,a0,v_int_const,s_int_const,tf,t1,nu_tilde]
Xopt = fsolve(@eqns,X0,opt,p); % Numerische Losung mit fsolve
p.l1_1=Xopt(1); p.l1_2=Xopt(2); p.a0=Xopt(3); p.v_int_const=Xopt(4); p.s_int_const=Xopt(5); p.tf=Xopt(6); p.t1=Xopt(7); p.nu_tilde=Xopt(8);% Losung
t = linspace(0,p.tf,1000); % Trajektorien
s = sfct(p.l1_1,p.l1_2,p.tf,t,p.a0,p.v_int_const,p.s_int_const,p);
v = vfct(p.l1_1,p.l1_2,p.tf,t,p.a0,p.v_int_const,p);
a = afct(p.l1_1,p.l1_2,p.tf,t,p.a0,p);
l1 = l1fct(p.l1_1,p.l1_2,p.tf,t,p);
l2 = l2fct(p.l1_1,p.l1_2,p.tf,t,p.a0,p);
end
% ----------------------------------------------------
function res = eqns(X,p) % Gleichungen in Residuenform
l1_1=X(1);l1_2=X(2);a0=X(3);v_int_const=X(4);s_int_const=X(5);tf=X(6);t1=X(7);nu_tilde=X(8);
H_tf = l1fct(l1_1,l1_2,tf,tf,p,2)*vfct(l1_1,l1_2,tf,tf,a0,v_int_const,p,2); %+ 1/2*p.f*afct(l1_1,l1_2,tf,tf,a0,p,2)^2 + l2fct(l1_1,l1_2,tf,tf,a0,p,2)*afct(l1_1,l1_2,tf,tf,a0,p,2);
H1 = 1/2*p.f*afct(l1_1,l1_2,tf,t1,a0,p,1)^2 + l1fct(l1_1,l1_2,tf,t1,p,1)*vfct(l1_1,l1_2,tf,t1,a0,v_int_const,p,1) + l2fct(l1_1,l1_2,tf,t1,a0,p,1)*afct(l1_1,l1_2,tf,t1,a0,p,1);
H2 = 1/2*p.f*afct(l1_1,l1_2,tf,t1,a0,p,2)^2 + l1fct(l1_1,l1_2,tf,t1,p,2)*vfct(l1_1,l1_2,tf,t1,a0,v_int_const,p,2) + l2fct(l1_1,l1_2,tf,t1,a0,p,2)*afct(l1_1,l1_2,tf,t1,a0,p,2);

res = [sfct(l1_1,l1_2,tf,tf,a0,v_int_const,s_int_const,p,2) - p.sf;
H_tf + 1;
sfct(l1_1,l1_2,tf,t1,a0,v_int_const,s_int_const,p,1) - p.s1;
sfct(l1_1,l1_2,tf,t1,a0,v_int_const,s_int_const,p,1) - sfct(l1_1,l1_2,tf,t1,a0,v_int_const,s_int_const,p,2);
vfct(l1_1,l1_2,tf,t1,a0,v_int_const,p,1) - vfct(l1_1,l1_2,tf,t1,a0,v_int_const,p,2);
afct(l1_1,l1_2,tf,t1,a0,p,1) - afct(l1_1,l1_2,tf,t1,a0,p,2);
l1fct(l1_1,l1_2,tf,t1,p,1) - l1fct(l1_1,l1_2,tf,t1,p,2) - 2*nu_tilde;
H1 - H2 + 1];
end

% ----------------------------------------------------
function s = sfct(l1_1,l1_2,tf,t,a0,v_int_const,s_int_const,p,varargin) % Funktionen fur s, v und a
if nargin == 9
    l1orl2 = varargin{1,1};
    if l1orl2 == 1
        s = l1_1/(6*p.f).*t.^3 + a0/2*t.^2 + p.v0*t + p.s0;
    elseif l1orl2 == 2
        s = l1_2/(6*p.f).*t.^3 - l1_2/(2*p.f)*tf.*t.^2 + v_int_const*t + s_int_const;
    end
elseif nargin == 8
    t0t1 = t(t<=p.t1);
    t1t2 = t(t>p.t1);
    st0t1 = l1_1/(6*p.f).*t0t1.^3 + a0/2*t0t1.^2 + p.v0*t0t1 + p.s0;
    st1t2 = l1_2/(6*p.f).*t1t2.^3 - l1_2/(2*p.f)*tf.*t1t2.^2 + v_int_const*t1t2 + s_int_const;
    s = [st0t1 st1t2];
end
end
% ----------------------------------------------------
function v = vfct(l1_1,l1_2,tf,t,a0,v_int_const,p,varargin)
if nargin == 8
    l1orl2 = varargin{1,1};
    if l1orl2 == 1
        v = l1_1/(2*p.f).*t.^2 + a0*t + p.v0;
    elseif l1orl2 == 2
        v = l1_2/(2*p.f).*t.^2 - l1_2/p.f*tf.*t + v_int_const;
    end
elseif nargin == 7
    t0t1 = t(t<=p.t1);
    t1t2 = t(t>p.t1);
    vt0t1 = l1_1/(2*p.f).*t0t1.^2 + a0*t0t1 + p.v0;
    vt1t2 = l1_2/(2*p.f).*t1t2.^2 - l1_2/p.f*tf.*t1t2 + v_int_const;
    v = [vt0t1 vt1t2];
end
end
% ----------------------------------------------------
function a = afct(l1_1,l1_2,tf,t,a0,p,varargin)
if nargin == 7
    l1orl2 = varargin{1,1};
    if l1orl2 == 1
        a = l1_1/p.f.*t + a0;
    elseif l1orl2 == 2
        a = l1_2/p.f.*t - l1_2/p.f.*tf;
    end
elseif nargin == 6
    t0t1 = t(t<=p.t1);
    t1t2 = t(t>p.t1);
    at0t1 = l1_1/p.f.*t0t1 + a0;
    at1t2 = l1_2/p.f.*t1t2 - l1_2/p.f.*tf;
    a = [at0t1 at1t2];
end
end
% ----------------------------------------------------
function l1 = l1fct(l1_1,l1_2,tf,t,p,varargin)
if nargin ==6
    l1orl2 = varargin{1,1};
    if l1orl2 == 1
        l1 = l1_1;
    elseif l1orl2 == 2
        l1 = l1_2;
    end
elseif nargin == 5
    t0t1 = t(t<=p.t1);
    t1t2 = t(t>p.t1);
    l1t0t1 = l1_1*ones(size(t0t1));
    l1t1t2 = l1_2*ones(size(t1t2));
    l1 = [l1t0t1 l1t1t2];
end
end
% ----------------------------------------------------
function l2 = l2fct(l1_1,l1_2,tf,t,a0,p,varargin)
if nargin == 7
    l1orl2 = varargin{1,1};
    if l1orl2 == 1
        l2 = -l1_1.*t - a0*p.f;
    elseif l1orl2 == 2
        l2 = -l1_2.*t + l1_2.*tf;
    end
elseif nargin == 6
    t0t1 = t(t<=p.t1);
    t1t2 = t(t>p.t1);
    l2t0t1 = -l1_1.*t0t1 - a0*p.f;
    l2t1t2 = -l1_2.*t1t2 + l1_2.*tf;
    l2 = [l2t0t1 l2t1t2];
end
end
