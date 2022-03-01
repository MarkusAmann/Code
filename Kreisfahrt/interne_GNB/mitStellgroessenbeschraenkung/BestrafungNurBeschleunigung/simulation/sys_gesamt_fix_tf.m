function f = sys_gesamt_fix_tf(t,X,region,p)
% delta_t1 = param;
x = X(1:4); 
u = uopt(X,p);
% switch region
%     case 1
%         f = delta_t1*[fsys(x,u,p);...% kanonische Gleichungen
%             kangln(X,u,p)];
%     case 2 
%         f = [fsys(x,u,p);...% kanonische Gleichungen
%             kangln(X,u,p)];
% end
f = [fsys(x,u,p);...% kanonische Gleichungen
    kangln(X,u,p)];
end