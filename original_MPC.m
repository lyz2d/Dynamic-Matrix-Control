%% initialization
n=50;
p=10;
m=1;
weight=0.0;
ysp=1;
timesp=1;
delt=0.1;
tfinal=6;
noise=0;

t=0:delt:tfinal;
kfinal=length(t);

%% System Formulation
ksp=fix(timesp/delt);
r=[zeros(ksp,1);ones(kfinal-ksp,1)*ysp];
% a = [0 1 0 0;
%     51.2 1 2.56 0;
%     0 0 0 1;
%     128.02 0 -6.4 10.2];
% b = [0;0;0;1];
% c1 = [1 0 0 0];
% d = [0];
% c2 = [1280.2 0 64.01 0];
% prompt1 = "What is the matrix A?\n";
% A = input(prompt1);

sys_mo = LTISystem();
% [A,B,C,D] = ssdata(sys_mod);
% sys_mod1=ss(a,b,c1,d);
% sys_mod2=ss(a,b,c2,d);
% model1 = c2d(sys_mod1,delt);
% model2 = c2d(sys_mod2,delt);
for i = size(D,1)
    [A,B,C{i},D{i}]=ssdata(sys_model{i});
end
% [A,B,C1,D]=ssdata(model1);
% [A,B,C2,D]=ssdata(model2);

%% Step Response
%1,产生阶跃响应
for i = size(D,1)
    [s{i}]=step(sys_mod{i},delt:delt:n*delt);
    plot(s{i},'o') 
    pause
end
% [s2]=step(sys_mod2,[delt:delt:n*delt]);
% plot(s2,'-')
% pause

%% DMC
%2，计算预测部分系数和已知部分系数
% for i = size(D,1)
%     [Sf{i},Sp{i},Kmat{i}]=original_smatgen(s{i},p,m,n,weight);
% end
[Sf1,Sp1,Kmat1]=original_smatgen(s1,p,m,n,weight);
[Sf2,Sp2,Kmat2]=original_smatgen(s2,p,m,n,weight);

%3，其他初始化设置
xinit=zeros(size(a,1),1);
uinit=0;
yinit=0;

u = ones(min(p,kfinal),1)*uinit;
dup = zeros(n-2,1);
sn1 = s1(n);
sn2 = s2(n);
x(:,1) = xinit;
y1(1) = yinit;
y2(1) = yinit;
dist(1) = 0;

%%initial
%start simulation
for k=1:kfinal
    du(k)=original_dmccalc(Sp1,Kmat1,sn1,dup,dist(k),r(k),u,k,n);
    if k>1
        u(k)=u(k-1)+du(k);
    else
        u(k)=uinit+du(k);
    end
    x(:,k+1)=A*x(:,k)+B*u(k);
    y1(k+1)=C1*x(:,k+1);
    if(k-n+1)>0
        ymod(k+1)=s1(1)*du(k)+Sp1(1,:)*dup+sn1*u(k-n+1);
    else
        ymod(k+1)=s1(1)*du(k)+Sp1(1,:)*dup;
    end
    dist(k+1)=y1(k+1)-ymod(k+1);
    
    dup=[du(k);dup(1:(n-3))];
end

for k=1:kfinal
    du(k)=original_dmccalc(Sp2,Kmat2,sn2,dup,dist(k),r(k),u,k,n);
    if k>1
        u(k)=u(k-1)+du(k);
    else
        u(k)=uinit+du(k);
    end
    x(:,k+1)=A*x(:,k)+B*u(k);
    y2(k+1)=C2*x(:,k+1);
    if(k-n+1)>0
        ymod(k+1)=s2(1)*du(k)+Sp2(1,:)*dup+sn2*u(k-n+1);
    else
        ymod(k+1)=s2(1)*du(k)+Sp2(1,:)*dup;
    end
    dist(k+1)=y2(k+1)-ymod(k+1);
    
    dup=[du(k);dup(1:(n-3))];
end

[tt,uu]=stairs(t,u);
[ttr,rr]=stairs(t,r);
figure(1)
subplot(2,1,1)
plot(ttr,rr,'o',t,y1(1:length(t)),'o')
ylabel('y')
xlabel('time')
title('time plant')
subplot(2,1,2)
plot(tt,uu,'o')
ylabel('u')
xlabel('time')
pause
[tt,uu]=stairs(t,u);
[ttr,rr]=stairs(t,r);
figure(1)
subplot(2,1,1)
plot(ttr,rr,'-',t,y2(1:length(t)))
ylabel('y')
xlabel('time')
title('time plant')
subplot(2,1,2)
plot(tt,uu)
ylabel('u')
xlabel('time')