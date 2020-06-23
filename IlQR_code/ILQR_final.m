clc;
clear all;

mc=1    ;  % mass of cart
mp=0.3  ;  % mass of pendulum
l=0.5   ;  % length of pendulum
g=9.806 ;  % gravity
T= 501  ;


[x0,u0]=Cart_Pole_Energy_Shaping();
xp=x0;
up=u0;

Q_cost = [10 0 0 0;...
    0 10 0 0;...
    0 0 10 0;...
    0 0 0 10];
R_cost = 10;
del_U_all=[];
del_x_all= zeros(T,4);
K_all = zeros(T,4);
k_all = zeros(T,1);
V = zeros(4,4);
v = zeros(4,1);

% for j=1:5
for i=T:-1:1
    
    [C,c]=Cost(Q_cost,R_cost,i,xp(i,:),up(i));
    [A,B]=Linearize_dynamics(mc,mp,l,g,xp(i,:),up(i));
    F=[A B];
    
    f = state_with_input(mc,mp,l,g,xp(i,:),up(i));
    
    i
    Qt = C + F'*V*F;
    qt = c + F'*V*f + F'*v;
    
    Quu = Qt(5,5);
    Qux = Qt(5,1:4);
    Qxx = Qt(1:4,1:4);
    Qxu =  Qt(1:4,5);
    q_ut = qt(5,1);
    q_xt = qt(1:4,1);
    
    k = -1*(Quu)\q_ut;
    K = -1*(Quu)\Qux;
    V = Qxx + Qxu*K + K'*Qux+K'*Quu*K;
    v = q_xt + Qxu*k +K'*q_ut+K'*Quu*k;
    
    K_all(i,:) = K;
    k_all(i,:) = k;
    
end






for i=1:T
    i

   [A,B]=Linearize_dynamics(mc,mp,l,g,xp(i,:),up(i));
    del_u = K_all(i,:)*del_x_all(i,:)'+0.2*k_all(i,1);
    
    u0(i) = up(i) + del_u;
    
    xg=xp(i+1,:);
    x0(i+1,:)=A*x0(i,:)'+B*u0(i);
    x0(i+1,2)=wrapToPi(x0(i+1,2));
    
    del_x_all(i+1,:)=x0(i+1,:)- xp(i+1,:);
   
   
%     del_x = A*del_x_all(i,:)'+B*del_u;
%     del_x_all(i+1,:) = del_x';
%     
%     x0(i+1,:) =x0(i+1,:) + del_x_all(i+1,:);
    
end
xp = x0;
up = u0;

% end