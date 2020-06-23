 function [X_main,tou]= Cart_Pole_Energy_Shaping()



mc=1;  % mass of cart
mp=0.3;% mass of pendulum
l=0.5; % length of pendulum
g=9.81;
dt=0.01;
X_main=[0 0 0 0]';
tou=[];
i=1;

X=[1 0.03 0 0.02]';
X_dot=[0 0 0 0]';
X_des=[0, pi,0,0]';
% initializing variables
for t=0:dt:5
    x1=X(1);
    theta1=wrapToPi(X(2));
    x1_dot=X(3);
    theta1_dot=X(4);
    
    
    % defining dynamics of a system
    
    H=[mc+mp , mp*l*cos(theta1);...
        mp*l*cos(theta1), mp*(l^2)];
    
    C=[0 , -mp*l*theta1_dot*sin(theta1);...
        0, 0];
    
    G=[0, mp*g*l*sin(theta1)]';
    
    B=[1 0]';
    
    
    % Applying energy shaping
    
    E_des=mp*g*l;
    E=(0.5*mp*(l^2)*(theta1_dot^2))-mp*g*l*cos(theta1);
    E_ref=E-E_des;
    ke=2;
    kp=1;
    kd=0.8;
    x2dot_des=ke*(theta1_dot*cos(theta1)*E_ref)-kp*x1-kd*(x1_dot);
    
    % Applying pfl
    
    % torque_input
    u=(H(1,1)*x2dot_des)+(H(1,2)*(H(2,2)\(-G(2,1)-(H(2,1)*x2dot_des))))+(C(1,2)*theta1_dot);
    torque=[u,0]';
    tou(i)=u;
    
    % Getting values for next iteration
    q=[x1, theta1]';
    qdot=[x1_dot, theta1_dot]';
    q2dot=H\(B*u-C*qdot-G);
    X_dot(3,:)=q2dot(1,1);
    X_dot(4,:)=q2dot(2,1);
    X(3:4,1)=X(3:4,1)+(X_dot(3:4,1)*dt);
    X(1:2,1)=X(1:2,1)+(X(3:4,1)*dt);
    i=i+1;
    X_main(i,1)=X(1,1);
    X_main(i,2)=wrapToPi(X(2,1));
    X_main(i,3)=X(3,1);
    X_main(i,4)=X(4,1);
    % X_main=X_main';
end
%   plot(X_main(:,2),X_main(:,4));


 end
% for k=1:40000
% %       drawcartpend_bw(optimal(k,:),m2,m1,2);
%        drawcartpend(X_main(:,k),mp,mc,0.5);
% end