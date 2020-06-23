function [A,B]= Linearize_dynamics(mc,mp,l,g,x0,u0)

syms x1 theta1 x1_dot theta1_dot u1;


% Defining dynamics of Cart_Pole
H=[mc+mp , mp*l*cos(theta1);...
    mp*l*cos(theta1), mp*(l^2)];

C=[0 , -mp*l*theta1_dot*sin(theta1);...
    0, 0];

G=[0, mp*g*l*sin(theta1)]';
b=[1 0]';
qdot=[x1_dot, theta1_dot]';
q2dot=H\(b*u1-C*qdot-G);
X_dot=[qdot ; q2dot];

 B=jacobian(X_dot,u1);

 A=jacobian(X_dot,[x1,theta1,x1_dot,theta1_dot]);
 x1=x0(1,1);
 theta1=wrapToPi(x0(1,2));
 x1_dot=x0(1,3);
 theta1_dot=x0(1,4);
 u1=u0(1,1);
 A = double(subs(A));
 B = double(subs(B));
 
end
% x1=x0(1);
% theta1=x0(1,2);
% x1_dot=x0(1,3);
% theta1_dot=x0(1,4);
% u=u0;
% A=subs(A);
% B=subs(B);
