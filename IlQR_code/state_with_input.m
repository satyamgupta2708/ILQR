function [X]= state_with_input(mc,mp,l,g,x0,u0)
syms x1 theta1 x1_dot theta1_dot u1
% X_main=[];
x=[];
u=[u1 ;0];

H=[mc+mp , mp*l*cos(theta1);...
    mp*l*cos(theta1), mp*(l^2)];

C=[0 , -mp*l*theta1_dot*sin(theta1);...
    0, 0];

G=[0, mp*g*l*sin(theta1)]';

B = [1 0] ;
% Getting values for next iteration
q=[x1, theta1]';
qdot = [x1_dot; theta1_dot];
q2dot = H\(B*u-C*qdot-G);
% X_dot(3,:)=q2dot(1,1);
% X_dot(4,:)=q2dot(2,1);
% X(3:4,1)=X(3:4,1)+(X_dot(3:4,1)*dt);
% X(1:2,1)=X(1:2,1)+(X(3:4,1)*dt);
% i=i+1;
% X_main(i,1)=X(1,1);
% X_main(i,2)=wrapToPi(X(2,1));
% X_main(i,3)=X(3,1);
% X_main(i,4)=X(4,1);
 X = [qdot; q2dot];
 x1 = x0(1,1);
 theta1 = wrapToPi(x0(1,2));
 x1_dot = x0(1,3);
 theta1_dot = x0(1,4);
 u1=u0(1,1);
 X = double(subs(X));
end
% X_main=X_main';