function [C,c]=Cost(Q,R,i,x0,u0)
syms x1 theta1 x1_dot theta1_dot u1;
x=[x1 theta1 x1_dot theta1_dot]';
% display(x0)
u=u1;
x_des=[0,pi,0,0]';
X=(x-x_des);
Cost=0.5*(X'*Q*X + u'*R*u);
dcdx=[diff(Cost,x1) diff(Cost,theta1) diff(Cost,x1_dot) diff(Cost,theta1_dot)]';
dcdu=diff(Cost,u1);
% d2cdx2=[diff(dcdx,x1) diff(dcdx,theta1) diff(dcdx,x1_dot) diff(dcdx,theta1_dot)];
d2cdx2=jacobian(dcdx,[x1,theta1,x1_dot,theta1_dot]);
d2cdu2=diff(dcdu,u1);
d2cdxdu=diff(dcdx,u1);
d2cdudx=[diff(dcdu,x1) diff(dcdu,theta1) diff(dcdu,x1_dot) diff(dcdu,theta1_dot)];
C=[d2cdx2,d2cdxdu ;...
    d2cdudx,d2cdu2];

c=[dcdx; dcdu];
 x1=x0(1,1);
 theta1=wrapToPi(x0(1,2));
 x1_dot=x0(1,3);
 theta1_dot=x0(1,4);
 u1=u0(1,1);
 C=double(subs(C));
 c=double(subs(c));
end