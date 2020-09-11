function dx = car4(~,x,Vc,T,phi)

dx = zeros(3,1);
dx(1) =-sin(x(3))*Vc;
dx(2) =cos(x(3))*Vc;
dx(3) = phi/T;
% dx(3)=-pi/2/T*phi*cos(pi/2/T*t);
