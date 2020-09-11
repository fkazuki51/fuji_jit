%data base

%固定するパラメータ
h=0.02;
Roff=0.5;
Rc=2;


%変化させるパラメータ
%input

phimin=-pi;
phih=pi/90;
phimax=pi;

%output
Tmin=0.04; %s
Th=0.02;
Tmax=2;

Vcmin=16; %m/s
Vch=0.05;
Vcmax=18;

data=zeros(floor((phimax-phimin+phih)/phih)*floor((Tmax-Tmin+Th)/Th)*floor((Vcmax-Vcmin+Vch)/Vch),5);
n=1;

for phi=phimin:phih:phimax
    for T=Tmin:Th:Tmax
        for Vc=Vcmin:Vch:Vcmax

            %微分方程式
            options=odeset('RelTol',1e-4,'AbsTol',[1e-4 1e-4 1e-4]);
            [T1,X1]=ode45(@car4, (0:h:T), [0;0;0], options, Vc, T, phi);


                data(n,1)=X1(floor(T/h)+1,1);
                data(n,2)=X1(floor(T/h)+1,2);
                data(n,3)=phi;
                data(n,4)=T;
                data(n,5)=Vc;
                n=n+1;
        end
    end
end

data=data(1:n-1,:);
csvwrite('JIT_database1.csv',data);
n=n-1;



