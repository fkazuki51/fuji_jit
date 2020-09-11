
% control

%% ŒÅ’è‚·‚éƒpƒ‰ƒ[ƒ^
h=0.03;
Roff=0.5;
Rc=2;

%% input
x=-5; %ŽÔ‚Ì‰ŠúÀ•W
y=-30; %ŽÔ‚Ì‰ŠúÀ•W
theta=pi/180*(160); %ŽÔ‚Ì‰ŠúÀ•W‚É‚¨‚¯‚éŒX‚«Šp

%% ƒf[ƒ^ƒtƒ@ƒCƒ‹‚Ì“Ç‚Ýž‚Ý
filename='JIT_database1.csv';
G=csvread(filename);
n1=length(G);
plot(G(:,1),G(:,2));
%% “¹˜Hî•ñ‚ðŠÜ‚Þƒtƒ@ƒCƒ‹‚Ì“Ç‚Ýž‚Ý
filename='road_information_fuji.xlsx';
B=xlsread(filename,'fuji2');
n2=length(B);

C=zeros(n2-1,1); %ƒtƒ‰ƒO@1‚Ífalse,0‚Ísuccess
D=zeros(n2-1,1);
%% E=zeros(n2-1,1);

H=zeros(100,3,n2-1); %À•W‚ÌŠi”[êŠ
Twhole=zeros(n2-1,1); %Tc‚ÌŠi”[êŠ
Vwhole=zeros(n2-1,1); %Vc‚ÌŠi”[êŠ

%% •`‰æ‚Ì”ÍˆÍ‚Ì•Ï”’è‹`
a1=min(B);
% b1=B(1,2);
A1=min(a1(1,1),x)-10;
B1=min(a1(1,2),y)-10;

%% o—Í—\‘ª‚Ì‚½‚ß‚Ì“ü—Í’Tõ”ÍˆÍ
f=0.1; %error of X
g=0.1; %error of Y
%% s=pi/60; %error of phi
phi(n2-1,1)=zeros;
P=zeros(n1,8); % ƒf[ƒ^Ši”[ŒÉ
w=0;
%% “¹˜Hî•ñ”‚Åfor‚Ü‚í‚·
for i=1:n2-2
    %À•W•ÏŠ·
    A=[cos(theta) sin(theta); -sin(theta) cos(theta)];    
    X3=x;
    Y3=y;
    X0=A*[X3; Y3]; %‘O‚ÌÀ•W
    X4=B(i+1,1);
    Y4=B(i+1,2);
    X1=A*[X4; Y4]; %–Ú•WÀ•W
    X42=B(i+1,1)-x;
    Y42=B(i+1,2)-y;    
    X12=A*[X42; Y42];
    X5=B(i+2,1);
    Y5=B(i+2,2);
    X2=A*[X5; Y5]; %ŽŸ‚ÌÀ•W
    %ŒX‚«‚ÌŒvŽZ
    x1=[X0(1,1), X1(1,1), X2(1,1)];
    y1=[X0(2,1), X1(2,1), X2(2,1)];
    p=polyfit(x1,y1,2); %2ŽŸ‹Èü‚Ö‚Ì‹ßŽ—
    
    m=2*p(1,1)*X1(1,1)+p(1,2); %ŒX‚«
    theta1=atan(m); %ŒX‚«‚©‚çŠp“x‚ð‹‚ß‚é(xŽ²Šî€)
    if theta1>0 %yŽ²Šî€‚ÌŠp“x‚É’¼‚·
        phi(i)=theta1-pi/2;
    else
        phi(i)=pi/2+theta1;
    end
    
    %o—Í—\‘ª‚Ì‚½‚ß‚ÌA“ü—Í‚Ì’Tõ
    V=0; %‰Šú‰»
    T1=0;
    n=1; %count
    
    %ƒf[ƒ^’Šo
    P(:,:)=[];
    for j=1:n1
    if (X12(1,1)+f)>=G(j,1) && (X12(1,1)-f)<=G(j,1) && (X12(2,1)+g)>=G(j,2) && (X12(2,1)-g)<=G(j,2) 
       P(n,1)=G(j,1); 
       P(n,2)=G(j,2);
       P(n,3)=G(j,3);
       P(n,4)=G(j,4); %T
       P(n,5)=G(j,5); %V
       
       %‹——£ŒvŽZ
       P(n,6)=realsqrt((X12(1,1)-P(n,1))^2+(X12(2,1)-P(n,2))^2+(phi(i)-P(n,3))^2);
       n=n+1;
    end
    end
    TF(i)=isempty(P);
     if TF(i)==1
        
         w=w+1;
         disp(w)
         continue %P‚ª‹ó‚È‚çŽŸ‚Ìƒ‹[ƒv‚Ö
     end
   
    P=P(1:n-1,:);
    %sort
    dd=sortrows(P,6);
     %d‚Ý•t‚«•½‹Ï‚ÅT‚ð‹‚ß‚é
    %nn=round(n1/100);
    Tc=0; %‰Šú‰»
    DD=0; %‹——£‚Ì‹t”‚Ì‘˜a‚Ì‰Šú‰»
    if dd(1,6)==0 %‹——£‚ª0
        Tc=dd(1,4);
    else
        for k=1:n-1 %ÅŒã‚ªn+1‚ÅI‚í‚é‚½‚ß
            Tc=Tc+dd(k,4)/dd(k,6);
            DD=DD+1/dd(k,6);
        end
    Tc=Tc/DD;
    end
    Twhole(i,1)=Tc;
    %d‚Ý•t‚«•½‹Ï‚ÅT‚ð‹‚ß‚é
    
     %‹——£ŒvŽZ
     for j=1:n-1
       P(j,7)=realsqrt((X12(1,1)-P(j,1))^2+(X12(2,1)-P(j,2))^2+(phi(i)-P(j,3))^2+(Tc-P(j,4))^2);
     end
     dd1=sortrows(P,7);
       
    Vc=0; %‰Šú‰»
    DD=0; %‹——£‚Ì‹t”‚Ì‘˜a‚Ì‰Šú‰»
    if dd1(1,7)==0 %‹——£‚ª0
        Vc=dd1(1,5);
    else
        for k=1:n-1
            Vc=Vc+dd1(k,5)/dd1(k,7);
            DD=DD+1/dd1(k,7);
        end
    Vc=Vc/DD;
    end
    Vwhole(i,1)=Vc;

    
   
    
%@@@”÷•ª•û’öŽ®
    options=odeset('RelTol',1e-4,'AbsTol', [1e-4 1e-4 1e-4]);
    [T,X]=ode45(@car4, (0:h:Tc), [0;0;0], options, Vc,Tc,phi(i));
    
    %”»’è
    
    D(i)=sqrt((X12(1,1)-X(floor(Tc/h)+1,1))^2+(X12(2,1)-X(floor(Tc/h)+1,2))^2); %‹——£
%     E(i)=abs(X(floor(Tc/h)+1,3)-phi(i));
    
%     A1=[cos(-theta) sin(-theta); -sin(-theta) cos(-theta)];
    
   if D(i)<=1.2
        C(i)=0;
        %‚à‚Æ‚ÌÀ•WŽ²‚É–ß‚µ‚ÄÀ•W‚ð•Û‘¶
        for j=1:floor(Tc/h)+1
            H(j,1:2,i)=(A\[X(j,1); X(j,2)]).'+[x y];
            H(j,3,i)=X(j,3)+theta;
        end
        
        %ŽÔ‚ÌÀ•WˆÊ’u‚Ì•ÏŠ·i‚à‚Æ
        x=H(floor(Tc/h)+1,1,i); 
        y=H(floor(Tc/h)+1,2,i);
        theta=H(floor(Tc/h)+1,3,i);
        
     w=0; %w‚ÌƒŠƒZƒbƒg
     
    else
        C(i)=1;
        break;
   end
    disp(i)
end

csvwrite('JIT_distance.csv',D);

a2=max(B);
%% b2=B(n2,2);
ac2=H(floor(Tc/h)+1,1,n2-1);
bc2=H(floor(Tc/h)+1,2,n2-1);
A2=max(a2(1,1),ac2)+10;
B2=max(a2(1,2),bc2)+10;

%% •`‰æ
figure
axis([A1,A2,B1,B2]);
plot(B(1:n2-1,1),B(1:n2-1,2),'-o','Color','r','LineWidth',1); %–Ú•W“_
hold on;




for i=1:n2-2
    axis([A1,A2,B1,B2]);
   
    if C(i)==0
    for j=1:floor(Twhole(i,1)/h)+1
         C1 = [cos(H(j,3,i)) -sin(H(j,3,i)); sin(H(j,3,i)) cos(H(j,3,i))]*[1/2*Rc;sqrt(3)/2*Rc] + [H(j,1,i);H(j,2,i)];
         C2 = [cos(H(j,3,i)) -sin(H(j,3,i)); sin(H(j,3,i)) cos(H(j,3,i))]*[-1/2*Rc;sqrt(3)/2*Rc] + [H(j,1,i);H(j,2,i)];
         C3 = [cos(H(j,3,i)) -sin(H(j,3,i)); sin(H(j,3,i)) cos(H(j,3,i))]*[-1/2*Rc;-sqrt(3)/2*Rc] + [H(j,1,i);H(j,2,i)];
         C4 = [cos(H(j,3,i)) -sin(H(j,3,i)); sin(H(j,3,i)) cos(H(j,3,i))]*[1/2*Rc;-sqrt(3)/2*Rc] + [H(j,1,i);H(j,2,i)];
         
         if i==1 && j==1 %‰ŠúˆÊ’u
         line([C1(1);C2(1)],[C1(2);C2(2)],'Color','g','LineWidth',1);
         line([C2(1);C3(1)],[C2(2);C3(2)],'Color','b','LineWidth',1);
         line([C3(1);C4(1)],[C3(2);C4(2)],'Color','b','LineWidth',1);
         line([C4(1);C1(1)],[C4(2);C1(2)],'Color','b','LineWidth',1);
         hold on;
         t=0:0.1:360;
         plot(0.3*cos(t)+H(j,1,i),0.3*sin(t)+H(j,2,i),'b','LineWidth',1);
         pause(0.1);
         end
         
         if j==floor(Twhole(i,1)/h)+1
         line([C1(1);C2(1)],[C1(2);C2(2)],'Color','g','LineWidth',1);
         line([C2(1);C3(1)],[C2(2);C3(2)],'Color','b','LineWidth',1);
         line([C3(1);C4(1)],[C3(2);C4(2)],'Color','b','LineWidth',1);
         line([C4(1);C1(1)],[C4(2);C1(2)],'Color','b','LineWidth',1);
         hold on;
         t=0:0.1:360;
         plot(0.3*cos(t)+H(j,1,i),0.3*sin(t)+H(j,2,i),'b','LineWidth',1);
         pause(0.1);
%          else
%          line([C1(1);C2(1)],[C1(2);C2(2)],'Color','g','LineWidth',1);
%          line([C2(1);C3(1)],[C2(2);C3(2)],'Color','g','LineWidth',1);
%          line([C3(1);C4(1)],[C3(2);C4(2)],'Color','g','LineWidth',1);
%          line([C4(1);C1(1)],[C4(2);C1(2)],'Color','g','LineWidth',1);
%          hold on;
         end
         %t=0:0.1:360;
%          plot(0.3*cos(t)+H(j,1,i),0.3*sin(t)+H(j,2,i),'r','LineWidth',1);
%          pause(0.1);
    end
    else
        break;
    end
end

