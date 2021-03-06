
% control

%% 固定するパラメータ
h=0.03;
Roff=0.5;
Rc=2;

%% input
x=-5; %車の初期座標
y=-30; %車の初期座標
theta=pi/180*(160); %車の初期座標における傾き角

%% データファイルの読み込み
filename='JIT_database1.csv';
G=csvread(filename);
n1=length(G);
plot(G(:,1),G(:,2));
%% 道路情報を含むファイルの読み込み
filename='road_information_fuji.xlsx';
B=xlsread(filename,'fuji2');
n2=length(B);

C=zeros(n2-1,1); %フラグ　1はfalse,0はsuccess
D=zeros(n2-1,1);
%% E=zeros(n2-1,1);

H=zeros(100,3,n2-1); %座標の格納場所
Twhole=zeros(n2-1,1); %Tcの格納場所
Vwhole=zeros(n2-1,1); %Vcの格納場所

%% 描画の範囲の変数定義
a1=min(B);
% b1=B(1,2);
A1=min(a1(1,1),x)-10;
B1=min(a1(1,2),y)-10;

%% 出力予測のための入力探索範囲
f=0.1; %error of X
g=0.1; %error of Y
%% s=pi/60; %error of phi
phi(n2-1,1)=zeros;
P=zeros(n1,8); % データ格納庫
w=0;
%% 道路情報数でforまわす
for i=1:n2-2
    %座標変換
    A=[cos(theta) sin(theta); -sin(theta) cos(theta)];    
    X3=x;
    Y3=y;
    X0=A*[X3; Y3]; %前の座標
    X4=B(i+1,1);
    Y4=B(i+1,2);
    X1=A*[X4; Y4]; %目標座標
    X42=B(i+1,1)-x;
    Y42=B(i+1,2)-y;    
    X12=A*[X42; Y42];
    X5=B(i+2,1);
    Y5=B(i+2,2);
    X2=A*[X5; Y5]; %次の座標
    %傾きの計算
    x1=[X0(1,1), X1(1,1), X2(1,1)];
    y1=[X0(2,1), X1(2,1), X2(2,1)];
    p=polyfit(x1,y1,2); %2次曲線への近似
    
    m=2*p(1,1)*X1(1,1)+p(1,2); %傾き
    theta1=atan(m); %傾きから角度を求める(x軸基準)
    if theta1>0 %y軸基準の角度に直す
        phi(i)=theta1-pi/2;
    else
        phi(i)=pi/2+theta1;
    end
    
    %出力予測のための、入力の探索
    V=0; %初期化
    T1=0;
    n=1; %count
    
    %データ抽出
    P(:,:)=[];
    for j=1:n1
    if (X12(1,1)+f)>=G(j,1) && (X12(1,1)-f)<=G(j,1) && (X12(2,1)+g)>=G(j,2) && (X12(2,1)-g)<=G(j,2) 
       P(n,1)=G(j,1); 
       P(n,2)=G(j,2);
       P(n,3)=G(j,3);
       P(n,4)=G(j,4); %T
       P(n,5)=G(j,5); %V
       
       %距離計算
       P(n,6)=realsqrt((X12(1,1)-P(n,1))^2+(X12(2,1)-P(n,2))^2+(phi(i)-P(n,3))^2);
       n=n+1;
    end
    end
    TF(i)=isempty(P);
     if TF(i)==1
        
         w=w+1;
         disp(w)
         continue %Pが空なら次のループへ
     end
   
    P=P(1:n-1,:);
    %sort
    dd=sortrows(P,6);
     %重み付き平均でTを求める
    %nn=round(n1/100);
    Tc=0; %初期化
    DD=0; %距離の逆数の総和の初期化
    if dd(1,6)==0 %距離が0
        Tc=dd(1,4);
    else
        for k=1:n-1 %最後がn+1で終わるため
            Tc=Tc+dd(k,4)/dd(k,6);
            DD=DD+1/dd(k,6);
        end
    Tc=Tc/DD;
    end
    Twhole(i,1)=Tc;
    %重み付き平均でTを求める
    
     %距離計算
     for j=1:n-1
       P(j,7)=realsqrt((X12(1,1)-P(j,1))^2+(X12(2,1)-P(j,2))^2+(phi(i)-P(j,3))^2+(Tc-P(j,4))^2);
     end
     dd1=sortrows(P,7);
       
    Vc=0; %初期化
    DD=0; %距離の逆数の総和の初期化
    if dd1(1,7)==0 %距離が0
        Vc=dd1(1,5);
    else
        for k=1:n-1
            Vc=Vc+dd1(k,5)/dd1(k,7);
            DD=DD+1/dd1(k,7);
        end
    Vc=Vc/DD;
    end
    Vwhole(i,1)=Vc;

    
   
    
%　　　微分方程式
    options=odeset('RelTol',1e-4,'AbsTol', [1e-4 1e-4 1e-4]);
    [T,X]=ode45(@car4, (0:h:Tc), [0;0;0], options, Vc,Tc,phi(i));
    
    %判定
    
    D(i)=sqrt((X12(1,1)-X(floor(Tc/h)+1,1))^2+(X12(2,1)-X(floor(Tc/h)+1,2))^2); %距離
%     E(i)=abs(X(floor(Tc/h)+1,3)-phi(i));
    
%     A1=[cos(-theta) sin(-theta); -sin(-theta) cos(-theta)];
    
   if D(i)<=1.2
        C(i)=0;
        %もとの座標軸に戻して座標を保存
        for j=1:floor(Tc/h)+1
            H(j,1:2,i)=(A\[X(j,1); X(j,2)]).'+[x y];
            H(j,3,i)=X(j,3)+theta;
        end
        
        %車の座標位置の変換（もと
        x=H(floor(Tc/h)+1,1,i); 
        y=H(floor(Tc/h)+1,2,i);
        theta=H(floor(Tc/h)+1,3,i);
        
     w=0; %wのリセット
     
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

%% 描画
figure
axis([A1,A2,B1,B2]);
plot(B(1:n2-1,1),B(1:n2-1,2),'-o','Color','r','LineWidth',1); %目標点
hold on;




for i=1:n2-2
    axis([A1,A2,B1,B2]);
   
    if C(i)==0
    for j=1:floor(Twhole(i,1)/h)+1
         C1 = [cos(H(j,3,i)) -sin(H(j,3,i)); sin(H(j,3,i)) cos(H(j,3,i))]*[1/2*Rc;sqrt(3)/2*Rc] + [H(j,1,i);H(j,2,i)];
         C2 = [cos(H(j,3,i)) -sin(H(j,3,i)); sin(H(j,3,i)) cos(H(j,3,i))]*[-1/2*Rc;sqrt(3)/2*Rc] + [H(j,1,i);H(j,2,i)];
         C3 = [cos(H(j,3,i)) -sin(H(j,3,i)); sin(H(j,3,i)) cos(H(j,3,i))]*[-1/2*Rc;-sqrt(3)/2*Rc] + [H(j,1,i);H(j,2,i)];
         C4 = [cos(H(j,3,i)) -sin(H(j,3,i)); sin(H(j,3,i)) cos(H(j,3,i))]*[1/2*Rc;-sqrt(3)/2*Rc] + [H(j,1,i);H(j,2,i)];
         
         if i==1 && j==1 %初期位置
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

