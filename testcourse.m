for i=0:0.
    


X=2;
y=2^x;





%% 距離
for i=1:n
    d(i,6) = sqrt((Vc-d(i,1))^2+(To-d(i,2))^2+(Xo-d(i,3))^2+(Vr-d(i,4))^2);
end

%% % 距離が小さい順に行をソートする %%%
dd = sortrows(d,6);