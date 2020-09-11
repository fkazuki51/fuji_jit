filename='road_information_fuji.xlsx';
B=xlsread(filename,'fuji2');

n2=length(B);
N = NaN(2*n2-2,2);
l=length(N);
n=1;
for s=1:n2-1
    
    N(n+1,1)=(B(s,1)+B(s+1,1))/2;
    N(n,1:2)=B(s,1:2);
    n=n+2;
    
end    

linfill = fillmissing(N,'spline');
plot(linfill(:,1),linfill(:,2))

csvwrite('JIT_dataloadfull2.csv',linfill);

%‚Q”{
for H=1:3

filename2='JIT_dataloadfull2.csv';
B=csvread(filename2);

n2=length(B);
N = NaN(2*n2-2,2);
l=length(N);
n=1;
for s=1:n2-1
    
    N(n+1,1)=(B(s,1)+B(s+1,1))/2;
    N(n,1:2)=B(s,1:2);
    n=n+2;
    
end    

linfill = fillmissing(N,'spline');
plot(linfill(:,1),linfill(:,2))

csvwrite('JIT_dataloadfull2.csv',linfill);

end