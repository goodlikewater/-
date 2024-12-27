%Cholesky分解法经典
function [X]=cholesky(A,b) 
[N,N]=size(A); 
X=zeros(N,1); 
Y=zeros(N,1); 
for i=1:N 
   A(i,i)=sqrt(A(i,i)-A(i,1:i-1)*A(i,1:i-1)'); 
   if A(i,i)==0 
   'A is singular. no unique solution' 
      break 
   end 
   for j=i+1:N
      A(j,i)=(A(j,i)-A(j,1:i-1)*A(i,1:i-1)')/A(i,i); 
   end 
end 
%前代法
for j=1:N 
   Y(j)=(b(j)-A(j,1:j-1)*Y(1:j-1))/A(j,j); 
end 
A=A' 
for k=N:-1:1 
   X(k)=(Y(k)-A(k,k+1:N)*X(k+1:N))/A(k,k); 
end
