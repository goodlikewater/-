% яе©к╠х
function J=jacobi(X,XK,ZK)
 
M=length(X);
N=length(XK);
J=zeros(M,N);
Y=zhengyan(X,XK,ZK);
for m=1:M
    X1=X;
    X1(m)=X1(m)+0.1;
    Y1=zhengyan(X1,XK,ZK);
    for n=1:N
        J(m,n)=(Y1(n)-Y(n))/0.1;
    end
end
end
