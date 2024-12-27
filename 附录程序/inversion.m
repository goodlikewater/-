% 反演计算
function X=fanyan(YS,X,XK,ZK,E,NS)
 
% 预设参数
M=length(X);            % 模型参数个数
Y0=zhengyan(X,XK,ZK);   % 初始模型异常
f0=mbhs(YS,Y0);         % 偏差平方和初值
count=0;                % 迭代次数计数器
D=100;                  % 阻尼因子初值
v=10;                   % 阻尼因子改变倍数
sec=0.5;                % 暂停时长
stop=false;             % 是否跳出最外层循环控制参数(stop=true则跳出)
h=figure(1);            % Create a figure window
set(h, 'position', get(0,'ScreenSize'));
 
% 迭代反演
while (1)
    count=count+1;
    D=D/v;
    J=jacobi(X,XK,ZK);
    A=J*J';
    g=J*(YS-zhengyan(X,XK,ZK))';
    % 规格化
    T=A;
    for i=1:M
        for j=1:M
            A(i,j)=T(i,j)/sqrt(T(i,i)*T(j,j));
        end
        g(i)=g(i)/sqrt(T(i,i));
    end
    
    % 计算修正量并判断该次迭代是否下降
    while (1)
        %乔列斯基分解法 经典
	  DX=cholesky(A+D.*eye(M),g);
        for i=1:M
            DX(i)=DX(i)/sqrt(T(i,i));
        end
        Y1=zhengyan(X+DX',XK,ZK);
        f1=mbhs(YS,Y1);
        if (f1>f0)
            D=D*v;
            if (D>1e5)
                stop=true;
                break
            end
        else
            break
        end
    end
    
    % 判断是否终止迭代
    if (stop || (f0-f1)/f1<E || count>NS)
        break
    end
    X=X+DX';
    f0=f1;
    
    % 绘图
    subplot(1,3,[1,2])
    hold on
    if (count==1)
        plot(XK,Y0,'ko');
        plot(XK,YS,'ro');
        pause(sec)
    else
        plot(XK,T1,'kx')
    end
    plot(XK,Y1,'bx');
    legend('初始值','实测值','反演值')
    xlabel('X / m')
    ylabel('重力异常Δg')
    title('二维板模型重力异常反演')
    
    subplot(1,3,3)
    hold on
    plot(count,f1,'r.','markersize',12)
    legend(num2str(f1))
    grid on
    xlabel('迭代次数')
    ylabel('偏差平方和')
    title('理论曲线与实测曲线的符合程度')
    
    pause(sec)
    T1=Y1;
end
saveas(h,'figure','tif')
end
