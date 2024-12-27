% ���ݼ���
function X=fanyan(YS,X,XK,ZK,E,NS)
 
% Ԥ�����
M=length(X);            % ģ�Ͳ�������
Y0=zhengyan(X,XK,ZK);   % ��ʼģ���쳣
f0=mbhs(YS,Y0);         % ƫ��ƽ���ͳ�ֵ
count=0;                % ��������������
D=100;                  % �������ӳ�ֵ
v=10;                   % �������Ӹı䱶��
sec=0.5;                % ��ͣʱ��
stop=false;             % �Ƿ����������ѭ�����Ʋ���(stop=true������)
h=figure(1);            % Create a figure window
set(h, 'position', get(0,'ScreenSize'));
 
% ��������
while (1)
    count=count+1;
    D=D/v;
    J=jacobi(X,XK,ZK);
    A=J*J';
    g=J*(YS-zhengyan(X,XK,ZK))';
    % ���
    T=A;
    for i=1:M
        for j=1:M
            A(i,j)=T(i,j)/sqrt(T(i,i)*T(j,j));
        end
        g(i)=g(i)/sqrt(T(i,i));
    end
    
    % �������������жϸôε����Ƿ��½�
    while (1)
        %����˹���ֽⷨ ����
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
    
    % �ж��Ƿ���ֹ����
    if (stop || (f0-f1)/f1<E || count>NS)
        break
    end
    X=X+DX';
    f0=f1;
    
    % ��ͼ
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
    legend('��ʼֵ','ʵ��ֵ','����ֵ')
    xlabel('X / m')
    ylabel('�����쳣��g')
    title('��ά��ģ�������쳣����')
    
    subplot(1,3,3)
    hold on
    plot(count,f1,'r.','markersize',12)
    legend(num2str(f1))
    grid on
    xlabel('��������')
    ylabel('ƫ��ƽ����')
    title('����������ʵ�����ߵķ��ϳ̶�')
    
    pause(sec)
    T1=Y1;
end
saveas(h,'figure','tif')
end
