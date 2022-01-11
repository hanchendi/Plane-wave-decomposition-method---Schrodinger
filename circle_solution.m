clear
clc

N=20;
M=20;
dx=0.01;
delta=10^(-13);
xx=0:dx:N;
n=length(xx);
for i=0:M
    s=1;
    yy=besselj(i,xx);
    for j=1:n-1
        if yy(j)*yy(j+1)<0
            dx1=dx/2;
            x0=xx(j);
            while dx1>delta
                if besselj(i,x0)*besselj(i,x0+dx1)<0
                    x0=x0;
                    dx1=dx1/2;
                else
                    x0=x0+dx1;
                    dx1=dx1/2;
                end
            end
            eigv(i+1+M,s)=x0;
            s=s+1;
        end
    end
    disp(i)
end

eigv1=eigv;
eigv(find(eigv==0))=NaN;
l1=length(eigv(1,:));
l2=length(eigv(:,1));
t=1;
s=1;
while t==1
    [a1,a2]=find(eigv==min(min(eigv)));
    if length(a1)~=1
        disp(a1)
        disp(a2)
    end
    for i=1:length(a1)
        number(s,1)=a1(i)-1-M;
        number(s,2)=a2(i);
        zero_point(s,1)=eigv(a1(i),a2(i));
        s=s+1;
        eigv(a1(i),a2(i))=NaN;
    end
    t=0;
    for i=1:l1
        for j=1:l2
            if isnan(eigv(j,i))
            else
                t=1;
            end
        end
    end
end
parameter(1)=length(zero_point);
parameter(2)=delta;

save([pwd,'/parameter.mat'],'parameter');
save([pwd,'/number.mat'],'number');
save([pwd,'/zero_point.mat'],'zero_point');