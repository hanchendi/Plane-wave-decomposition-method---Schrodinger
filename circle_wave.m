
mode=3;
load('zero_point.mat');
k=zero_point(mode);

dx=1/N;
x=-1:dx:1;
y=-1:dx:1;
for i=1:length(x)
    xx(i,:)=x;
    yy(:,i)=x;
end
for i=1:length(x)
    for j=1:length(x)
        if xx(i,j)^2+yy(i,j)^2>=1;
            xx(i,j)=NaN;
            yy(i,j)=NaN;
        end
    end
end
b(1:N-1,1)=0;
b(N,1)=1;
for i=1:N-1;
    for j=1:N
        A(i,j)=cos(number(i,1)*k*cos(theta(j))+number(i,2)*k*sin(theta(j))+phi(j));
    end
end
for j=1:N
    A(N,j)=cos(test_inside(1)*k*cos(theta(j))+test_inside(2)*k*sin(theta(j))+phi(j));
end
[U,S,V] = svd(A);
T=S;
T(find(S~=0)) = 1./S(find(S~=0));
G = V * T' * U';
a=G*b;
zz=zeros(length(x),length(x));
for i=1:N
    zz=zz+a(i)*cos(xx.*k*cos(theta(i))+yy.*k*sin(theta(i))+phi(i));
end
figure()
mesh(xx,yy,zz);
axis([-1 1 -1 1])
axis off

load('number.mat');
rr=linspace(0,1,N);
theta1=linspace(0,2*pi,N);
theta1=theta1';
xx1=cos(theta1)*rr;
yy1=sin(theta1)*rr;
zz1=cos(number(mode,1)*theta1)*besselj(number(mode,1),rr*zero_point(mode));
figure()
mesh(xx1,yy1,zz1);
axis([-1 1 -1 1])
axis off
