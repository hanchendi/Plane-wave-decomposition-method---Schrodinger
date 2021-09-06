clear
clc

N=500;
dk=0.001;
kmin=3.8;
kmax=3.9;
L=N;
delta_theta=2*pi/L;
theta=0:delta_theta:2*pi-delta_theta;
phi=2*pi*rand(N,1);

t=1;
for i=1:N
    number(t,1)=cos(theta(i));
    number(t,2)=sin(theta(i));
    t=t+1;
end

test_inside(1)=0.1234;
test_inside(2)=0.4321;
o=2*pi*rand(3);
for i=1:3
    test_boundary(i,1)=sin(o(i));test_boundary(i,2)=cos(o(i));
end

t=1;
b(1:N-1,1)=0;
b(N,1)=1;
for k=kmin:dk:kmax
    result(t,1)=k;
    result(t,2)=0;
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
    for i=1:3
        o=0;
        for j=1:N
            o=o+a(i)*cos(test_boundary(i,1)*k*cos(theta(j))+test_boundary(i,2)*k*sin(theta(j))+phi(j));
        end
        result(t,2)=result(t,2)+abs(o);
    end
    t=t+1;
    disp(k)
end
load('zero_point.mat')
plot(result(:,1),result(:,2),'b');hold on;
plot(zero_point,0,'*r');axis([kmin kmax 0 1*10^(12)])