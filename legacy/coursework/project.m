function [w,v1,v2,p]=project(n,v1bound,v2bound)
L = 1;
h = L/n;
xgr = linspace(0, L, n+1);
ygr = linspace(0, L, n+1);
[x,y] = meshgrid(xgr, ygr);
%system is [M C1 C2 0;   [w      [th
%           A1 0 0 G1;    v1  =  fh1
%           A2 0 0 G2;    v2     fh2
%           0 D1 D2 0;]   p]     gh]

%stencil is: [c3 c6 c3;
%             c2 c5 c2;
%             c1 c4 c1];

%---------------------------------------------------------------------------------------------------------%

%construct the vorticity matrix M
coeff = h^2/36;
c1ex = coeff;
c2ex = 4*coeff;
c3ex = coeff;
c4ex = 4*coeff;
c5ex = 16*coeff;
c6ex = 4*coeff;
il = diag(ones(1,n),-1); 
iu = diag(ones(1,n),1);
Mc = diag(c2ex*ones(1,n),1) + diag(c5ex*ones(1,n+1)) + diag(c2ex*ones(1,n),-1);
Ml = diag(c1ex*ones(1,n),1) + diag(c4ex*ones(1,n+1)) + diag(c1ex*ones(1,n),-1);
Mu = diag(c3ex*ones(1,n),1) + diag(c6ex*ones(1,n+1)) + diag(c3ex*ones(1,n),-1);
M = kron(iu,Mu)+kron(il,Ml)+kron(eye(n+1),Mc);

%construct the velocity 1 matrix C1
%stencil is: [c3 c6 c3;
%             c2 c5 c2;
%             c1 c4 c1];
coeff = h/6;
c1ex = 0;
c2ex = -coeff;
c3ex = coeff;
c4ex = 0;
c5ex = -4*coeff;
c6ex = 4*coeff;
il = diag(ones(1,n),-1); 
iu = diag(ones(1,n),1);
C1c = diag(c2ex*ones(1,n),1) + diag(c5ex*ones(1,n+1)) + diag(c2ex*ones(1,n),-1);
C1l = diag(c1ex*ones(1,n),1) + diag(c4ex*ones(1,n+1)) + diag(c1ex*ones(1,n),-1);
C1u = diag(c3ex*ones(1,n),1) + diag(c6ex*ones(1,n+1)) + diag(c3ex*ones(1,n),-1);
C1 = kron(iu,C1u)+kron(il,C1l)+kron(eye(n+1),C1c);
%eliminating the boundary conditions+removing extra basis-functions
C1f = C1;
C1 = C1(:,n+2:end); %eliminate v1 for i=1,...,n+1, j=0
C1s = C1;
k = 0;
for j=1:n
    C1(:,1+(j-1)*(n+1)-k)=[]; %eliminate v1 for i=1, j=2,...,n+1
    k = k+1;
end
k = 0;
for j=1:n
    C1(:,n+(j-1)*n - k)=[]; %eliminate v1 for i=n+1, j=2,...,n+1
    k = k+1;
end

%construct the velocity 2 matrix C2
%stencil is: [c3 c6 c3;
%             c2 c5 c2;
%             c1 c4 c1];
coeff = h/6;
c1ex = -coeff;
c2ex = -4*coeff;
c3ex = -coeff;
c4ex = coeff;
c5ex = 4*coeff;
c6ex = coeff;
il = diag(ones(1,n),-1); 
iu = diag(ones(1,n),1);
C2c = diag(c5ex*ones(1,n+1)) + diag(c2ex*ones(1,n),+1);
C2l = diag(c4ex*ones(1,n+1)) + diag(c1ex*ones(1,n),+1);
C2u = diag(c6ex*ones(1,n+1)) + diag(c3ex*ones(1,n),+1);
C2 = kron(iu,C2u)+kron(il,C2l)+kron(eye(n+1),C2c);
%eliminating the boundary conditions+removing extra basis-functions
C2f = C2;
k = 0;
for j=1:(n+1)
    C2(:,1+(j-1)*(n+1)-k)=[]; %eliminate v2 for i=1, j=1,...,n+1
    k = k+1;
end
C2s = C2;
C2 = C2(:,(n+1):end); %eliminate v2 for i=2,...,n+1, j=1
C2(:,end-n+1:end)=[]; %eliminate v2 for i=2,...,n+1, j=n+1

%construct the vorticity matrix A1
%stencil is: [c3 c6 c3;
%             c2 c5 c2;
%             c1 c4 c1];
coeff = h/6;
c1ex = -coeff;
c2ex = coeff;
c3ex = 0;
c4ex = -4*coeff;
c5ex = 4*coeff;
c6ex = 0;
il = diag(ones(1,n),-1); 
iu = diag(ones(1,n),1);
A1c = diag(c2ex*ones(1,n),1) + diag(c5ex*ones(1,n+1)) + diag(c2ex*ones(1,n),-1);
A1l = diag(c1ex*ones(1,n),1) + diag(c4ex*ones(1,n+1)) + diag(c1ex*ones(1,n),-1);
A1u = diag(c3ex*ones(1,n),1) + diag(c6ex*ones(1,n+1)) + diag(c3ex*ones(1,n),-1);
A1 = kron(iu,A1u)+kron(il,A1l)+kron(eye(n+1),A1c);
%eliminating the boundary conditions+removing extra basis-functions
A1f = A1;
A1 = A1(n+2:end,:); %eliminate w1 for i=1,...,n+1, j=0
A1s = A1;
k = 0;
for j=1:n
    A1(1+(j-1)*(n+1)-k,:)=[]; %eliminate w1 for i=1, j=2,...,n+1
    k = k+1;
end
k = 0;
for j=1:n
    A1(n+(j-1)*n - k,:)=[]; %eliminate w1 for i=n+1, j=2,...,n+1
    k = k+1;
end

%construct the vorticity matrix A2
%stencil is: [c3 c6 c3;
%             c2 c5 c2;
%             c1 c4 c1];
coeff = h/6;
c1ex = coeff;
c2ex = 4*coeff;
c3ex = coeff;
c4ex = -coeff;
c5ex = -4*coeff;
c6ex = -coeff;
il = diag(ones(1,n),-1); 
iu = diag(ones(1,n),1);
A2c = diag(c5ex*ones(1,n+1)) + diag(c2ex*ones(1,n),-1);
A2l = diag(c4ex*ones(1,n+1)) + diag(c1ex*ones(1,n),-1);
A2u = diag(c6ex*ones(1,n+1)) + diag(c3ex*ones(1,n),-1);
A2 = kron(iu,A2u)+kron(il,A2l)+kron(eye(n+1),A2c);
%eliminating the boundary conditions+removing extra basis-functions
A2f = A2;
k = 0;
for j=1:(n+1)
    A2(1+(j-1)*(n+1)-k,:)=[]; %eliminate w2 for i=1, j=1,...,n+1
    k = k+1;
end
A2s = A2;
A2 = A2((n+1):end,:); %eliminate w2 for i=2,...,n+1, j=1
A2(end-n+1:end,:)=[]; %eliminate w2 for i=2,...,n+1, j=n+1

%construct the pressure matrix G1
%stencil is: [c3 c6 c3;
%             c2 c5 c2;
%             c1 c4 c1];
coeff = h;
c1ex = 0;
c2ex = coeff;
c3ex = 0;
c4ex = 0;
c5ex = -coeff;
c6ex = 0;
il = diag(ones(1,n),-1); 
iu = diag(ones(1,n),1);
G1c = diag(c2ex*ones(1,n),1) + diag(c5ex*ones(1,n+1));
G1l = diag(c1ex*ones(1,n),1) + diag(c4ex*ones(1,n+1)) + diag(c1ex*ones(1,n),-1);
G1u = diag(c3ex*ones(1,n),1) + diag(c6ex*ones(1,n+1)) + diag(c3ex*ones(1,n),-1);
G1 = kron(iu,G1u)+kron(il,G1l)+kron(eye(n+1),G1c);
%eliminating the boundary conditions+removing extra basis-functions
G1f = G1;
G1 = G1(n+2:end,:); %eliminate w1 for i=1,...,n+1, j=0
G1s = G1;
k = 0;
for j=1:n
    G1(1+(j-1)*(n+1)-k,:)=[]; %eliminate w1 for i=1, j=2,...,n+1
    k = k+1;
end
k = 0;
for j=1:n
    G1(n+(j-1)*n - k,:)=[]; %eliminate w1 for i=n+1, j=2,...,n+1
    k = k+1;
end
G1t = G1;
G1 = G1(:, n+2:end); %eliminate p for j=1, i=1,...,n+1
k = 0;
for j=1:n
    G1(:,1+(j-1)*(n+1) - k)=[]; %eliminate p for i=1, j=2,...,n+1
    k = k+1;
end

%construct the pressure matrix G2
%stencil is: [c3 c6 c3;
%             c2 c5 c2;
%             c1 c4 c1];
coeff = h;
c1ex = 0;
c2ex = 0;
c3ex = 0;
c4ex = 0;
c5ex = -coeff;
c6ex = coeff;
il = diag(ones(1,n),-1); 
iu = diag(ones(1,n),1);
G2c = diag(c2ex*ones(1,n),1) + diag(c5ex*ones(1,n+1)) + diag(c2ex*ones(1,n),-1);
G2l = diag(c1ex*ones(1,n),1) + diag(c4ex*ones(1,n+1)) + diag(c1ex*ones(1,n),-1);
G2u = diag(c3ex*ones(1,n),1) + diag(c6ex*ones(1,n+1)) + diag(c3ex*ones(1,n),-1);
G2 = kron(iu,G2u)+kron(il,G2l)+kron(eye(n+1),G2c);
%eliminating the boundary conditions+removing extra basis-functions
G2f = G2;
k = 0;
for j=1:(n+1)
    G2(1+(j-1)*(n+1)-k,:)=[]; %eliminate w2 for i=1, j=1,...,n+1
    k = k+1;
end
G2s = G2;
G2 = G2((n+1):end,:); %eliminate w2 for i=2,...,n+1, j=1
G2(end-n+1:end,:)=[]; %eliminate w2 for i=2,...,n+1, j=n+1
G2t = G2;
G2 = G2(:, n+2:end); %eliminate p for j=1, i=1,...,n+1
k = 0;
for j=1:n
    G2(:,1+(j-1)*(n+1) - k)=[]; %eliminate p for i=1, j=2,...,n+1
    k = k+1;
end

%construct the velocity 1 matrix D1
%stencil is: [c3 c6 c3;
%             c2 c5 c2;
%             c1 c4 c1];
coeff = h;
c1ex = 0;
c2ex = -coeff;
c3ex = 0;
c4ex = 0;
c5ex = coeff;
c6ex = 0;
il = diag(ones(1,n),-1); 
iu = diag(ones(1,n),1);
D1c = diag(c5ex*ones(1,n+1)) + diag(c2ex*ones(1,n),-1);
D1l = diag(c1ex*ones(1,n),1) + diag(c4ex*ones(1,n+1)) + diag(c1ex*ones(1,n),-1);
D1u = diag(c3ex*ones(1,n),1) + diag(c6ex*ones(1,n+1)) + diag(c3ex*ones(1,n),-1);
D1 = kron(iu,D1u)+kron(il,D1l)+kron(eye(n+1),D1c);
%eliminating the boundary conditions+removing extra basis-functions
D1f = D1;
D1 = D1(:,n+2:end); %eliminate v1 for i=1,...,n+1, j=0
D1s = D1;
k = 0;
for j=1:n
    D1(:,1+(j-1)*(n+1)-k)=[]; %eliminate v1 for i=1, j=2,...,n+1
    k = k+1;
end
k = 0;
for j=1:n
    D1(:,n+(j-1)*n - k)=[]; %eliminate v1 for i=n+1, j=2,...,n+1
    k = k+1;
end
D1t = D1;
D1 = D1(n+2:end, :); %eliminate p for j=1, i=1,...,n+1
k = 0;
for j=1:n
    D1(1+(j-1)*(n+1) - k,:)=[]; %eliminate p for i=1, j=2,...,n+1
    k = k+1;
end

%construct the velocity 2 matrix D2
%stencil is: [c3 c6 c3;
%             c2 c5 c2;
%             c1 c4 c1];
coeff = h;
c1ex = 0;
c2ex = 0;
c3ex = 0;
c4ex = -coeff;
c5ex = coeff;
c6ex = 0;
il = diag(ones(1,n),-1); 
iu = diag(ones(1,n),1);
D2c = diag(c2ex*ones(1,n),1) + diag(c5ex*ones(1,n+1)) + diag(c2ex*ones(1,n),-1);
D2l = diag(c1ex*ones(1,n),1) + diag(c4ex*ones(1,n+1)) + diag(c1ex*ones(1,n),-1);
D2u = diag(c3ex*ones(1,n),1) + diag(c6ex*ones(1,n+1)) + diag(c3ex*ones(1,n),-1);
D2 = kron(iu,D2u)+kron(il,D2l)+kron(eye(n+1),D2c);
%eliminating the boundary conditions+removing extra basis-functions
D2f = D2;
k = 0;
for j=1:(n+1)
    D2(:,1+(j-1)*(n+1)-k)=[]; %eliminate v2 for i=1, j=1,...,n+1
    k = k+1;
end
D2s = D2;
D2 = D2(:,(n+1):end); %eliminate v2 for i=2,...,n+1, j=1
D2(:,end-n+1:end)=[]; %eliminate v2 for i=2,...,n+1, j=n+1
D2t = D2;
D2 = D2(n+2:end,:); %eliminate p for j=1, i=1,...,n+1
k = 0;
for j=1:n
    D2(1+(j-1)*(n+1) - k,:)=[]; %eliminate p for i=1, j=2,...,n+1
    k = k+1;
end
%---------------------------------------------------------------------------------------------------------%
%check that A=-C^T, D=-G^T
C = [C1 C2];
A = [A1; A2];
G = [G1; G2];
D = [D1 D2];
ans1 = any(A+C', 'all'); %check whether there are non-zero elements in A+C^T
ans2 = any(D+G', 'all'); %check whether there are non-zero elements in D+G^T

%construct the global matrix B of size 4n^2+1 x 4n^2+1
Z1=zeros((n+1)^2,(n)^2);
Z2=zeros(2*n*(n-1),2*n*(n-1));
Z3=zeros(n^2,(n+1)^2);
Z4=zeros(n^2,n^2);

B=[M C Z1;A Z2 G;Z3 D Z4];
%eliminating the line and column corresponding to p(1,1) to solve it for a
%full rank matrix B of size 4n^2 x 4n^2
L=B;
L(end-n^2+1,:)=[];
L(:,end-n^2+1)=[];

%construct the rhs
f1=zeros(4*n^2,1);
for i=1:n+1
    for j=1:n+1
        k=globalindex_n1(i,j,n);
        if (i==1 && j==1)
        th(k)= h*(v1bound(0,0)-v2bound(0,0));
        elseif (i==n+1 && j==1)
            th(k)= h*(v1bound(1,0)+v2bound(1,0));
        elseif (i==1 && j==n+1)
            th(k)= h*(-v1bound(0,1)-v2bound(0,1));
            elseif (i==n+1 && j==n+1)
            th(k)= h*(-v1bound(1,1)+v2bound(1,1));
            elseif (i==1 && j~=1 && j~=n+1)
            th(k)= -h*v2bound(0,j*h);
             elseif (i==n+1 && j~=1 && j~=n+1)
            th(k)= h*v2bound(1,j*h);
             elseif (j==1 && i~=1 && i~=n+1)
            th(k)= h*(v1bound(i*h,0));
             elseif (j==n+1 && i~=1 && i~=n+1)
            th(k)= -h*v1bound(i*h,1);
        else th(k)=0;
        end
    end
end
for i=1:(n+1)^2
    f1(i)=th(i);
end
% for i=(n+1)^2-(n):(n+1)^2
%     f1(i)=-h;
% end
% f1((n+1)^2)=-h/2;
% f1((n+1)^2-(n))=-h/2;

%compute the solution of the reduced system Lx=f
%solution=L\f;
solution2=L\f1;
%---------------------------------------------------------------------------------------------------------%

%take the vorticity
vorticity=solution2(1:(n+1)^2);

%take v1 and apply boundary conditions
velocity1=solution2((n+1)^2+1:n*(n-1)+(n+1)^2);
vel1=[];
 for k=1:n-1:n^2-n
    b=[0; velocity1(k:n+k-2); 0];
    vel1=[vel1;b];
 end

%take v2 and apply boundary conditions
velocity2=solution2(n*(n-1)+(n+1)^2+1:end-n^2+1);
vel2=[zeros(n,1); velocity2; zeros(n,1)];

%take the pressure and set first component to 0
pressure=solution2(end-n^2+2:end);
pressure=[0;pressure];

%from vector to matrix notation
for i = 1:n+1  
      for j = 1:n
          k = globalindex_n1(i,j,n);
          v1(j,i) = vel1(k);
      end
end
for i = 1:n
      for j = 1:n+1
          k = globalindex_n(i,j,n);
          v2(j,i) = vel2(k);
      end
end
for i = 1:n+1
      for j = 1:n+1
          k = globalindex_n1(i,j,n);
          w(j,i) = vorticity(k);
      end
end
for i = 1:n
      for j = 1:n
          k = globalindex_n(i,j,n);
          p(j,i) = pressure(k);
      end
end

%fix the dimension of v1 and v2 to make them coincide with the domain size
v1extratry=(v1(n-1,:)+v1(n,:));
v2extratry=(v2(:,n-1)+v2(:,n));
v1extratry(n)=0;
v1extratry(n+1)=0;
v1try=[v1; v1extratry];
v2try=[v2 v2extratry];


%fix the dimension of p to make it coincide with the domain size
%by adding the average of all rows/columns
pcolextra=p(:,n-1)+p(:,n);
prowextra=p(n-1,:)+p(n,:); 
pelextra=p(n-1,n-1)+p(n,n-1)+p(n-1,n)+p(n,n);
p=[p zeros(n,1);zeros(1,(n+1))];
p(:,n+1)=[pcolextra;pelextra];
p(n+1,1:end-1)=prowextra;

% for i=1:n
%     p1extra=p1extra+p(:,i);
%     p2extra=p2extra+p(i,:);
% end
% p(:,n+1)=p1extra;
% p(n+1,:)=p2extra;

%---------------------------------------------------------------------------------------------------------%
%plot the velocity field
figure(1)
quiver(x,y,v1try,v2try,1,'LineWidth',1.3);
axis equal
xlabel('x')
ylabel('y')

%plot solution
figure(2)
subplot(2,2,1);
surf(x,y,v1try);
legend('velocity 1');

subplot(2,2,2);
surf(x,y,v2try);
legend('velocity 2');

subplot(2,2,3);
surf(x,y,w);
legend('vorticity');

subplot(2,2,4);
surf(x,y,p);
legend('pressure');

%compute and plot divergence field
div=divergence(x,y,v1try,v2try);
figure(3);
surf(x,y,div);
legend('divergence of v');

%compute and plot postprocessed stream function
figure(4);
phi=v1try.*y-v2try.*x;
surf(x,y,phi);
legend('postprocessed stream function');

figure(5);
subplot(1,2,1);
contourf(x,y,v1try);
subplot(1,2,2);
contourf(x,y,v2try);

%global index function
function k = globalindex_n1(i, j, n)
    k = i + (j-1)*(n+1); 
end
function k=globalindex_n(i,j,n)
k=i+(j-1)*n;
end
%boundary functions for the velocity


end