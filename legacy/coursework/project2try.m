L = 1;
n = 3;
h = L/n;
xm=0:h:L;
ym=0:h:L;
v=0.1;
lambda=0.000001;
[x,y]=meshgrid(xm,ym);
f=@(x) 4*x.*(1-x);
uboundary=f(0:h:1);

%---------------------------------------------------------------------------------------------------%
%construct the matrix for (A), with elimination of bc-> (n-1)^2*(n-1)^2;
%stencil is: [c3 c6 c3;
%             c2 c5 c2;
%             c1 c4 c1];
coeff = v/6;
c1ex = -coeff;
c2ex = -4*coeff;
c3ex = -coeff;
c4ex = 2*coeff;
c5ex = 8*coeff;
c6ex = 2*coeff;
il = diag(ones(1,n),-1); 
iu = diag(ones(1,n),1); 
Mc = diag(c2ex*ones(1,n),1) + diag(c5ex*ones(1,n+1)) + diag(c2ex*ones(1,n),-1);
Ml = diag(c1ex*ones(1,n),1) + diag(c4ex*ones(1,n+1)) + diag(c1ex*ones(1,n),-1);
Mu = diag(c3ex*ones(1,n),1) + diag(c6ex*ones(1,n+1)) + diag(c3ex*ones(1,n),-1);
K11 = kron(iu,Mu)+kron(il,Ml)+kron(eye(n+1),Mc);


for j = 1:n+1
    K11(1+(j-1)*(n+1),:)=0;
        K11(:,1+(j-1)*(n+1))=0;

    K11(1+(j-1)*(n+1),1+(j-1)*(n+1)) = 1; 
   % z = z+1;
   K11(n+1+(j-1)*(n+1),:)=0;
   K11(:,n+1+(j-1)*(n+1))=0;
    K11(n+1+(j-1)*(n+1),n+1+(j-1)*(n+1)) = 1;
end
K11(end-(n+1)+1:end,:)=0;
K11(:,end-(n+1)+1:end)=0;
for i=1:n+1
K11(end-i+1,end-i+1) = 1;
end

%------------------------------------------------------------------------------%
%---------------------------------------------------------------------------------------------------%
%construct the matrix for (A), with elimination of bc-> (n-1)^2*(n-1)^2;
%stencil is: [c3 c6 c3;
%             c2 c5 c2;
%             c1 c4 c1];
coeff = v/6;
c1ex = -coeff;
c2ex = 2*coeff;
c3ex = -coeff;
c4ex = -4*coeff;
c5ex = 8*coeff;
c6ex = -4*coeff;
il = diag(ones(1,n),-1); 
iu = diag(ones(1,n),1); 
Mc = diag(c2ex*ones(1,n),1) + diag(c5ex*ones(1,n+1)) + diag(c2ex*ones(1,n),-1);
Ml = diag(c1ex*ones(1,n),1) + diag(c4ex*ones(1,n+1)) + diag(c1ex*ones(1,n),-1);
Mu = diag(c3ex*ones(1,n),1) + diag(c6ex*ones(1,n+1)) + diag(c3ex*ones(1,n),-1);
K12 = kron(iu,Mu)+kron(il,Ml)+kron(eye(n+1),Mc);
for j = 1:n+1
    K12(1+(j-1)*(n+1),:)=0;
        K12(:,1+(j-1)*(n+1))=0;

    K12(1+(j-1)*(n+1),1+(j-1)*(n+1)) = 1; 
   % z = z+1;
   K12(n+1+(j-1)*(n+1),:)=0;
   K12(:,n+1+(j-1)*(n+1))=0;
    K12(n+1+(j-1)*(n+1),n+1+(j-1)*(n+1)) = 1;
end
K12(end-(n+1)+1:end,:)=0;
K12(:,end-(n+1)+1:end)=0;
for i=1:n+1
K12(end-i+1,end-i+1) = 1;
end
%.---------------------------------------------------------------------------%
K21=K11;
K22=K12;
%-----------------------------------------------------------------------------%

%---------------------------------------------------------------------------------------------------%
%construct the matrix for (E), without elimination of bc-> (n-1)^2*(n+1)^2;
%stencil is: [c3 c6 c3ex;
%             c2 c5 c2ex;
%             c1 c4 c1ex];
coeff = h/12;
c1=coeff;
c1ex = -coeff;
c2=4*coeff;
c2ex = -4*coeff;
c3=coeff;
c3ex = -coeff;
c4ex = 0;
c5ex =0;
c6ex = 0;
il = diag(ones(1,n),-1); 
iu = diag(ones(1,n),1); 
Mc = diag(c2ex*ones(1,n),1) + diag(c5ex*ones(1,n+1)) + diag(c2*ones(1,n),-1);
Ml = diag(c1ex*ones(1,n),1) + diag(c4ex*ones(1,n+1)) + diag(c1*ones(1,n),-1);
Mu = diag(c3ex*ones(1,n),1) + diag(c6ex*ones(1,n+1)) + diag(c3*ones(1,n),-1);
B1 = kron(iu,Mu)+kron(il,Ml)+kron(eye(n+1),Mc);
for j = 1:n+1
    B1(1+(j-1)*(n+1),:)=0;
        B1(:,1+(j-1)*(n+1))=0;

    B1(1+(j-1)*(n+1),1+(j-1)*(n+1)) = 1; 
   % z = z+1;
   B1(n+1+(j-1)*(n+1),:)=0;
   B1(:,n+1+(j-1)*(n+1))=0;
    B1(n+1+(j-1)*(n+1),n+1+(j-1)*(n+1)) = 1;
end
B1(end-(n+1)+1:end,:)=0;
B1(:,end-(n+1)+1:end)=0;
for i=1:n+1
B1(end-i+1,end-i+1) = 1;
end%--------------------------------------------------------------------------------------------------%


%---------------------------------------------------------------------------------------------------%
%construct the matrix for (E), without elimination of bc-> (n-1)^2*(n+1)^2;
%stencil is: [c3 c6 c3ex;
%             c2 c5 c2ex;
%             c1 c4 c1ex];
coeff = h/12;
c1ex = coeff;
c2ex =0.01*coeff;
c3ex = -coeff;
c4ex = 2*coeff;
c5ex =0;
c6ex = -2*coeff;
il = diag(ones(1,n),-1); 
iu = diag(ones(1,n),1); 
Mc = diag(c2ex*ones(1,n),1) + diag(c5ex*ones(1,n+1)) + diag(c2ex*ones(1,n),-1);
Ml = diag(c1ex*ones(1,n),1) + diag(c4ex*ones(1,n+1)) + diag(c1ex*ones(1,n),-1);
Mu = diag(c3ex*ones(1,n),1) + diag(c6ex*ones(1,n+1)) + diag(c3ex*ones(1,n),-1);
B2 = kron(iu,Mu)+kron(il,Ml)+kron(eye(n+1),Mc);
for j = 1:n+1
    B2(1+(j-1)*(n+1),:)=0;
        B2(:,1+(j-1)*(n+1))=0;

    B2(1+(j-1)*(n+1),1+(j-1)*(n+1)) = 1; 
   % z = z+1;
   B2(n+1+(j-1)*(n+1),:)=0;
   B2(:,n+1+(j-1)*(n+1))=0;
    B2(n+1+(j-1)*(n+1),n+1+(j-1)*(n+1)) = 1;
end
B2(end-(n+1)+1:end,:)=0;
B2(:,end-(n+1)+1:end)=0;
for i=1:n+1
B2(end-i+1,end-i+1) = 1;
end
%---------------------------------------------------------------------------------------------------%
%construct the matrix for (I), without elimination of bc-> (n+1)^2*(n+1)^2;
%stencil is: [c3 c6 c3ex;
%             c2 c5 c2ex;
%             c1 c4 c1ex];
coeff = h^2/36;
c1ex = coeff;
c2ex = 4*coeff;
c3ex = coeff;
c4ex = 4*coeff;
c5ex =16*coeff;
c6ex = 4*coeff;
il = diag(ones(1,n),-1); 
iu = diag(ones(1,n),1); 
Mc = diag(c2ex*ones(1,n),1) + diag(c5ex*ones(1,n+1)) + diag(c2ex*ones(1,n),-1);
Ml = diag(c1ex*ones(1,n),1) + diag(c4ex*ones(1,n+1)) + diag(c1ex*ones(1,n),-1);
Mu = diag(c3ex*ones(1,n),1) + diag(c6ex*ones(1,n+1)) + diag(c3ex*ones(1,n),-1);
M = kron(iu,Mu)+kron(il,Ml)+kron(eye(n+1),Mc);
M=-lambda*M;


%---------------------------------------------------------------------------------------------------%
%construct the global matrix and system
G=[K11 K12 B1; K21 K22 B2;B1' B2' M];
F=zeros(size(G,2),1);
% F((n-1)^2-(n-2):(n-1)^2)=2/3*h;
z=2;
for i=(n-1)^2-(n-2):(n-1)^2   
F(i)=v/6*(-2*uboundary(z-1)-2*uboundary(z)-2*uboundary(z+1));
z=z+1;
end
z=2;
% for i=2*(n-1)^2-(n-2):2*(n-1)^2  
% F(i)=v/6*(-uboundary(z-1)+2*uboundary(z)-uboundary(z+1));
% z=z+1;
% end
% z=2;
% for i=length(F)-(n-2):length(F)
% F(i)=v/6*(-uboundary(z-1)+2*uboundary(z)-uboundary(z+1));
% z=z+1;
% end
sol=G\F;

u1=sol(1:(n-1)^2);
u2=sol((n-1)^2+1:(size(G,2))-(n+1)^2);
p=sol(end-(n+1)^2+1:end);

U1=zeros((n-1),(n-1));
for i=1:n-1
    for j=1:n-1
        k=globalindex_n1(i,j,n);
        U1(i,j)=u1(k);
    end
end
U1b=zeros((n+1),(n+1));
U1b(2:end-1,2:end-1)=U1;
U1b(n+1,:)=f(x(n+1,:));

U2=zeros((n-1),(n-1));
for i=1:n-1
    for j=1:n-1
        k=globalindex_n1(i,j,n);
        U2(i,j)=u2(k);
    end
end
U2b=zeros((n+1),(n+1));
U2b(2:end-1,2:end-1)=U2;

quiver(x,y,U1b,U2b);
axis([-0.05 1.05 -0.05 1.05]);











function k = globalindex_n1(i, j, n)
    k = i + (j-1)*(n-1); 
end
