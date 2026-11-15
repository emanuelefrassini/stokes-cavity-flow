L = 1;
n = 32;
h = L/n;
v=0.1;
lambda=1e-6;
xm=0:h:L;
ym=0:h:L;
%v=0.1;
%lambda=1e-7;
[x,y]=meshgrid(xm,ym);
f=@(x) 4*x.*(1-x);
uboundary=f(0:h:1);

%-----------------------------------------------00----------------------------------------------------%
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
%eliminate boundary conditions for K11
K11=K11(n+2:end,n+2:end);
z = 0;
for j = 1:n-1
    K11(:,1+(j-1)*(n+1)-z) = []; 
    K11(1+(j-1)*(n+1)-z,:) = []; 
   % z = z+1;
    K11(:,n+1+(j-1)*(n+1)-z) = [];
    K11(n+1+(j-1)*(n+1)-z,:) = [];
    z = z+1;
end
K11(end-(n+1)+1:end,:) = [];
K11(:,end-(n+1)+1:end) = [];


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
%eliminate boundary conditions for K11
K12=K12(n+2:end,n+2:end);
z = 0;
for j = 1:n-1
    K12(:,1+(j-1)*(n+1)-z) = []; 
    K12(1+(j-1)*(n+1)-z,:) = []; 
   % z = z+1;
    K12(:,n+1+(j-1)*(n+1)-z) = [];
    K12(n+1+(j-1)*(n+1)-z,:) = [];
    z = z+1;
end
K12(end-(n+1)+1:end,:) = [];
K12(:,end-(n+1)+1:end) = [];


%---------------------------------------------------------------------------------------------------%
%construct the matrix for (E), without elimination of bc-> (n-1)^2*(n+1)^2;
%stencil is: [c3 c6 c3ex;
%             c2 c5 c2ex;
%             c1 c4 c1ex];
coeff = -h/12;
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
% %eliminate boundary conditions for B1
% B1=B1(n+2:end,:);
% z = 0;
% for j = 1:n-1
%     B1(1+(j-1)*(n+1)-z,:) = []; 
%     B1(n+1+(j-1)*(n+1)-z,:) = [];
%     z = z+1;
% end
% B1(end-(n+1)+1:end,:) = [];
%eliminate boundary conditions for B1
B1=B1(:,n+2:end);
z = 0;
for j = 1:n-1
    B1(:,1+(j-1)*(n+1)-z) = []; 
    B1(:,n+1+(j-1)*(n+1)-z) = [];
    z = z+1;
end
B1(:,end-(n+1)+1:end) = [];
%----------------------------
% %--------------------------------------------------------------------------------------------------%


%---------------------------------------------------------------------------------------------------%
%construct the matrix for (E), without elimination of bc-> (n-1)^2*(n+1)^2;
%stencil is: [c3 c6 c3ex;
%             c2 c5 c2ex;
%             c1 c4 c1ex];
coeff = -h/12;
c1ex = coeff;
c2ex =0;
c3ex = -coeff;
c4ex = 4*coeff;
c5ex =0;
c6ex = -4*coeff;
il = diag(ones(1,n),-1); 
iu = diag(ones(1,n),1); 
Mc = diag(c2ex*ones(1,n),1) + diag(c5ex*ones(1,n+1)) + diag(c2ex*ones(1,n),-1);
Ml = diag(c1ex*ones(1,n),1) + diag(c4ex*ones(1,n+1)) + diag(c1ex*ones(1,n),-1);
Mu = diag(c3ex*ones(1,n),1) + diag(c6ex*ones(1,n+1)) + diag(c3ex*ones(1,n),-1);
B2 = kron(iu,Mu)+kron(il,Ml)+kron(eye(n+1),Mc);
%eliminate boundary conditions for B1
% B2=B2(n+2:end,:);
% z = 0;
% for j = 1:n-1
%     B2(1+(j-1)*(n+1)-z,:) = []; 
%     B2(n+1+(j-1)*(n+1)-z,:) = [];
%     z = z+1;
% end
% B2(end-(n+1)+1:end,:) = [];
B2=B2(:,n+2:end);
z = 0;
for j = 1:n-1
    B2(:,1+(j-1)*(n+1)-z) = []; 
    B2(:,n+1+(j-1)*(n+1)-z) = [];
    z = z+1;
end
B2(:,end-(n+1)+1:end) = [];

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
K=K11+K12;
Z=zeros((n-1)^2);
G=[K Z B1'; Z K B2';B1 B2 M]; %this should be correct


% G=[K Z -B1'; Z K -B2';-B1 -B2 M];
F=zeros(size(G,2),1);
% F((n-1)^2-(n-2):(n-1)^2)=2/3*h;
z=2;
for i=(n-1)^2-(n-2):(n-1)^2   
F(i)=-v/6*(-2*uboundary(z-1)-2*uboundary(z)-2*uboundary(z+1));
z=z+1;
end
z=2;
for i=length(F)-(n+1)-(n-2)-1:length(F)-(n+1)-1
F(i)=h/12*(uboundary(z-1)-uboundary(z+1));
z=z+1;
end
K2=[K Z;Z K];
Ktilde=diag(diag(K2));
B=[B1 B2];
R=M-B*(inv(Ktilde))*B';
DG=B*B';
D=B;
Z=zeros(size(B));
Du =diag(diag(K2));
Eu=-tril(K2,-1);
Fu=-triu(K2,1);
Dp=diag(diag(DG));
Ep=-tril(DG,-1);
Fp=-triu(DG,1);
rhop=max(abs(eig(Dp+Ep)));
rhou=max(abs(eig(Du+Eu)));
wp=2/(1+sqrt(1-rhop^2));
wu=2/(1+sqrt(1-rhou^2));

I=eye(size(B',1));
prec=[I, -B';Z,DG];
%prec=[I -B';Z DG];

%R=diag(diag(M))-
%% ITERATIVE GAUSS SEIDEL METHOD
%-------------------------------------------------------------------------------------%
res=1;
alpha=1.4; 
toll=1e-1;
xold_u=zeros(2*(n-1)^2,1);
xold_p=zeros((n+1)^2,1);
x0=zeros(length(F),1);
u0=x0(1:2*(n-1)^2);
p0=x0(2*(n-1)^2+1:end);
iter=1;
tic;
while res>toll
    %step 1
    r=F-G*x0;
    r1=r(1:2*(n-1)^2);
    r2=r(2*(n-1)^2+1:end);

    %step 2

    uprel = (Du-wu*Eu)\(Fu*xold_u)+(Du-wu*Eu)\r1;
   % uprel=(Du-w*Eu)\(Du*xold_u)+(Du-w*Eu)\(Fu*xold_u+w*F(1:2*(n-1)^2));
    %step 3
   pprel=(Dp-wp*Ep)\(Fp*xold_p)+(Dp-wp*Ep)\(D*uprel);
  %      pprel=(Dp-w*Ep)\(Dp*xold_p)+(Dp-w*Ep)\(Fp*xold_p+w*F(2*(n-1)^2+1:end));

    %step 4
    deltas=prec*[uprel;pprel];
    deltau=deltas(1:2*(n-1)^2);
    deltap=deltas(2*(n-1)^2+1:end);

    %step 5
    u=u0+alpha*deltau;
    p=p0+alpha*deltap;

    u0=u;
    p0=p;
    xold_u=uprel;
    xold_p=pprel;
    x0=[u0;p0];
    res(iter)=norm(r)/norm(F);
    iter=iter+1;
end
toc;
disp('Number of iterations is:');
disp(iter);
disp('H is:');
disp(h);
% PLOTS
%%




sol=x0;

u1=sol(1:(n-1)^2);
u2=sol((n-1)^2+1:(size(G,2))-(n+1)^2);
p=sol(end-(n+1)^2+1:end);

U1=zeros((n-1),(n-1));
for i=1:n-1
    for j=1:n-1
        k=globalindex_n1(i,j,n);
        U1(j,i)=u1(k);
    end
end
U1b=zeros((n+1),(n+1));
U1b(2:end-1,2:end-1)=U1;
U1b(n+1,:)=f(x(n+1,:));

U2=zeros((n-1),(n-1));
for i=1:n-1
    for j=1:n-1
        k=globalindex_n1(i,j,n);
        U2(j,i)=u2(k);
    end
end
U2b=zeros((n+1),(n+1));
U2b(2:end-1,2:end-1)=U2;
P=zeros(n+1);
for i=1:n+1
    for j=1:n+1
        k=globalindex_nplus(i,j,n);
        P(j,i)=p(k);
    end
end


figure(1);
subplot(2,2,1);
quiver(x,y,U1b,U2b);
axis([-0.05 1.05 -0.05 1.05]);
subplot(2,2,2);
contourf(x,y,U1b);
colorbar("eastoutside");

title('Velocity 1');
subplot(2,2,3);
contourf(x,y,U2b);
colorbar("eastoutside");

title('Velocity 2');
subplot(2,2,4);
surf(x,y,P);
colorbar("eastoutside");

title('Pressure');









    function k = globalindex_nplus(i, j, n)
    k = i + (j-1)*(n+1); 
end

function k = globalindex_n1(i, j, n)
    k = i + (j-1)*(n-1); 
end