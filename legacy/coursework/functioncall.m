n=8;
v1bound = @(x,y) (y==1);
v2bound=@(x,y) 0*x;
[w,v1,v2,p]=project(n,v1bound,v2bound);

