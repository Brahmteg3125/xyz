%Simplex method
clc
clear all
cost=[4 10 0 0 0 0];
 a=[2 1 1 0 0;2 5 0 1 0;2 3 0 0 1];
 b=[50;100;90];
A=[a b];
 Var={'x1','x2','s1','s2','s3','sol'}
 bv=[3 4 5];

zjcj=cost(bv)*A-cost;
% Display initial simplex table
simplex_table=[zjcj; A];
array2table(simplex_table,'VariableNames',Var)
RUN=true;
 while RUN
if any(zjcj(1:end-1)<0) % check for negative value
 fprintf(' The current BFS is not optimal \n');
 zc=zjcj(1:end-1);
 [Enter_val, pvt_col]= min(zc) ;
 if all(A(:,pvt_col)<=0)
  error('LPP is Unbounded');
 else
 sol=A(:,end);
 column=A(:,pvt_col);
  for i=1:size(A,1)
 if column(i)>0   % pivot column value positive
 ratio(i)= sol (i)./column(i);
 else
 ratio(i)=inf;
 end
  end 
  [leaving_value,pvt_row]=min(ratio);
 end
 bv(pvt_row)=pvt_col;    % replaced leaving variable with entering variable
 pvt_key=A(pvt_row, pvt_col);
 A(pvt_row,:)=A (pvt_row,:)./pvt_key;
 % row operation 
for i=1:size(A,1)
 if i~=pvt_row
 A(i,:)=A(i,:)-A (i, pvt_col).*A(pvt_row,:);
 end
end
 zjcj=cost(bv)*A-cost;
 next_table=[zjcj; A];
array2table(next_table,'VariableNames',Var)

else
    RUN=false;
    fprintf('The table is optimal \n');
    z=input(' Enter 0 for minimization and 1 for max \n');
    if z==0
        Obj_value=-zjcj(end);
    else
        Obj_value=zjcj(end);
    end
    fprintf('The final optimal value is % f \n',Obj_value);
end
 end







 
% Big M method
clc
clear all
M=1000;

cost=[5 3 0 0 0 -M 0];
a=[1 1 1 0 0 0;5 2 0 1 0 0;2 8 0 0 -1 1];
b=[2;10;12];
artifical_var=[6];
bv=[3 4 6];

A=[a b];
Var={'x1','x2','s1','s2','s3','A1','sol'};

zjcj=cost(bv)*A-cost;
% Display initial simplex table
simplex_table=[zjcj; A];
array2table(simplex_table,'VariableNames',Var)
RUN=true;
 while RUN
if any(zjcj(1:end-1)<0) % check for negative value
 fprintf(' The current BFS is not optimal \n');
 zc=zjcj(1:end-1);
 [Enter_val, pvt_col]= min(zc) ;
 if all(A(:,pvt_col)<=0)
  error('LPP is Unbounded');
 else
 sol=A(:,end);
 column=A(:,pvt_col);
  for i=1:size(A,1)
 if column(i)>0
 ratio(i)= sol (i)./column(i);
 else
 ratio(i)=inf;
 end
  end 
  [leaving_value,pvt_row]=min(ratio);
 end
 bv(pvt_row)=pvt_col;
 pvt_key=A(pvt_row, pvt_col);
 A(pvt_row,:)=A (pvt_row,:)./pvt_key;
 % row operation 
for i=1:size(A,1)
 if i~=pvt_row
 A(i,:)=A(i,:)-A (i, pvt_col).*A(pvt_row,:);
 end
end
 zjcj=cost(bv)*A-cost;
 next_table=[zjcj; A];
array2table(next_table,'VariableNames',Var)

else
    RUN=false;
    if any(bv==artifical_var(1))
        error('Infeasible solution');
    else
    fprintf('The table is optimal \n');
    end
    z=input(' Enter 0 for minimization and 1 for max \n');
   
    if z==0
        Obj_value=-zjcj(end);
    else
        Obj_value=zjcj(end);
    end
    fprintf('The final optimal value is % f \n',Obj_value);
end
 end









least cost
clc
clear all
c=[2 10 4 5;6 12 8 11;3 9 5 7]
a=[12 25 20]
b=[25 10 15 5]
m=size(c,1)
n=size(c,2)
if sum(a)==sum(b)
    fprintf('Given problem is balanced')
else
    fprintf('Given problem is Unbalanced')
    if sum(a)<sum(b)
        c(end+1,:)=zeros(1,n)
        a(end+1)=sum(b)-sum(a)
        m=m+1
    else
        c(:,end+1)=zeros(m,1)
        b(end+1)=sum(a)-sum(b)
        n=n+1
    end
end
X=zeros(m,n)
OrinalC=c
for i=1:m
    for j=1:n
        cpq=min(min(c))
        if cpq==inf
            break
        end
        % or min(c(:))
        [p1,q1]=find(cpq==c)
        p=p1(1)
        q=q1(1)
        X(p,q)=min(a(p),b(q))
        if min(a(p),b(q))==a(p)
            b(q)=b(q)-a(p)
            a(p)=a(p)-X(p,q)
            c(p,:)=inf
        else
            a(p)=a(p)-b(q)
            b(q)=b(q)-X(p,q)
            c(:,q)=inf
        end
    end
end
sum=0
for i=1:m
    for j=1:n
        sum=sum+X(i,j)*OrinalC(i,j)
    end
end






steepset descent

clc
clear all
%Define Objective function
syms x1 x2 %symbolic 
%f(x1,x2)= x1-x2+2*x1^2+2*x1*x2+x2^2
%if we use @ it can't find derivative that's why we are using inline
f= x1-x2+2*x1^2+2*x1*x2+x2^2
fx=inline(f)
fobj=@(X) fx(X(1),X(2))
%f[0 0]
%X=[0 0]
%Find gradient of f
grad=gradient(f)
G=inline(grad)
Gx=@(X) G(X(1),X(2))
%Hessian matrix
H=hessian(f)
X0=[0 ,0]
tol=10^-3
maxiter=6
itr=0
%a=[a1,a2]
%norm(a)=sqrt((a1^2)+(a2^2))
while norm(Gx(X0))>tol && itr<maxiter
    S=-Gx(X0) 
    lambda=(S'*S)/(S'*H*S)
    Xnew=X0+S*lambda
    X0=Xnew
    itr=itr+1
end
fprintf('Optimal sol is %f',X0)



