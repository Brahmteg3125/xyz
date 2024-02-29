clear all
clc
A=[2 3 -1 4;1 -2 6 -7];
B=[8;-3];
C=[2 3 4 7];
n=size(A,2);
m=size(A,1);
nCm=nchoosek(n,m);
p=nchoosek(1:n,m);
sol=[];
for i=1:nCm
    X=zeros(n,1);
    A1=A(:,p(i,:));
    if det(A1)~=0
        X1=inv(A1)*B;
        if all(X1>=0)
            X(p(i,:))=X1;
            sol=[sol X];
        end
    end
end
Z=C*sol
[maxz indz]=max(Z)
ans=sol(:,indz)
fprintf('Objective function value is %f', max(Z))
display(ans)
