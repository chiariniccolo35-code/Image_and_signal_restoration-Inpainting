function M = div(px,py)

% Compute the divergence of a vector (px,py) where
% - px has the same size as py
% BW differences with Dirichlet b.c.

[m,n]=size(px);

M=zeros(m,n);
Mx=M;
My=M;

Mx(2:m-1,1:n)=px(2:m-1,1:n)-px(1:m-2,1:n);
Mx(1,:)=px(1,:);
Mx(m,:)=-px(m-1,:);

My(1:m,2:n-1)=py(1:m,2:n-1)-py(1:m,1:n-2);
My(:,1)=py(:,1);
My(:,n)=-py(:,n-1);

M=Mx+My;

