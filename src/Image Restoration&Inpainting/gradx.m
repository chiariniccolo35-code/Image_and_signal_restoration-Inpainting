function M=gradx(I)

% Compute the gradient of an image I along horizontal direction
% FW differences with Neumann b.c.


[m,n]=size(I);
M=zeros(m,n);
%FW
M(1:m-1,1:n)=-I(1:m-1,:)+I(2:m,:);
M(m,1:n)=zeros(1,n);

%BW
% M(2:m,1:n)=-I(1:m-1,:)+I(2:m,:);
% M(1,1:n)=zeros(1,n);
