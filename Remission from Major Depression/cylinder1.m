function [X,Y,Z]=cylinder1(P1,P2,R,N)
% 
% P1 - center of the first base of cylinder
% P2 - center of the second base of cylinder
% R - radius of cylinder
% N - number of sectors of the base of cylinder
 
if nargin<4
N=100;
end

P1=P1(:);
P2=P2(:);
P12=P2-P1;
A=[0 0 1;1 0 0;0 1 0]';
B=[];
z1=P12/norm(P12);
if z1==A(:,1)
    B=A;
else
     x1=cross(A(:,1),z1);
     x1=x1(:)/norm(x1);
     y1=cross(z1,x1);
     y1=y1(:)/norm(y1);
     B=[z1 x1 y1];
 end
 
 E=B*A^(-1);
 
 H=norm(P12);
 
alf=linspace(0,2*pi,N+1);
 r=[R;R];
 x=r*cos(alf);
 y=r*sin(alf);
 z=[0;H]*ones(1,N+1);
 
 xyz1=E*[x(1,:);y(1,:);z(1,:)];
 xyz2=E*[x(2,:);y(2,:);z(2,:)];
 xe=P1(1)*ones(2,N+1)+[xyz1(1,:);xyz2(1,:)];
 ye=P1(2)*ones(2,N+1)+[xyz1(2,:);xyz2(2,:)];
 ze=P1(3)*ones(2,N+1)+[xyz1(3,:);xyz2(3,:)];
 
 %{
 surf(xe,ye,ze)
 shading interp
 camlight
 grid on
 hold on
 axis equal

 
 X1=xe(1,:); Y1=ye(1,:); Z1=ze(1,:);
 fill3(X1,Y1,Z1,Z1)
 
 X2=xe(2,:); Y2=ye(2,:); Z2=ze(2,:);
 fill3(X2,Y2,Z2,Z2)
  %}
 X=xe; Y=ye; Z=ze;
 
 %------------------------------------------------
 
 % Example:
 %{
 ------------------------------------------------
 
 [X,Y,Z]=cylinder1([1 1 1],[10 5 5],2);
 surf(X,Y,Z)
 shading interp
 camlight
 grid on
 
 %}
