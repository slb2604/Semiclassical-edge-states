function [] = edgestate_bent()
L = 60;
tsteps = 61; % number of time steps 
epsilon =0.1;
t = linspace(-2*pi,2*pi,tsteps); %time-steps
% number of sampling points
N = 100;
% sample rate
dl = L / N;
% highest frequency detectable
kmax= 1 / ( 2 * dl );

theta = atan(sech(t).^2);

% array of x values
xl = linspace( -L/2, L/2, N );
% array of k values
kl = linspace( -kmax, kmax, N );
xl2 = linspace( -epsilon*L/2, epsilon*L/2, N );
[x,y]=meshgrid(xl2,xl2);
%Initial state
Gauss = 1/(sqrt(sqrt(pi)))*exp(-xl.^2/2);
Gauss2 = kron(Gauss,Gauss);
fun1 = @(t) 1./sqrt(1+sech(t).^4);
fun2 = @(t) sech(t).^2./sqrt(1+sech(t).^4);

for k=1:tsteps
    y1(k) = integral(fun1,0,t(k));
    y2(k) = integral(fun2,0,t(k));
for i=1:N
    for j=1:N
        Gaussmom(i+(j-1)*(N)) = (xl(i)*sech(y1(k))^2*(xl(i)+2*xl(j)*sech(y1(k))^2-xl(i)*sech(y1(k))^4)*tanh(y1(k)))/(1+sech(y1(k))^4)^(3/2)*1/(sqrt(pi)).*exp(-(xl(i)^2+xl(j)^2)/2);
    end
end
            initstatea(k,:)=Gauss2*exp(-1i*theta(k)/2);% This is the Fourier representation of the inital state, 1st component
            initstateb(k,:)=-Gauss2*exp(1i*theta(k)/2);% This is the Fourier representation of the inital state, 2nd component
            initstatea2(k,:)=(sinh(t(k))/(1+cosh(t(k))^4)*Gauss2+Gaussmom )*exp(-1i*theta(k)/2);% This is the Fourier representation of T1 applied to the initial state, 1st component
            initstateb2(k,:)=(sinh(t(k))/(1+cosh(t(k))^4)*Gauss2+Gaussmom )*exp(1i*theta(k)/2);% This is the Fourier representation of T1 applied to the initial state, 2nd component
end

exponent =  -1i * 2 * pi * kl' * xl ; 
A = exp( exponent ) / sqrt( N );
invA = inv(A);
for i=1:N
    Derivative_Operator(i,:)= 2*pi*1i*kl(i)*A(i,:);
    Position_Operator(:,i) = xl(i)*A(:,i);
end
Derivative_Operator = invA*Derivative_Operator;
Position_Operator = invA*Position_Operator;
One = eye(size(Derivative_Operator)); 
XDer = kron(Derivative_Operator,One);
YDer = kron(One,Derivative_Operator);
XPos = kron(Position_Operator,One);
YPos = kron(One,Position_Operator);


for i=1:tsteps
    L11 = -1i* XDer*cos(theta(i)) -1i*YDer*sin(theta(i)) + -sech(y1(i))^2/sqrt(1+sech(y1(i))^4)*XPos +1/sqrt(1+sech(y1(i))^4)*YPos ;
    L12 = -1i* XDer - YDer ;
    L21 = -1i* XDer + YDer ;
    L22 = -1i*XDer*cos(theta(i)) -1i* YDer*sin(theta(i)) +sech(y1(i))^2/sqrt(1+sech(y1(i))^4)*XPos -1/sqrt(1+sech(y1(i))^4) * YPos ;
    L=[L11, L12 ; L21 , L22];
    initstate2 = [initstatea2(i,:),initstateb2(i,:)];
    A=-reshape(lsqr(L,transpose(initstate2),0.0001,500),[],2);
    B=A(:,1);
    C=A(:,1);
    a1a(i,:,:) = reshape(B,N,N);
    a1b(i,:,:) = reshape(C,N,N);
end

teta=-2*pi:0.01:2*pi;
for i=1:length(teta)
    w1(i) = integral(fun1,0,teta(i));
    w2(i) = integral(fun2,0,teta(i));
end
figure(1)
hold on
plot3(w1,w2,zeros(1,numel(w1)),'LineWidth',5,'Color',[0,0.7,0.9])
grid on
XY = [x(:) y(:)];   
for i=1:tsteps % Create Matrix Of Vectors
R=[cos(theta(i)) -sin(theta(i)); sin(theta(i)) cos(theta(i))]; %CREATE THE MATRIX
rotXY=XY*R'; %MULTIPLY VECTORS BY THE ROT MATRIX 
Xqr = reshape(rotXY(:,1), size(x,1), []);
Yqr = reshape(rotXY(:,2), size(y,1), []);
Xqrs = Xqr+y1(i);
Yqrs = Yqr+y2(i);
s=surf(Xqrs,Yqrs,squeeze(abs(a1a(i,:,:)).^2+abs(a1b(i,:,:)).^2));
shading interp
alpha(s,'z')
end
xlim([-5,5])
ylim([-1.5,1.5])
hold off
figure(2)
hold on
plot3(w1,w2,zeros(1,numel(w1)),'LineWidth',5,'Color',[0,0.7,0.9])
grid on
XY = [x(:) y(:)];   
for i=1:tsteps % Create Matrix Of Vectors
R=[cos(theta(i)) -sin(theta(i)); sin(theta(i)) cos(theta(i))]; %CREATE THE MATRIX
rotXY=XY*R'; %MULTIPLY VECTORS BY THE ROT MATRIX 
Xqr = reshape(rotXY(:,1), size(x,1), []);
Yqr = reshape(rotXY(:,2), size(y,1), []);
Xqrs = Xqr+y1(i);
Yqrs = Yqr+y2(i);
b1a(i,:,:) = reshape(initstatea(i,:),N,N);
b1b(i,:,:) = reshape(initstateb(i,:),N,N);
s=surf(Xqrs,Yqrs,squeeze(abs(epsilon*a1a(i,:,:)+b1a(i,:,:)).^2+abs(epsilon*a1b(i,:,:)+b1b(i,:,:) ).^2));
shading interp
alpha(s,'z')
end
xlim([-5,5])
ylim([-1.5,1.5])
hold off
