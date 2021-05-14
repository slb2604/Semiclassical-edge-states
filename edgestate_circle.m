function [] = edgestate()
L = 60;
tsteps = 4; % number of time steps 
epsilon =0.1;
t = linspace(-pi,pi/2,tsteps); %time-steps
% number of sampling points
N = 130;
% sample rate
dl = L / N;
% highest frequency detectable
kmax= 1 / ( 2 * dl );

theta = t-pi/2;

% array of x values
xl = linspace( -L/2, L/2, N );
% array of k values
kl = linspace( -kmax, kmax, N );
xl2 = linspace( -epsilon*L/2, epsilon*L/2, N );
[x,y]=meshgrid(xl2,xl2);
%Initial state
Gauss = 1/(sqrt(sqrt(pi)))*exp(-xl.^2/2);
Gauss2 = kron(Gauss,Gauss);

for i=1:N
    for j=1:N
        Gaussmom(i+(j-1)*(N)) = (xl(i)^2+xl(j)^2)/2*1/(sqrt(pi)).*exp(-(xl(i)^2+xl(j)^2)/2);
    end
end

for i=1:tsteps
            initstatea(i,:)=Gauss2*exp(-1i*theta(i)/2);% This is the Fourier representation of the inital state, 1st component
            initstateb(i,:)=-Gauss2*exp(1i*theta(i)/2);% This is the Fourier representation of the inital state, 2nd component
            initstatea2(i,:)=(-Gauss2/2+Gaussmom )*exp(-1i*theta(i)/2);% This is the Fourier representation of T1 applied to the initial state, 1st component
            initstateb2(i,:)=(-Gauss2/2+Gaussmom )*exp(1i*theta(i)/2);% This is the Fourier representation of T1 applied to the initial state, 2nd component
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
    y1(i) = cos(t(i));
    y2(i) = sin(t(i));
    L11 = -1i* XDer*cos(theta(i)) -1i*YDer*sin(theta(i)) + y1(i)*XPos + y2(i) * YPos ;
    L12 = -1i* XDer - YDer ;
    L21 = -1i* XDer + YDer ;
    L22 = -1i*XDer*cos(theta(i)) -1i* YDer*sin(theta(i)) - y1(i)*XPos - y2(i) * YPos ;
    L=[L11, L12 ; L21 , L22];
    initstate2 = [initstatea2(i,:),initstateb2(i,:)];
    A=-reshape(lsqr(L,transpose(initstate2),1e-5,500),[],2);
    B=A(:,1);
    C=A(:,1);
    a1a(i,:,:) = reshape(B,N,N);
    a1b(i,:,:) = reshape(C,N,N);
end
hold on
teta=-pi:0.01:pi;
    w1 = cos(teta);
    w2 = sin(teta);
plot3(w1,w2,zeros(1,numel(w1)),'LineWidth',5,'Color',[0,0.7,0.9])
grid on
XY = [x(:) y(:)];   
for i=1:tsteps % Create Matrix Of Vectors
R=[cos(theta(i)) -sin(theta(i)); sin(theta(i)) cos(theta(i))]; %CREATE THE MATRIX
theta(i)
rotXY=XY*R'; %MULTIPLY VECTORS BY THE ROT MATRIX 
Xqr = reshape(rotXY(:,1), size(x,1), []);
Yqr = reshape(rotXY(:,2), size(y,1), []);
Xqrs = Xqr+y1(i);
Yqrs = Yqr+y2(i);
s=surf(Xqrs,Yqrs,squeeze(abs(a1a(i,:,:)).^2+abs(a1b(i,:,:)).^2));
%b1a(i,:,:) = reshape(initstatea(i,:),N,N);
%b1b(i,:,:) = reshape(initstateb(i,:),N,N);
%s=surf(Xqrs,Yqrs,squeeze(abs(epsilon*a1a(i,:,:)+b1a(i,:,:)).^2+abs(epsilon*a1b(i,:,:)+b1b(i,:,:) ).^2));
shading interp
alpha(s,'z')
end
xlim([-2,2])
ylim([-2,2])
hold off
