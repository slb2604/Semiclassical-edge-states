function [] = CrankNicolsonbent_many()
L = 125;
tsteps = 1001; % number of time steps 
epsilon =0.1;
time_step=0.01;
% number of sampling points
N = 100;
% sample rate
dl = L / N;
% highest frequency detectable
kmax= 1 / ( 2 * dl );

% array of x values
xl = linspace( -L/2, L/2, N );

% array of k values
kl = linspace( -kmax, kmax, N );
[x,y]=meshgrid(xl,xl);
epsi = linspace(0,20,7);
%Initial state
r=1;
for k=1:7
for i=1:N
    for j=1:N
        Gauss(i+(j-1)*(N)) = sqrt(r)/(sqrt(pi)).*exp(-(r*(epsilon*xl(j)-4.5)^2/(2*epsilon)+(epsilon*xl(i))^2/(2*epsilon)));
        %Gauss(i+(j-1)*(N)) = sqrt(r)/(sqrt(pi)).*exp(-(r*(epsilon*xl(j))^2/(2*epsilon)+epsilon*r*xl(i)^2/(2)));
    end
end


exponent =  -1i * 2 * pi * kl' * xl ; 
A = exp( exponent ) / sqrt( N );
invA = inv(A);
AA = kron(A,A);
invAA = inv(AA);
for i=1:N
    Derivative_Operator(i,:)= 2*pi*1i*kl(i)*A(i,:);
    Position_Operator(:,i) = xl(i)*A(:,i);
end   
Derivative_Operator = invA * Derivative_Operator ;
One = eye(size(Derivative_Operator));
Position_Operator = invA * Position_Operator;
XDer = kron(Derivative_Operator,One);
YDer = kron(One,Derivative_Operator);
XPos = kron(Position_Operator,One);
YPos = kron(One,Position_Operator);
for j=1:N
for i=1:N
        kappa(:,i+(j-1)*N) = (epsilon*xl(i)-epsilon*sqrt(xl(j)^2+epsi(k)^2)+sqrt(4.5^2+epsilon^2*epsi(k)^2))*AA(:,i+(j-1)*N);
end
end
kappa = invAA*kappa;
Hamiltonian = 1/epsilon*[kappa, -1i*XDer - YDer; -1i*XDer+YDer,-kappa];
%Hamiltonian = 1/epsilon*[epsilon*YPos, -1i*XDer - YDer+Be(k)*epsilon*kron(One,Position_Operator); -1i*XDer+YDer+Be(k)*epsilon*kron(One,Position_Operator),-epsilon*YPos];
Hamiltonian = (Hamiltonian + Hamiltonian')/2;
Time_Evol_Op_Up  = [kron(One,One),zeros(size(kron(One,One)));zeros(size(kron(One,One))),kron(One,One)] -1/2*1i*Hamiltonian*time_step;
Time_Evol_Op_Down  =[kron(One,One),zeros(size(kron(One,One)));zeros(size(kron(One,One))),kron(One,One)]+1/2*1i*Hamiltonian*time_step;
CrankNicolson = Time_Evol_Op_Up*inv(Time_Evol_Op_Down);
theta = atan(4.5/sqrt(4.5^2+epsi(k)^2));
w= transpose([exp(-1i*theta),-exp(1i*theta)]);
Z = transpose([Gauss*w(1),Gauss*w(2)]);
%Crank-Nicolson scheme
X1 = 0;
Y1 = 0;
t(1)=0;
Exp(k,1)=0;
for i=1:tsteps 
    Z=CrankNicolson*Z;
   H=reshape(Z,[],2);
   B=H(:,1);
    nor = norm(B)^2;
   X2 = B'*XPos*B/nor; 
   Y2 = B'*YPos*B/nor;
   t(i+1)=i*time_step;
     X1 =X2;
     Y1 = Y2;
   C=H(:,2);
   if k==1
    a1a(i,:,:) = reshape(abs(B).^2+abs(C).^2,N,N);
    Exp(k,i+1) = max(max(a1a(i,:,:)));
   end
   if k==2
    a2a(i,:,:) = reshape(abs(B).^2+abs(C).^2,N,N);
    Exp(k,i+1) = max(max(a2a(i,:,:)));
   end
   if k==3
    a3a(i,:,:) = reshape(abs(B).^2+abs(C).^2,N,N);
    Exp(k,i+1) = max(max(a3a(i,:,:)));
   end
   if k==4
    a4a(i,:,:) = reshape(abs(B).^2+abs(C).^2,N,N);
    Exp(k,i+1) = max(max(a4a(i,:,:)));
   end
    if k==5
    a5a(i,:,:) = reshape(abs(B).^2+abs(C).^2,N,N);
    Exp(k,i+1) = max(max(a5a(i,:,:)));
    end
    if k==6
    a6a(i,:,:) = reshape(abs(B).^2+abs(C).^2,N,N);
    Exp(k,i+1) = max(max(a6a(i,:,:)));
    end
    if k==7
    a7a(i,:,:) = reshape(abs(B).^2+abs(C).^2,N,N);
    Exp(k,i+1) = max(max(a7a(i,:,:)));
   end
end
end
figure(1)
hold on
plot(t,Exp(1,:))
plot(t,Exp(4,:))
plot(t,Exp(7,:))
hold off

w1=linspace(-6/epsilon,6/epsilon);
w21 =sqrt(w1.^2+epsi(1)^2)-sqrt(4.5^2+epsilon^2*epsi(1)^2)/epsilon;
w22 =sqrt(w1.^2+epsi(2)^2)-sqrt(4.5^2+epsilon^2*epsi(2)^2)/epsilon;
w23 =sqrt(w1.^2+epsi(3)^2)-sqrt(4.5^2+epsilon^2*epsi(3)^2)/epsilon;
w24 =sqrt(w1.^2+epsi(4)^2)-sqrt(4.5^2+epsilon^2*epsi(4)^2)/epsilon;
w25 =sqrt(w1.^2+epsi(5)^2)-sqrt(4.5^2+epsilon^2*epsi(5)^2)/epsilon;
w26 =sqrt(w1.^2+epsi(6)^2)-sqrt(4.5^2+epsilon^2*epsi(6)^2)/epsilon;
w27 =sqrt(w1.^2+epsi(7)^2)-sqrt(4.5^2+epsilon^2*epsi(7)^2)/epsilon;
figure(2)
hold on
for i=1:tsteps % Create Matrix Of Vectors
    if mod(i-1,50)==0
   plot3(w1,w21-30,zeros(1,numel(w1)),'LineWidth',5,'Color',[0,0.7,0.9])
   plot3(w1,w22-20,zeros(1,numel(w1)),'LineWidth',5,'Color',[0,0.7,0.9])
   plot3(w1,w23-10,zeros(1,numel(w1)),'LineWidth',5,'Color',[0,0.7,0.9])
   plot3(w1,w24,zeros(1,numel(w1)),'LineWidth',5,'Color',[0,0.7,0.9])
   plot3(w1,w25+10,zeros(1,numel(w1)),'LineWidth',5,'Color',[0,0.7,0.9])
   plot3(w1,w26+20,zeros(1,numel(w1)),'LineWidth',5,'Color',[0,0.7,0.9])
   plot3(w1,w27+30,zeros(1,numel(w1)),'LineWidth',5,'Color',[0,0.7,0.9])
   grid on 
   colorMap = jet(256);
colormap(colorMap);   % Apply the colormap
colorbar;
s=surf(x,y-30,squeeze(a1a(i,:,:)));
s2=surf(x,y-20,squeeze(a2a(i,:,:)));
s3=surf(x,y-10,squeeze(a3a(i,:,:)));
s4=surf(x,y,squeeze(a4a(i,:,:)));
s5=surf(x,y+10,squeeze(a5a(i,:,:)));
s6=surf(x,y+20,squeeze(a6a(i,:,:)));
s7=surf(x,y+30,squeeze(a7a(i,:,:)));
shading interp
alpha(s,'z')
alpha(s2,'z')
alpha(s3,'z')
alpha(s4,'z')
alpha(s5,'z')
alpha(s6,'z')
alpha(s7,'z')
xlim([-6/epsilon,6/epsilon])
ylim([-8/epsilon,4/epsilon])
    end
end
hold off
end