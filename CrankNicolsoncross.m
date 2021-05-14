function [] = CrankNicolsoncross()
L = 250;
tsteps = 901; % number of time steps 
epsilon =0.02;
time_step=0.005;
% number of sampling points
N = 135;
% sample rate
dl = L / N;
% highest frequency detectable
kmax= 1 / ( 2 * dl );

% array of x values
xl = linspace( -L/2, L/2, N );

% array of k values
kl = linspace( -kmax, kmax, N );

[x,y]=meshgrid(xl,xl);
%Initial state
r=1;

for i=1:N
    for j=1:N
        Gauss(i+(j-1)*(N)) = sqrt(r)/(sqrt(pi)).*exp(-(r*(epsilon*xl(j)-2)^2/(2*epsilon)+(epsilon*xl(i))^2/(2*epsilon)));
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
end   
Derivative_Operator = invA * Derivative_Operator ;
One = eye(size(Derivative_Operator));
XDer = kron(Derivative_Operator,One);
YDer = kron(One,Derivative_Operator);
for j=1:N
for i=1:N
        kappa(:,i+(j-1)*N) =epsilon^2*xl(i)*xl(j)*AA(:,i+(j-1)*N);
end
end
kappa = invAA*kappa;
Hamiltonian = 1/epsilon*[kappa, -1i*XDer - YDer; -1i*XDer+YDer,-kappa];
Hamiltonian = (Hamiltonian + Hamiltonian')/2;
Time_Evol_Op_Up  = [kron(One,One),zeros(size(kron(One,One)));zeros(size(kron(One,One))),kron(One,One)] -1/2*1i*Hamiltonian*time_step;
Time_Evol_Op_Down  =[kron(One,One),zeros(size(kron(One,One)));zeros(size(kron(One,One))),kron(One,One)]+1/2*1i*Hamiltonian*time_step;
CrankNicolson = Time_Evol_Op_Up*inv(Time_Evol_Op_Down);
theta=0;
Z = transpose([Gauss*exp(-1i*theta/2),-Gauss*exp(1i*theta/2)]);
%Crank-Nicolson scheme
for i=1:tsteps 
    i
    Z=CrankNicolson*Z;
   H=reshape(Z,[],2);
   B=H(:,1);
   C=H(:,2);
    a1a(i,:,:) = reshape(abs(B).^2,N,N);
    a1b(i,:,:) = reshape(abs(C).^2,N,N);
end
w1=linspace(-4/epsilon,4/epsilon);
w2 = 0*w1;
    hold on
for i=1:tsteps % Create Matrix Of Vectors
    if mod(i-1,30)==0
   plot3(w1,w2,zeros(1,numel(w1)),'LineWidth',5,'Color',[0,0.7,0.9])
    plot3(w2,w1,zeros(1,numel(w1)),'LineWidth',5,'Color',[0,0.7,0.9])
   grid on 
s=surf(x,y,squeeze(a1b(i,:,:)));
shading interp
alpha(s,'z')
xlim([-2.5/epsilon,2.5/epsilon])
ylim([-2.5/epsilon,2.5/epsilon])
    end
end
hold off
end