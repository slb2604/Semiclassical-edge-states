function [] = CrankNicolsonfigeigth()
L = 230;
tsteps = 1800*pi; % number of time steps 
epsilon =0.02;
time_step=0.001;
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
        Gauss(i+(j-1)*(N)) = sqrt(r)/(sqrt(pi)).*exp(-(r*(epsilon*xl(j)-1)^2/(2*epsilon)+(epsilon*xl(i)-1)^2/(2*epsilon)));
        %Gauss(i+(j-1)*(N)) = sqrt(r)/(sqrt(pi)).*exp(-(r*(epsilon*xl(j))^2/(2*epsilon)+epsilon*r*xl(i)^2/(2)));
    end
end


exponent =  -1i * 2 * pi * kl' * xl ; 
A = exp( exponent ) / sqrt( N );
invA = inv(A);
for i=1:N
    Derivative_Operator(i,:)= 2*pi*1i*kl(i)*A(i,:);
    Position_Operatorm(:,i) = (epsilon*xl(i)-1)^2*A(:,i);
    Position_Operatorp(:,i) = (epsilon*xl(i)+1)^2*A(:,i);
    Position_Operator2(:,i) = xl(i)^2*A(:,i);
end   
Derivative_Operator = invA * Derivative_Operator;
Position_Operatorm = invA * Position_Operatorm;
Position_Operatorp = invA * Position_Operatorp;
Position_Operator2 = invA * Position_Operator2;
One = eye(size(Derivative_Operator));
XDer = kron(Derivative_Operator,One);
YDer = kron(One,Derivative_Operator);
kappa = r/2*(kron(Position_Operatorm,One)+epsilon^2*kron(One,Position_Operator2)-kron(One,One))...
    *r/2*(kron(Position_Operatorp,One)+epsilon^2*kron(One,Position_Operator2)-kron(One,One));
Hamiltonian = 1/epsilon*[kappa, -1i*XDer - YDer; -1i*XDer+YDer,-kappa];
Hamiltonian = (Hamiltonian + Hamiltonian')/2;
Time_Evol_Op_Up  = [kron(One,One),zeros(size(kron(One,One)));zeros(size(kron(One,One))),kron(One,One)] -1/2*1i*Hamiltonian*time_step;
Time_Evol_Op_Down  =[kron(One,One),zeros(size(kron(One,One)));zeros(size(kron(One,One))),kron(One,One)]+1/2*1i*Hamiltonian*time_step;
CrankNicolson = Time_Evol_Op_Up*inv(Time_Evol_Op_Down);
theta=0;
Z = transpose([Gauss*exp(-1i*theta/2),-Gauss*exp(1i*theta/2)]);

%Crank-Nicolson scheme
argu = linspace(0,0,tsteps);
argu2 = linspace(0,0,tsteps);
for i=1:tsteps 
    Z=CrankNicolson*Z;
   H=reshape(Z,[],2);
   B=H(:,1);
   nor = norm(B)^2;
   for j=1:length(H(:,1))
     argu(i) = argu(i)+angle(B(j))*abs(B(j))^2/nor;
   end
   C=H(:,2);
      nor2 = norm(C)^2;
   for j=1:length(H(:,2))
     argu2(i) = argu2(i)+angle(C(j))*abs(C(j))^2/nor2;
   end
    a1a(i,:,:) = reshape(abs(B).^2,N,N);
    a1b(i,:,:) = reshape(abs(C).^2,N,N);
end

teta=-pi:0.01:pi;
w1 = 1/epsilon+cos(teta)/epsilon;
w1b = -1/epsilon+cos(teta)/epsilon;
    w2 = sin(teta)/epsilon;
    hold on
for i=1:tsteps % Create Matrix Of Vectors
    if mod(i-1,300)==0
   plot3(w1,w2,zeros(1,numel(w1)),'LineWidth',5,'Color',[0,0.7,0.9])
   plot3(w1b,w2,zeros(1,numel(w1)),'LineWidth',5,'Color',[0,0.7,0.9])
   grid on 
s=surf(x,y,squeeze(a1b(i,:,:)));
shading interp
alpha(s,'z')
xlim([-2.25/epsilon,2.25/epsilon])
ylim([-2.25/epsilon,2.25/epsilon])
    end
end
hold off
end