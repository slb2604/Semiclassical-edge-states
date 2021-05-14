function [] = CrankNicolsonbent_line()
epsilon =0.1;
L = 11/epsilon;
tsteps =1001; % number of time steps 
time_step=0.01;
times = linspace(0,tsteps*time_step,tsteps);
% number of sampling points
N =100;
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
Nb=4;
Be= linspace(-0.9,0,Nb);
for p=1:Nb
    B=0;
    B_strength=Be(p)
sigma = sqrt(1+B^2);
for i=1:N
    for j=1:N
        Gauss(i+(j-1)*(N)) = exp(-((epsilon*xl(j)-3*pi/2)^2/(2*epsilon)+sigma*(epsilon*xl(i))^2/(2*epsilon)));
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
        kappa(:,i+(j-1)*N) = (1+B_strength*sin(xl(j)*epsilon))*epsilon*xl(i)*AA(:,i+(j-1)*N);
        A_mag_1(:,i+(j-1)*N) = 0*AA(:,i+(j-1)*N);
        A_mag_2(:,i+(j-1)*N) = 0*AA(:,i+(j-1)*N);
        %*exp(-1/(1-epsilon^2*xl(j)^2))
end
end
kappa = invAA*kappa;
A_mag_1 = invAA*A_mag_1;
A_mag_2 = invAA*A_mag_2;
Hamiltonian = 1/epsilon*[kappa, -1i*XDer-A_mag_1 - YDer+1i*A_mag_2; -1i*XDer-A_mag_1+YDer-1i*A_mag_2,-kappa];
Hamiltonian = (Hamiltonian + Hamiltonian')/2;
Time_Evol_Op_Up  = [kron(One,One),zeros(size(kron(One,One)));zeros(size(kron(One,One))),kron(One,One)] -1/2*1i*Hamiltonian*time_step;
Time_Evol_Op_Down  =[kron(One,One),zeros(size(kron(One,One)));zeros(size(kron(One,One))),kron(One,One)]+1/2*1i*Hamiltonian*time_step;
CrankNicolson = Time_Evol_Op_Up*inv(Time_Evol_Op_Down);
theta=0;
Z = transpose([Gauss*exp(-1i*theta/2),-Gauss*exp(1i*theta/2)]);
%Crank-Nicolson scheme
for i=1:tsteps 
    Z=CrankNicolson*Z;
   H=reshape(Z,[],2);
   B=H(:,1);
   C=H(:,2);
   w0=linspace(-5/epsilon,5/epsilon);

   if p==1
    a1a(i,:,:) = reshape(abs(B).^2,N,N);
    a1b(i,:,:) = reshape(abs(C).^2,N,N);
    w1 = 0*w0-1.5/epsilon;
   end
   if p==2
    a2a(i,:,:) = reshape(abs(B).^2,N,N);
    a2b(i,:,:) = reshape(abs(C).^2,N,N);
    w2 = 0*w0-0.5/epsilon;
   end
   if p==3
    a3a(i,:,:) = reshape(abs(B).^2,N,N);
    a3b(i,:,:) = reshape(abs(C).^2,N,N);
    w3 = 0*w0+0.5/epsilon;
   end
   if p==4
    a4a(i,:,:) = reshape(abs(B).^2,N,N);
    a4b(i,:,:) = reshape(abs(C).^2,N,N);
    w4 = 0*w0+1.5/epsilon;
   end
end
end
figure(1)
    hold on
for i=1:tsteps % Create Matrix Of Vectors
    if mod(i-1,60)==0
   plot3(w0,w1,zeros(1,numel(w1)),'LineWidth',5,'Color',[0,0.7,0.9])
   plot3(w0,w2,zeros(1,numel(w1)),'LineWidth',5,'Color',[0,0.7,0.9])
   plot3(w0,w3,zeros(1,numel(w1)),'LineWidth',5,'Color',[0,0.7,0.9])
   plot3(w0,w4,zeros(1,numel(w1)),'LineWidth',5,'Color',[0,0.7,0.9])
   grid on 
s1=surf(x,y-1.5/epsilon,sqrt(squeeze(a1b(i,:,:))+squeeze(a1a(i,:,:))));
s2=surf(x,y-0.5/epsilon,sqrt(squeeze(a2b(i,:,:))+squeeze(a2a(i,:,:))));
s3=surf(x,y+0.5/epsilon,sqrt(squeeze(a3b(i,:,:))+squeeze(a3a(i,:,:))));
s4=surf(x,y+1.5/epsilon,sqrt(squeeze(a4b(i,:,:))+squeeze(a4a(i,:,:))));
shading interp
alpha(s1,'z')
alpha(s2,'z')
alpha(s3,'z')
alpha(s4,'z')
xlim([-5/epsilon,5/epsilon])
ylim([-5/epsilon,5/epsilon])
    end
end
hold off
end