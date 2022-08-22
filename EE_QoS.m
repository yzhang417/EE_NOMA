clc;
clear;
%close all;

%circulation
LOOP = 5000;
%Power of Noise
sigmadBm = -70;
sigma = 10^(sigmadBm/10)*10^(-3); 
%path loss factor
a = 3.0;
%dm
dm = 80;
%de
de = 80;
%M
M_array = [2 3 4];
%ptot
PtotdBm = 20; 
Pc_dBm = 30;
Pc = 10^(Pc_dBm/10)*10^(-3);
%Q
Qmin = 0;
Qmax = 5;
dQ = 0.25;
%Result Array
len = length(Qmin:dQ:Qmax);
Rs = zeros(LOOP,len);
AverageRs = zeros(4,len);

%% NOMA
for index=1:1:3
    index = index
    % Number of licit user
    M = M_array(index);
    %Ptot
    Ptot = 10.^(PtotdBm/10)*10^(-3);  
    for loop = 1:1:LOOP     
       %% Channel Information
        h = 1/sqrt(2)*randn(1,M)+1/sqrt(2)*randn(1,M)*1j;
        h = h.*conj(h);
        nh = sort(h)./(1+dm.^a);                                     
        Qindex = 0;
        for QoS = Qmin : dQ : Qmax
            Qindex = Qindex+1;  
            Q = QoS*ones(1,M);
            P = zeros(1,M);                   
            % Pmin        
            B = (2.^(Q)-1)./nh;
            for m = M : -1 : 1
                P(m) = B(m)*(nh(m)*sum(P(m+1:M)) + sigma); 
            end           
            Pmin = sum(P);
            PmindBm = 10*log10(Pmin*10^3);  
            if PmindBm > PtotdBm
                PmindBm;
                Rs(loop,Qindex) = 0;
            else
                ropt = zeros(1,M);
                A = (2.^(Q)-1)./(nh*Ptot);
                for m = 1 : 1 : M-1
                    ropt(m) = A(m)*(Ptot*nh(m)*(1-sum(ropt(1:m-1))) + sigma)/(2^(QoS)); 
                end
                ropt(M) = 1-sum(ropt(1:M-1));
                Rbm = zeros(1,M);
                for m = 1:1:M
                    Rbm(m) = log2( 1 + (Ptot*nh(m)*ropt(m))/(Ptot*nh(m)*sum(ropt(m+1:M))+sigma) );
                end
                Rs(loop,Qindex) = sum(Rbm);
            end
        end       
    end
    AverageRs(index,:) = mean(Rs)/(Ptot+Pc);
end

%% TDMA
for loop=1:1:LOOP
    Qindex = 0;
    % Channel gain of licit user 
    h = 1/sqrt(2)*randn(1,1)+1/sqrt(2)*randn(1,1)*1j;
    nh = h.*conj(h)./(1+dm.^a);                                    
    
    for QoS = Qmin : dQ : Qmax
        QoS;
        Qindex = Qindex + 1;
        Pmin = (2^QoS-1)*sigma/nh;
        if Pmin>Ptot
            Rs(loop,Qindex) = 0;
        else
            Rs(loop,Qindex) = log2(1+Ptot*nh/sigma);
        end
    end
end
AverageRs(4,:) = nanmean(Rs)/(Ptot+Pc);

%% Plot
figure;
Ptotrange = Qmin : dQ : Qmax;
plot(Ptotrange,AverageRs(3,:),'-<k'); hold on;
plot(Ptotrange,AverageRs(2,:),'-ob'); hold on;
plot(Ptotrange,AverageRs(1,:),'-*r'); hold on;
plot(Ptotrange,AverageRs(4,:),'-<g'); hold on;
grid on;
xlabel('$R_m^\textrm{Min}$ (bits/s/Hz)');ylabel('Average EE (bits/Joule/Hz)');
legend('NOMA, K=4, Energy-Efficient PA', ...
             'NOMA, K=3, Energy-Efficient PA', ...
             'NOMA, K=2, Energy-Efficient PA', ...
             'Conventional OMA, Energy-Efficient PA');
