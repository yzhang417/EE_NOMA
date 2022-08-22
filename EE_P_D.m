clc;
clear;
close all;

%circulation
LOOP = 500;
%Power of Noise
sigmadBm = -70;
sigma = 10^(sigmadBm/10)*10^(-3); 
%path loss factor
a = 3.0;
%dm
%dm = [ 60 50 40; 70 60 40; 60 55 50; 80 80 80;];
dm = [ 60 50 40; 70 55 40; 60 55 50; 80 80 80;];
%M
M_array = 3*ones(1,4);
%Q
QoS = 1.0; 
%Ptot
PtotmindBm = 5;
PtotmaxdBm = 35;
dPdBm = 1;
Ptotrange = PtotmindBm:dPdBm:PtotmaxdBm;
Ptotrange2 = 10.^(Ptotrange/10)*10^(-3);
Pc_dBm = 30;
Pc = 10^(Pc_dBm/10)*10^(-3);
%Result Array
len = length(PtotmindBm:dPdBm:PtotmaxdBm);
Rs = zeros(LOOP,len);
Rs_TDMA = zeros(LOOP,len);
AverageRs = zeros(length(M_array),len);
AverageRs_TDMA = zeros(3,len);

%% NOMA
for index=1:1:length(M_array)
    index = index
    % Number of licit user
    M = M_array(index);
    %Qm
    Q = QoS*ones(1,M);
    for loop = 1:1:LOOP
       %% Channel Information
        h = 1/sqrt(2)*randn(1,M)+1/sqrt(2)*randn(1,M)*1j;
        h = h.*conj(h);
        nh = sort(h)./(1+dm(index,:).^a);                                       
        P = zeros(1,M);
        B = (2.^(Q)-1)./nh;
        for m = M : -1 : 1
            P(m) = B(m)*(nh(m)*sum(P(m+1:M)) + sigma); 
        end
        Pmin = sum(P);
        PmindBm = 10*log10(Pmin*10^3);  
        compter=0;
        while PmindBm > PtotmindBm
            compter = compter+1;
            if compter>20 
                index = index;
                loop = loop
                break;
            end          
            index;
            loop = loop;
            % Channel gain of licit user 
            h = 1/sqrt(2)*randn(1,M)+1/sqrt(2)*randn(1,M)*1j;
            h = h.*conj(h);
            nh = sort(h)./(1+dm(index,:).^a);                                        
            P = zeros(1,M);
            B = (2.^(Q)-1)./nh;
            for m = M : -1 : 1
                P(m) = B(m)*(nh(m)*sum(P(m+1:M)) + sigma); 
            end
            Pmin = sum(P);
            PmindBm = 10*log10(Pmin*10^3);
        end
        Ptotindex = 0;
        for PtotdBm = PtotmindBm : dPdBm : PtotmaxdBm
            Ptotindex = Ptotindex+1;
            Ptot = 10.^(PtotdBm/10)*10^(-3);  
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
            Rs(loop,Ptotindex) = sum(Rbm);
        end       
    end
    AverageRs(index,:) = nanmean(Rs)./(Ptotrange2+Pc);
end

%% Data EE
EE = AverageRs;
for i=1:1:4
    [val pos] = max(EE(i,:));
    EE(i,pos:end) = val;
end

%% Plot
figure;
Ptotrange = PtotmindBm : dPdBm : PtotmaxdBm;
plot(Ptotrange,EE(1,:),'-k'); hold on;
plot(Ptotrange,AverageRs(1,:),'--k'); hold on;
plot(Ptotrange,EE(2,:),'-b'); hold on;
plot(Ptotrange,AverageRs(2,:),'--b'); hold on;
plot(Ptotrange,EE(3,:),'-r'); hold on;
plot(Ptotrange,AverageRs(3,:),'--r'); hold on; 
plot(Ptotrange,EE(4,:),'-g'); hold on;
plot(Ptotrange,AverageRs(4,:),'--g'); hold on;
grid on;
xlabel('$P$ (dBm)');ylabel('Average EE (bits/Joule/Hz)');
legend(['Distance: ',num2str(dm(1,:)),', Energy-Efficient PA'],['Distance: ',num2str(dm(1,:)),', Full Power Consumption'], ...
            ['Distance: ',num2str(dm(2,:)),', Energy-Efficient PA'],['Distance: ',num2str(dm(2,:)),', Full Power Consumption'], ...
            ['Distance: ',num2str(dm(3,:)),', Energy-Efficient PA'],['Distance: ',num2str(dm(3,:)),', Full Power Consumption'], ...
            ['Distance: ',num2str(dm(3,:)),', Energy-Efficient PA'],['Distance: ',num2str(dm(3,:)),', Full Power Consumption'] )
