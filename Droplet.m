%% 
% Simulation for the Producer Enrichment Using Micro-droplet System
% For any questions, please contact 
% Gaoyang Fan at bridget.gfan@gmail.com

clear;
%% Parameters
Drop = 1;             % Droplet or Well-mixed
T = 12;               % Single Simulation Cycle Time (unit: hr)
k = 1;                % Total No. of Simulation Cycles
d = 50;               % Dilution Factor

par.g4 = 0.5;         % AHL Production Rate Low/High: 0.5/10;

% IPTG Dependent Producer Growth Rates
if par.g4 > 1
    par.g1 = 0.5404;
    par.rp = 1.6;
else
    par.g1 = 0.6078;
    par.rp = 1;
end

par.g2 = 0.6452;      % Non-Producer Growth Rate
par.g3 = 0.5131;      % Reciever Growth Rate

if Drop == 1
    ini_V = .5;       % Initial Vol (unit: mL)
    par.S_max = 150;  % Droplet Carrying Capacity (unit: #/drop)
    pd = 33.5*1e-9;   % Droplet Size (unit: mL)
else
    ini_V = 5;        % Initial Vol (unit: mL)
    par.S_max = 400*sum(ini_V); % Tube Carrying Capacity (unit: M)
end

% QS signal Induced Killing
par.q = 2;            % Hill Coeff.
par.a = 0.65;         % Max Kill Rate

par.QS_in = 0.5;      % Percentage of Signal Produced Stay In-droplet
par.gamma=0;          % Signal Degradation

theta = 1;            % Threshold of Killing per mL

if Drop == 1
    par.theta_Q = theta * 10^6 * pd;
else
    par.theta_Q = theta * ini_V;
end

sorter = 1;           % If Sorter is Added
sorter_amount = 1;    % Total # of Droplet Sorted
sorter_theta = 0.95;  % Sorting Threshold

% Initial Pop. Density for Each Strain (unit: M/mL)
des_P = 0.06;
des_N = 6;
des_R = 300;

des = [des_P, des_N, des_R];

V_ratio = [1, 1, 2];            % Volume Ratio P:N:R
V = V_ratio*ini_V/sum(V_ratio); % Volume of P, N, R (unit: mL)

ini = des.*V;                   % Initial Pop. Count for P, N, R (unit: M)
R_add = ini(3);                 % # of Recievers Spiked-in Each Round

if Drop == 1
    drops = ini_V/pd/10^6;      % Number of Droplet
    m = 50;                     % Max. # of Cells/Strain/Droplet Considered
    N = length(des);            % # of Strains
    M_all = zeros((m+1)^N,N);   % Initialization of the I.C. Matrix

    % Construct the I.C. Matrix
    for i=1:N
        v_i = repmat(0:m,1,(m+1)^(i-1));
        v_ij = repmat(v_i,(m+1)^(N-i),1);
        M_all(:,i) = v_ij(:);
    end
    M_all(:,4) = 0;             % Zero Out AHL I.C. for All Droplets
    p = zeros(size(M_all,1),1); % Initialize the Probability Matrix
end

t_all = [];
dt = 0.01; % Plot Time Interval

%% Numerical Solution
Nt = T/dt+1;
Y_all = zeros(Nt*k,3);
Q_all = zeros(Nt*k,1);

for i = 1:k
    if Drop == 1
        % Calculation of the Prob. for each I.C.
        lambda = ini/drops;
        pois1  = poisspdf(M_all(:,1),lambda(1));
        pois2  = poisspdf(M_all(:,2),lambda(2));
        pois3  = poisspdf(M_all(:,3),lambda(3));
        p_all = pois1.* pois2.*pois3;

        % Truncate I.C. Matrix by Prob.<eps
        [p_descend, ind] = sort(p_all,'descend');
        ind_sig = sum(p_descend>=1/(drops*10^7));
        p = p_descend(1:ind_sig);
        M = M_all(ind(1:ind_sig),:);
        p = p/sum(p);
    end
    % Solving ODE for Droplet or Well-mixed System
    if Drop == 1
        ini_all = reshape(M',[],1);
        [t,Y] = ode45(@(t,X) Incub(t,X,par,p), (0:dt:T)+(i-1)*T, ini_all);
        Y_P = Y(:,1:4:end);
        Y_N = Y(:,2:4:end);
        Y_R = Y(:,3:4:end);
        Y_Q = Y(:,4:4:end)*p*drops;
        Y_count=[Y_P*p,Y_N*p,Y_R*p]*drops;
    else
        [t,Y] = ode15s(@(t,X) Incub(t,X,par),(0:dt:T)+(i-1)*T, [ini 0]);
        Y_count=Y(:,1:3);
        Y_Q=Y(:,4);
    end
    Y_all((i-1)*Nt+1:i*Nt,:)=Y_count;
    Q_all((i-1)*Nt+1:i*Nt,:)=Y_Q;
    t_all((i-1)*Nt+1:i*Nt)=t;

    ini = Y_count(end,:)/d + [0, 0, R_add];

    if sorter == 1
        Y_sorted_all=zeros(k,3);
        p_In=p;

        tspan = 0:dt:T;
        [~,G] = ode45(@(t,y) GFP(t,y,tspan,Y_R), tspan, 0.*Y_R(end,:));
        [Y_ordered,ranks] = sort(G(end,:),2,'descend');
        ind_Out = cumsum(p(ranks))<sorter_theta;
        p_In(ranks(ind_Out)) = 0;
        sort_R = Y_ordered(sum(ind_Out)+1);

        % Sorting out exactly threshold percentage of droplets
        if sum(ind_Out)
            p_n_theta=p_In(ranks(sum(ind_Out)+1))-(sorter_theta-sum(p(ranks(ind_Out))));
            p_In(ranks(sum(ind_Out)+1))=p_n_theta;
        end
        Y_sorted=[Y_P(end,:)*p_In,Y_N(end,:)*p_In,Y_R(end,:)*p_In]*sorter_amount;
        Y_sorted_all(i,:)=Y_sorted;
        ini = Y_sorted/d + [0,0,R_add];
    end
end
