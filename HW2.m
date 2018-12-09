N = 2000; % length of sequence
M = 11; %taps
sim_times = 200; % simulation times
tau = 7; % delays
mu = 0.075; % step size
d = double(randn(sim_times,N)>0); % input signal
case_num = 4;
n = [1 2 3];
x = 1:1:2000;
%%
%W = 2.9; %Channel parameter
W = 2.9:0.2:3.5;
U = zeros(1,M); 
mse = zeros(case_num,N); % mean squared error
y = zeros(1,N); % output of equalizer
e = zeros(1,N); % error
%u_rev = fliplr(u);

for k = 1:4
    h_tmp = 1/2*(1+cos(2*pi*(n-2)/W(k))); 
    h = [0 h_tmp]; % raised cosine channel
    for j = 1:sim_times
        u_tmp = conv(d(j,:),h);
        u = [zeros(1,M-1) u_tmp(1:N)];
        %u = [zeros(1,M-1) u_tmp(1:N)] + randn(1,M+N-1)*sqrt(0.001);
        w = zeros(1,M);

        for i=tau+1:N
            U = fliplr(u(i:i+M-1));
            y(i) = U * transpose(w);
            e(i) = d(j,i-tau)-y(i);
            w = w + mu*e(i)*U;
        end
        mse(k,:) = mse(k,:) + e.^2;
    end
end

mse = mse/sim_times;
%%
figure,plot(x, mse(1,:));
set(gca, 'YScale', 'log')
xlabel('Number of adaptation cyckes, n');
ylabel('Mean squarred error');
title('Learning curve');

%%
figure,plot(x, mse(1,:),x, mse(2,:),x, mse(3,:), x, mse(4,:));
set(gca, 'YScale', 'log')
xlabel('Number of adaptation cyckes, n');
ylabel('Mean squarred error');
title('Learning curve');
legend({'W = 2.9','W = 3.1', 'W = 3.3', 'W = 3.5'},'Location','northeast')