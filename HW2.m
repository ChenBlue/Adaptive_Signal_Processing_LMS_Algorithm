N = 2000; % length of sequence
M = 11; %taps
sim_times = 200; % simulation times
tau = 7; % delays
mu = 0.075; % step size
W = 2.9; %Channel parameter
%d = [zeros(sim_times, tau) double(randn(sim_times,N)>0)]; % input signal
d = double(randn(sim_times,N)>0); % input signal
%d = [zeros(sim_times, M) double(randn(sim_times,N)>0)]; % input signal
%d = [double(randn(sim_times,N+tau)>0)]; % input signal
d_ext = [double(randn(sim_times, tau)>0) d];
%d_ext = [zeros(sim_times, tau) d];
%d_ext = [double(randn(sim_times, tau)>0) d];
n = [1 2 3];
x = 1:1:2000;
a = 1;

%%
% Different W
W_test = 2.9:0.2:3.5;
mse = zeros(length(W_test),N); % mean squared error
y = zeros(1,N); % output of equalizer
e = zeros(1,N); % error

for k = 1:length(W_test)
    h_tmp = 1/2*(1+cos(2*pi*(n-2)/W_test(k))); 
    h = [0 h_tmp]; % raised cosine channel
    for j = 1:sim_times
        %u_tmp = conv(d(j,tau+1:end),h);
        %u_tmp = filter(h, a, d(j,tau+1:end));
        u_tmp = conv(d_ext(j,:),h);
        u = [zeros(1,M-tau-1) u_tmp];
        %u_tmp = conv(d(j,:),h);
        %u = [zeros(1,M-1) u_tmp];
        %u = [zeros(1,M-1) u_tmp(1:N)] + randn(1,M+N-1)*sqrt(0.001);
        %u = [u_tmp(1:N) zeros(1,M-1)] + randn(1,M+N-1)*sqrt(0.001);
        w = randn(1,M);

        for i=1:N
            U = fliplr(u(i:i+M-1));
            y(i) = U * transpose(w);
            e(i) = d_ext(j,i)-y(i);
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
ylabel('Mean squared error');
title('Learning curve');

%%
figure,plot(x, mse(1,:),x, mse(2,:),x, mse(3,:), x, mse(4,:));
set(gca, 'YScale', 'log')
xlabel('Number of adaptation cyckes, n');
ylabel('Mean squared error');
title('Learning curve for different W');
legend({'W = 2.9','W = 3.1', 'W = 3.3', 'W = 3.5'},'Location','northeast')

%%
%%%%% Different mu
mu_test = [0.0075 0.025 0.075 0.1 0.15];
mse_mu = zeros(length(mu_test),N); % mean squared error
y = zeros(1,N); % output of equalizer
e = zeros(1,N); % error
h_tmp = 1/2*(1+cos(2*pi*(n-2)/W));  % W=2.9
h = [0 h_tmp]; % raised cosine channel

for k = 1:length(mu_test)
    for j = 1:sim_times
        %u_tmp = conv(d(j,tau+1:end),h);
        %u_tmp = filter(h, a, d(j,tau+1:end));
        %u_tmp = conv(d(j,:),h);
        u_tmp = conv(d_ext(j,:),h);
        u = [zeros(1,M-tau-1) u_tmp];
        %u = [zeros(1,M-1) u_tmp];
        %u = [zeros(1,M-1) u_tmp(1:N)] + randn(1,M+N-1)*sqrt(0.001);
        %u = [u_tmp(1:N) zeros(1,M-1)] + randn(1,M+N-1)*sqrt(0.001);
        w = randn(1,M);

        for i=1:N
            U = fliplr(u(i:i+M-1));
            y(i) = U * transpose(w);
            e(i) = d_ext(j,i)-y(i);
            w = w + mu_test(k)*e(i)*U;
        end
        mse_mu(k,:) = mse_mu(k,:) + e.^2;
    end
end
mse_mu = mse_mu/sim_times;

%%
figure,plot(x, mse_mu(1,:),x, mse_mu(2,:),x, mse_mu(3,:), x, mse_mu(4,:), x, mse_mu(5,:));
set(gca, 'YScale', 'log')
xlabel('Number of adaptation cyckes, n');
ylabel('Mean squarred error');
title('Learning curve for different step size £g');
legend({'£g = 0.0075','£g = 0.025', '£g = 0.075', '£g = 0.1' , '£g = 0.15'},'Location','northeast')

%%
%%%%% Different delay
tau_test = [4 6 7 8 9];
mse_delay = zeros(length(tau_test),N); % mean squared error
y = zeros(1,N); % output of equalizer
e = zeros(1,N); % error
h_tmp = 1/2*(1+cos(2*pi*(n-2)/W(1)));  % W(1)=2.9
h = [0 h_tmp]; % raised cosine channel

for k = 1:length(tau_test)
    d_test = [zeros(sim_times, tau_test(k)) d];
    for j = 1:sim_times
        %u_tmp = conv(d(j,tau+1:end),h);
        %u_tmp = filter(h, a, d(j,tau+1:end));
        %u_tmp = conv(d(j,:),h);
        u_tmp = conv(d_ext(j,:),h);
        u = [zeros(1,M-tau_test(k)-1) u_tmp];
        %u = [zeros(1,M-1) u_tmp];
        %u = [zeros(1,M-1) u_tmp(1:N)] + randn(1,M+N-1)*sqrt(0.001);
        %u = [u_tmp(1:N) zeros(1,M-1)] + randn(1,M+N-1)*sqrt(0.001);
        w = randn(1,M);

        for i=1:N
            U = fliplr(u(i:i+M-1));
            y(i) = U * transpose(w);
            e(i) = d_test(j,i)-y(i);
            w = w + mu*e(i)*U;
        end
        mse_delay(k,:) = mse_delay(k,:) + e.^2;
    end
end
mse_delay = mse_delay/sim_times;

%%
figure,plot(x, mse_delay(1,:),x, mse_delay(2,:),x, mse_delay(3,:), x, mse_delay(4,:), x, mse_delay(5,:));
set(gca, 'YScale', 'log')
xlabel('Number of adaptation cyckes, n');
ylabel('Mean squared error');
title('Learning curve for different dalay £n');
legend({'£n = 4','£n = 6', '£n = 7', '£n = 8' , '£n = 9'},'Location','northeast')

%%
mse_noise = zeros(length(W_test),N); % mean squared error
y = zeros(1,N); % output of equalizer
e = zeros(1,N); % error

for k = 1:length(W_test)
    h_tmp = 1/2*(1+cos(2*pi*(n-2)/W_test(k))); 
    h = [0 h_tmp]; % raised cosine channel
    for j = 1:sim_times
        %u_tmp = conv(d(j,tau+1:end),h);
        %u_tmp = filter(h, a, d(j,tau+1:end));
        %u_tmp = conv(d(j,:),h);
        u_tmp = conv(d_ext(j,:),h);
        u = [zeros(1,M-tau-1) u_tmp(1:N+tau)] + randn(1,M+N-1)*sqrt(0.001);
        %u = [zeros(1,M-1) u_tmp];
        %u = [zeros(1,M-1) u_tmp(1:N)] + randn(1,M+N-1)*sqrt(0.001);
        %u = [u_tmp(1:N) zeros(1,M-1)] + randn(1,M+N-1)*sqrt(0.001);
        w = randn(1,M);

        for i=1:N
            U = fliplr(u(i:i+M-1));
            y(i) = U * transpose(w);
            e(i) = d_ext(j,i)-y(i);
            w = w + mu*e(i)*U;
        end
        mse_noise(k,:) = mse_noise(k,:) + e.^2;
    end
end
mse_noise = mse_noise/sim_times;

%%
figure,plot(x, mse_noise(1,:),x, mse_noise(2,:),x, mse_noise(3,:), x, mse_noise(4,:));
set(gca, 'YScale', 'log')
xlabel('Number of adaptation cyckes, n');
ylabel('Mean squared error');
title('Learning curve for different W with noise');
legend({'W = 2.9','W = 3.1', 'W = 3.3', 'W = 3.5'},'Location','northeast')

%%
noise_test = [sqrt(0.00001) sqrt(0.0001) sqrt(0.001) sqrt(0.01) sqrt(0.1)];
mse_n = zeros(length(noise_test),N); % mean squared error
y = zeros(1,N); % output of equalizer
e = zeros(1,N); % error
h_tmp = 1/2*(1+cos(2*pi*(n-2)/W_test(1))); 
h = [0 h_tmp]; % raised cosine channel

for k = 1:length(noise_test)
    for j = 1:sim_times
        %u_tmp = conv(d(j,tau+1:end),h);
        %u_tmp = filter(h, a, d(j,tau+1:end));
        %u_tmp = conv(d(j,:),h);
        u_tmp = conv(d_ext(j,:),h);
        u = [zeros(1,M-tau-1) u_tmp(1:N+tau)] + randn(1,M+N-1)*noise_test(k);
        %u = [zeros(1,M-1) u_tmp];
        %u = [zeros(1,M-1) u_tmp(1:N)] + randn(1,M+N-1)*sqrt(0.001);
        %u = [u_tmp(1:N) zeros(1,M-1)] + randn(1,M+N-1)*sqrt(0.001);
        w = randn(1,M);

        for i=1:N
            U = fliplr(u(i:i+M-1));
            y(i) = U * transpose(w);
            e(i) = d_ext(j,i)-y(i);
            w = w + mu*e(i)*U;
        end
        mse_n(k,:) = mse_n(k,:) + e.^2;
    end
end
mse_n = mse_n/sim_times;

%%
figure,plot(x, mse_n(1,:),x, mse_n(2,:),x, mse_n(3,:), x, mse_n(4,:), x, mse_n(5,:));
set(gca, 'YScale', 'log')
xlabel('Number of adaptation cyckes, n');
ylabel('Mean squared error');
title('Learning curve for different noise power');
legend({'Noise power = 0.00001','Noise power = 0.0001', 'Noise power = 0.001', 'Noise power = 0.01', 'Noise power = 0.1'},'Location','northeast')
