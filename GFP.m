function dydt = GFP(t,y,t_f,f)
% ODE for the GFP production

R = interp1(t_f,f,t);       % Interpolate the data set (t_f,f) at time t
p_gfp = 1;                  % GFP Production Rate
gamma = 0;                  % GFP Degradate Rate
dydt = p_gfp.*R' - gamma*y; % Evaluate ODE at time t
end
