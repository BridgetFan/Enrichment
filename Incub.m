function rhs = Incub(~,X,par,varargin)
% ODE system for the droplet/well-mixed system

rhs = X;
n = length(X);

P = X(1:4:n);
N = X(2:4:n);
R = X(3:4:n);
Q = X(4:4:n);

S = par.rp*P+N+R;

a = par.a;
q = par.q;
theta_Q = par.theta_Q;

f_QS = a*Q.^q./(theta_Q^q+Q.^q);

g4 = par.g4;
gamma = par.gamma;

rhs(1:4:n) = par.g1*P.*(par.S_max-S)/par.S_max;
rhs(2:4:n) = par.g2*N.*(par.S_max-S)/par.S_max;
rhs(3:4:n) = par.g3*R.*(par.S_max-S)/par.S_max-f_QS(:).*R;

if nargin>3
    p = varargin{1};
    rhs(4:4:end) = par.QS_in*g4*P+(1-par.QS_in)*g4*P'*p-gamma*Q;
else
    rhs(4:4:end) = g4*P-gamma*Q;
end

end
