close all;
%% INVERTED PENDULUM
% Using the Hybrid Zonotope approximation of Sin/Cos as the dynamics of a unicycle:
%     x1_{k+1} = x1_k + x2_k
%     x2_{k+1} = x2_k + (g/l)*sin(x1_k) + u_k 
% (the dt term could be absorbed into u_k.)

% Initial State
Z0_x = zono(0.001,0.5);
Z0_v = zono(0.001,0);
U = zono(0.001,0);

% System parameters
m = 1;
l = 1;
g = 10;

%% Trig Functions as Hybrid Zonotopes
x_max = pi/4;
x_min = -x_max;
n_points = 3; % Number of points (should be an odd number)

sinZono = makeSinX(n_points,[-pi pi]);
%% Successor Set method
N = 10;
Z = evolveUnicycleMemZono(sinZono,U,Z0_x,Z0_v,N);

%% Plotting
% merge all time steps together
% to do so we need to relabel the dimKeys
% Z{i}({'x','y'}) ensures that the dimKeys are in the expected order for relabeling
% i = 1;
% liftZ = copy(Z{i}({'x','v'}),{sprintf('x_%i',i),sprintf('v_%i',i)});
% for i = 2:N
%     liftZ = merge(liftZ, copy(Z{i}({'x','v'}),{sprintf('x_%i',i),sprintf('v_%i',i)}));
% end

% imagine that new data is received that informs the y value at time step 3
% Z_newdata = memZono(zono(0.005,0.025),'z_pin');
% Z_newdata.dimKeys = 'v_3';
% liftZ_newdata = merge(liftZ,Z_newdata,'pin_merge');


for i = 1:N
    plot(Z{i},{'x','v'},'b',0.1)
    title('Memory Approach')
end


%% Helper Functions
function [Z] = evolveUnicycleMemZono(sinZono,U,Z0_x,Z0_v,N)
    dt = 0.1;

    Zx = {memZono(Z0_x,sprintf('x_0'))};                                      
    Zv = {memZono(Z0_v,sprintf('v_0'))};
    % Give the state the correct dimension labels
    Zx{1}.dimKeys = {'x'};
    Zv{1}.dimKeys = {'v'};
    Z{1} = combine(Zx{1},Zv{1});
    for k = 2:N
        %===== These need to be inside the loop so that each usage of sin/cos/U is memory independent
        sinMZ = memZono(sinZono,sprintf('sin_k%i',k));
        % sin theta updates the v state
        sinMZ.dimKeys = {'x','v'};
      
        UMZ = memZono(U,sprintf('U_%i',k-1));
        UMZ.dimKeys = {'v'};
    
        % does x_{k+1}  = x_{k} + v_{k}
        %      v_{k+1}  = v_{k} + g/l*sin(x_{k}) + U_{k}
        Zx{k} = combine(Zx{k-1},dt*Zv{k-1});

        % does X_{k+1} = ... + v*sin(theta)
        % only use Z{k-1}({'theta'})
        %dv = minSum(10*merge(sinMZ,Z{k-1}({'x'}),sprintf('sin_v_%i',k)),UMZ{k-1});
        temp = 10*merge(sinMZ,Zx{k-1},sprintf('sin_x%i',k));
        Zv{k} = combine(dt*combine(temp({'v'}),UMZ),Zv{k-1});   
        %Zv{k} = combine(Zv{k},Zv{k-1});
        Z{k} = combine(Zx{k},Zv{k});
    end
end


function sinx = makeSinX(n,bd)

th = linspace(bd(1), bd(2), 2*n+1);
xi = th;
yi = sin(th);

V = [xi; yi];
nc = 2*n+1;
nb = nc-1;
naux = nc;
c = zeros(2,1);
Gc = V;

M = zeros(nc,nb);
for i = 1:nb
    M(i,i) = 1;
    M(i+1,i) = 1;
end

c = c+0.5*Gc*ones(nc,1);
Gc = 0.5*Gc;
Gc = [Gc,zeros(2,naux)];
Gb = zeros(2,nb);

Ac = [0.5*ones(1,nc) zeros(1,naux)]; %(1)
Ab = zeros(1,nb);
b = [1-0.5*nc];

Ac = [Ac; 0.5*eye(nc) 0.5*eye(nc)]; %(2)
Ab = [Ab; -0.5*M];
b = [b; 0.5*M*ones(nb,1)-ones(nc,1)];

Ac = [Ac; zeros(1,nc+naux)]; %(3)
Ab = [Ab; 0.5*ones(1,nb)];
b = [b; 1-0.5*nb];

sinx = hybZono(Gc,Gb,c,Ac,Ab,b);

end