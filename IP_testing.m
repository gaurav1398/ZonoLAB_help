close all; clear all; clc
%% INVERTED PENDULUM
% Using the Hybrid Zonotope approximation of Sin/Cos as the dynamics of a unicycle:
%     x1_{k+1} = x1_k + dt*x2_k
%     x2_{k+1} = x2_k + dt*[(g/l)*sin(x1_k) + (1/m*l^2)*u_k] 
% (the dt term could be absorbed into u_k.)

% Initial State
initialAngle = 12; % degrees
initialAngle = initialAngle * (2 * pi / 360); % convert to radians
initAngleDelta = 2; % degrees
initAngleDelta = initAngleDelta * (2 * pi / 360); % convert to radians
initVelocity = 0; % m/s
Z0 = zono([initAngleDelta 0;0 0.01],[initialAngle;initVelocity]);
U = zono(0,0);

% System parameters
m = 5;
l = 1;
g = 9.81;
mu = 0.4;
k0 = 1 - mu;

resFreq = -g/l;
J = 1/(m*l^2);

%% Trig Functions as Hybrid Zonotopes
n_points = 3; % Number of points (should be an odd number)

sinZono = makeSinX(n_points,[-initialAngle - initAngleDelta,initialAngle + initAngleDelta]);
%% Successor Set method
N = 20;
Z = evolveUnicycleMemZono(sinZono,U,Z0,N,resFreq,J,k0);

%% Plotting

for i = 1:N
    plot(Z{i},{'x','v'},'b',0.1)
    title('Memory Approach')
end

% drawnow;
% for k = 1:N
%     plot(Z_inter, {sprintf('x_%d',k),sprintf('v_%d',k)}, selectColor(k), 0.6);
%     plot(Z{k}, 'all', selectColor(k), 0.2);
%     drawnow;
% end
%% Helper Functions
function [Z] = evolveUnicycleMemZono(sinZono,U,Z0,N,resFreq,J,k0 ...
    )
dt = 1;

%     Zx = {memZono(Z0_x,sprintf('x_0'))};                                      
%     Zv = {memZono(Z0_v,sprintf('v_0'))};
%     % Give the state the correct dimension labels
%     Zx{1}.dimKeys = {'x'};
%     Zv{1}.dimKeys = {'v'};
%     Z{1} = combine(Zx{1},Zv{1});

Z = {memZono(Z0,sprintf('x_1'))};
Z{1}.dimKeys = {'x','v'};
Z_all = Z{1};
X_F = zono(0.25*eye(2),ones(2,1));
X_F = memZono(X_F,sprintf('x_%d',N));


switch 'transform&merge' %'Combine' 'transform&merge' 'linSysTrans&merge' 
    case 'Combine'
    
    for k = 1:N-1
        sinMZ = memZono(sinZono,sprintf('sin_k%i',k));
        % sin theta updates the v state
        sinMZ.dimKeys = {'x','v'};

        UMZ = memZono(U,sprintf('U_%i',k-1));
        UMZ.dimKeys = {'v'};
        %  x_{k+1}  = x_{k} + v_{k}
        %  v_{k+1}  = v_{k} + (g/l)*sin(x_{k}) + U_{k})
        
        dX = Z{k}({'v'});
        Z{k+1} = combine(dX.copy('x'),dt*Z{k});

        dV1 = resFreq*merge(sinMZ,Z{k}({'x'}),sprintf('sin_v_%i',k));
        dV2 = combine(UMZ,dV1.copy('v'));
        Z{k+1} = combine(mu*Z{k}.copy('v'),dV2);
    end


    case 'transform&merge'
        for k=2:N
            sinMZ = memZono(sinZono,sprintf('sin_k%i',k));
            % sin theta updates the v state
            sinMZ.dimKeys = {'x','v'};
            
            UMZ = memZono(U,sprintf('U_%i',k));
            UMZ.dimKeys = {'v'};

            % individual functions
            w1 = Z{k-1};

            w2 = transform(Z{k-1}({'v'}),[],dt,{},{'x'});

            w3 = transform(merge(sinMZ,Z{k-1}({'x'}),sprintf('sin_v_%i',k)),[],resFreq,{},{'v'});

            w4 = UMZ;
            w4.dimKeys = {'v'};

            w5 = transform(Z{k-1}({'v'}),[],k0,{},{'v'});

            Z{k} = combine(w1,w2);
            Z{k} = dt*combine(combine(w3,w4),w5);
        end


    case 'linSysTrans&merge'
        A = [0 1;resFreq mu];
        B = [0;J];
        for k = 1:N-1
            % Current Input
            U_{k} = memZono(U,sprintf('u_%d',k));

            % Step Update
            newDims = {sprintf('x_%d',k+1),sprintf('v_%d',k+1)};
            Z{k+1} = Z{k}.transform(U_{k}.transform([],B,{},newDims),A,{},newDims); %<== transform has affine A x + B

            % Save Data
            Z_all = Z_all.merge(U_{k});
            Z_all = Z_all.merge(Z{k+1});
        end
        Z_inter = Z_all.merge(X_F,'terminal_cons'); % <--- intersect common dimensions
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