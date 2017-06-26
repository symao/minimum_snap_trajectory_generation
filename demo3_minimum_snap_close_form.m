function demo3_minimum_snap_close_form()
clear,clc;

%% condition
waypts = [0,0;
          1,2;
          2,-1;
          4,8;
          5,2]';
v0 = [0,0];
a0 = [0,0];
v1 = [0,0];
a1 = [0,0];
T = 5;
ts = arrangeT(waypts,T);
n_order = 5;

%% trajectory plan
polys_x = minimum_snap_single_axis_close_form(waypts(1,:),ts,n_order,v0(1),a0(1),v1(1),a1(1));
polys_y = minimum_snap_single_axis_close_form(waypts(2,:),ts,n_order,v0(2),a0(2),v1(2),a1(2));

%% result show
figure(1)
plot(waypts(1,:),waypts(2,:),'*r');hold on;
plot(waypts(1,:),waypts(2,:),'b--');
title('minimum snap trajectory');
color = ['grc'];
for i=1:size(polys_x,2)
    tt = ts(i):0.01:ts(i+1);
    xx = polys_vals(polys_x,ts,tt,0);
    yy = polys_vals(polys_y,ts,tt,0);
    plot(xx,yy,color(mod(i,3)+1));
end

end

function polys = minimum_snap_single_axis_close_form(wayp,ts,n_order,v0,a0,v1,a1)
n_coef = n_order+1;
n_poly = length(wayp)-1;
% compute Q
Q_all = [];
for i=1:n_poly
    Q_all = blkdiag(Q_all,computeQ(n_order,3,ts(i),ts(i+1)));
end

% compute Tk   Tk(i,j) = ts(i)^(j-1)
tk = zeros(n_poly+1,n_coef);
for i = 1:n_coef
    tk(:,i) = ts(:).^(i-1);
end

% compute A (n_continuous*2*n_poly) * (n_coef*n_poly)
n_continuous = 3;  % 1:p  2:pv  3:pva  4:pvaj  5:pvajs
A = zeros(n_continuous*2*n_poly,n_coef*n_poly);
for i = 1:n_poly
    for j = 1:n_continuous
        for k = j:n_coef
            if k==j
                t1 = 1;
                t2 = 1;
            else %k>j
                t1 = tk(i,k-j+1);
                t2 = tk(i+1,k-j+1);
            end
            A(n_continuous*2*(i-1)+j,n_coef*(i-1)+k) = prod(k-j+1:k-1)*t1;
            A(n_continuous*2*(i-1)+n_continuous+j,n_coef*(i-1)+k) = prod(k-j+1:k-1)*t2;
        end
    end
end

% compute M
M = zeros(n_poly*2*n_continuous,n_continuous*(n_poly+1));
for i = 1:n_poly*2
    j = floor(i/2)+1;
    rbeg = n_continuous*(i-1)+1;
    cbeg = n_continuous*(j-1)+1;
    M(rbeg:rbeg+n_continuous-1,cbeg:cbeg+n_continuous-1) = eye(n_continuous);
end

% compute C
num_d = n_continuous*(n_poly+1);
C = eye(num_d);
df = [wayp,v0,a0,v1,a1]';% fix all pos(n_poly+1) + start va(2) + end va(2) 
fix_idx = [1:3:num_d,2,3,num_d-1,num_d];
free_idx = setdiff(1:num_d,fix_idx);
C = [C(:,fix_idx) C(:,free_idx)];

AiMC = inv(A)*M*C;
R = AiMC'*Q_all*AiMC;

n_fix = length(fix_idx);
Rff = R(1:n_fix,1:n_fix);
Rpp = R(n_fix+1:end,n_fix+1:end);
Rfp = R(1:n_fix,n_fix+1:end);
Rpf = R(n_fix+1:end,1:n_fix);

dp = -inv(Rpp)*Rfp'*df;

p = AiMC*[df;dp];

polys = reshape(p,n_coef,n_poly);

end

