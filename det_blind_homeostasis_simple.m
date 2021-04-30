% parameter values from monolix fit & literature
r = 0.06;
beta = 10^(-4)*10^(-6);
dI = 1;
p = 10^4;
c = 23;
K = (2.8);%*10^(6);
f = 1e-4;
a = 5.7e-5;
g = 0.05;

% splining the viral load
times = VLP1186(1,:):1e-3:VLP1186(end,:);
slope1 = (VLP1186(2,2)-VLP1186(1,2))/(VLP1186(2,1)-VLP1186(1,1));
slope2 = 0;
VIRUS = spline(VLP1186(:,1).',[slope1 VLP1186(:,2).' slope2],times);
% plot(times,VIRUS)
% hold on
% plot(VLP1186(:,1),VLP1186(:,2),'*')

% CD4 Interpolates
CD4 = zeros(size(times));
% linear fit initial decline of cd4
% moved the pre-point to same as virus
%CD4_part1 = polyfit(CD4P1186(1:3,1),CD4P1186(1:3,2),1);
CD4_part1 = spline(CD4P1186(1:3,1).',10.^CD4P1186(1:3,2).',times(1:find(times>=11,1)));%polyval(CD4_part1,times(1:find(times>=11,1)));
% plot(times(1:find(times>=11,1)),CD4_part1)
% hold on
% plot(CD4P1186(1:3,1),10.^CD4P1186(1:3,2),'*')

CD4(1:size(CD4_part1,2)) = CD4_part1;

% logisitic growth of recovery
[t,C] = ode45(@(t, C) logistic(t,C,K,r),...
    [1.54 11.1],(10^0.312));
CD4(find(times>=11,1)+1) = C(end);

for i=(find(times>=11,1)+2):size(times,2)
    [t,C] = ode45(@(t, C) logistic(t,C,K,r),...
        [times(i-1) times(i)],CD4(i-1));
    CD4(i)=C(end);
end

CD4((size(CD4_part1,2)+1):end) = 10.^CD4((size(CD4_part1,2)+1):end);
% plot(times,CD4)
% hold on
% plot(CD4P1186(:,1),10.^CD4P1186(:,2),'*')
% 
% K=10^K*10^6;
% CD4 = 10^6*CD4;
% % Deterministic latent reservoir and infected cells compartment
% L = zeros(2,size(times,2));
% I = zeros(2,size(times,2));
% % LI(2,1)=1;
% for i=2:size(times,2)
%     [t,C] = ode45(@(t, C) LandI(t,C,K,r,f,g,beta,a,dI,CD4(i-1),10^VIRUS(i-1)*10^3),...
%         [times(i-1) times(i)],[LI(1,i-1) LI(2,i-1)]);
%     LI(1,i)=C(end,1);
%     LI(2,i)=C(end,2);
% end

K=10^K*10^6;
CD4 = 10^6*CD4;
% Deterministic latent reservoir and infected cells compartment
L = zeros(4,size(times,2));
I = zeros(4,size(times,2));
% LI(2,1)=1;
r = 0.06*[0 0.25 1 2];
for h = 1:4
for i = 2:size(times,2)
    [t,C] = ode45(@(t, C) LandI(t,C,K,r(h),f,g,beta,a,dI,CD4(i-1),10^VIRUS(i-1)*10^3),...
        [times(i-1) times(i)],[L(h,i-1) I(h,i-1)]);
    L(h,i)=C(end,1);
    I(h,i)=C(end,2);
end
end



