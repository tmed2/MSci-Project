%This script implements the reaction-diffusion model of cancer using the
%finite element modelling methods in PDE toolbox.
%
%The system being modelled is that of a organ-on-a-chip and a cuboidal
%geometry has been chosen to model the system with appropoarte boundary
%conditions
%

clear variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Parameters
global D_c
global D_O2
global D_g
global alpha
global K_o2
global v_o2
global O2_0
global K_g
global v_g
global g_0

%this will print the values to the console
                                        %default
D_c = 100%um^2/day (can be up to 2000)  %100
D_O2 = 8.64e7 %um^2/day                 %8.64e7
D_g = 3.648e6 %um^2/day                 %3.648e6
alpha = 2 %per day                      %2
K_o2 = 4.3 %umol                        %4.3
v_o2 = 100 %umol/day                    %100
O2_0 = 260 %umol                        %260
K_g = 0.5 %mmol                         %0.5
v_g = 7 %mmol/day                       %7
g_0 = 5 %mmol                           %5
c_coeffs = [D_c;0;0;0;D_c;0;0;0;D_c;D_O2;0;0;0;D_O2;0;0;0;D_O2;D_g;0;0;0;D_g;0;0;0;D_g];
 
%Geometic parameters
global xmax
global ymax
global zmax
global tmax
global allcentres
global Nclusters
global tstep

xmax = 1000 %um                         %1000
ymax = 1000 %um                         %1000
zmax = 100  %um                         %100
tmax = 30   %days                       %51
tstep = 3   %days                       %3

%must initalise random centres outside of initial function scope as it is
%called multiple times, and randomness within will confuse the solver
Nclusters = 3;
%allcentres = (0.25 + 0.5*rand(Nclusters,3)).*[xmax,ymax,zmax]

%we use the same pregenerated random centres for consistency
allcentres = [553.7475, 578.2829, 39.4964; 597.1873, 640.3779, 65.8522; 415.7227, 253.0553, 43.7209]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
datetime('now')
tic %start timer

%create a pde object, and define the geometry of the problem
model = createpde(3);
[x,y,z] = meshgrid(0:xmax:xmax,0:ymax:ymax,0:zmax:zmax); %distances in um
t = 0:tstep:tmax; %time in days
x = x(:);
y = y(:);
z = z(:);
K = convhull(x,y,z);
nodes = [x';y';z'];
elements = K';
geometryFromMesh(model,nodes,elements);

% figure;
% pdegplot(model,'FaceLabels','on','FaceAlpha',0.5); %plots the geometry for checking
% ax = gca;
% ax.XLim = [0,xmax];
% ax.YLim = [0,ymax];
% ax.ZLim = [0,zmax];
% xlabel('x / um');
% ylabel('y / um');
% zlabel('z / um');
% view(-45,45);

%No flux of cells and nutrients through closed boundaries
applyBoundaryCondition(model,'neumann',...
                             'Face',1:6,...
                             'g',0,...
                             'q',0);

%oxygen and glucose through sides
applyBoundaryCondition(model,'dirichlet',...
                             'Face',[3,5],...
                             'h',[0,0,0;0,1,0;0,0,1],...
                             'r',[0,O2_0,g_0]);

%oxygen through top
applyBoundaryCondition(model,'dirichlet',...
                             'Face',6,...
                             'h',[0,0,0;0,1,0;0,0,0],...
                             'r',[0,O2_0,0]);

%specify the coeffients to make our system of equations
specifyCoefficients(model,'m',0,...
                          'd',[1;1;1],...
                          'c',c_coeffs,...
                          'a',0,...
                          'f',@source_func);

setInitialConditions(model,@init_func);
mesh = generateMesh(model,'Hmax',30);
result = solvepde(model,t);

toc %stop timer

tic %start timer

%prints out images of plots at given times
[X,Y,Z] = meshgrid(0:xmax/100:xmax,0:ymax/100:ymax,0:zmax/10:zmax);
U1 = interpolateSolution(result,X,Y,Z,1,1:length(t));
U1 = squeeze(U1);

for j=1:length(t)
t_fig = t(j);
fig_title = 'Cancer Cell Concentration at ';
fig_title = fig_title + string(t_fig) + ' days';
U1_t = U1(:,j);
U1_t = reshape(U1_t,size(X));

h1 = figure('visible', 'off');
contourslice(X,...
             Y,...
             Z,...
             U1_t,...
             0:xmax/10:xmax,...
             0:ymax/10:ymax,...
             [0.25*zmax,0.5*zmax,0.75*zmax]);
ax1 = gca;
ax1.XLim = [0,xmax];
ax1.YLim = [0,ymax];
ax1.ZLim = [0,zmax];
xlabel('x / um');
ylabel('y / um');
zlabel('z / um');
title(fig_title)
colormap jet;
c1 = colorbar;
c1.Label.String = 'u';
caxis([0 1]);

view(-45,45);
axis equal;
t_fig_str = strcat('u_t',num2str(t_fig));
print(h1,'-bestfit',t_fig_str,'-dpdf');
print(h1,t_fig_str,'-dpng');
close
end

toc %stop timer
tic %start timer

U2 = interpolateSolution(result,X,Y,Z,2,1:length(t));
U2 = squeeze(U2);

for j=1:length(t)
t_fig = t(j);
fig_title = 'Oxygen Concentration at ';
fig_title = fig_title + string(t_fig) + ' days';
U2_t = U2(:,j);
U2_t = reshape(U2_t,size(X));

h2 = figure('visible', 'off');
contourslice(X,...
             Y,...
             Z,...
             U2_t,...
             0:xmax/10:xmax,...
             0:ymax/10:ymax,...
             [0.25*zmax,0.5*zmax,0.75*zmax]);

ax2 = gca;
ax2.XLim = [0,xmax];
ax2.YLim = [0,ymax];
ax2.ZLim = [0,zmax];
xlabel('x / um');
ylabel('y / um');
zlabel('z / um');
title(fig_title)
colormap jet;
c2 = colorbar;
c2.Label.String = 'O2 / umol';
caxis([0, 1.1*O2_0]);



view(-45,45);
axis equal;
t_fig_str = strcat('O2_t',num2str(t_fig));
print(h2,'-bestfit',t_fig_str,'-dpdf');
print(h2,t_fig_str,'-dpng');
close
end

toc %stop timer

tic %start timer

U3 = interpolateSolution(result,X,Y,Z,2,1:length(t));
U3 = squeeze(U3);

for j=1:length(t)
t_fig = t(j);
fig_title = 'Glucose Concentration at ';
fig_title = fig_title + string(t_fig) + ' days';
U3_t = U3(:,j);
U3_t = reshape(U3_t,size(X));

h3 = figure('visible', 'off');
contourslice(X,...
             Y,...
             Z,...
             U3_t,...
             0:xmax/10:xmax,...
             0:ymax/10:ymax,...
             [0.25*zmax,0.5*zmax,0.75*zmax]);

ax3 = gca;
ax3.XLim = [0,xmax];
ax3.YLim = [0,ymax];
ax3.ZLim = [0,zmax];
xlabel('x / um');
ylabel('y / um');
zlabel('z / um');
title(fig_title)
colormap jet;
c3 = colorbar;
c3.Label.String = 'Glucose / mmol';
caxis([0, 1.1*g_0]);



view(-45,45);
axis equal;
t_fig_str = strcat('g_t',num2str(t_fig));
print(h3,'-bestfit',t_fig_str,'-dpdf');
print(h3,t_fig_str,'-dpng');
close
end

toc %stop timer

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Source and inital condition functions defined here

function S = source_func(region,state)
%This defines the quadratic loss term for logistic proliferation with
%nutrition
global alpha
global K_o2
global v_o2
global K_g
global v_g
%

N=3;
nr = length(region.x);
S = zeros(N,nr);

%cf Michaelis Menten enzyme dynamics
O2_rate_prop = ((state.u(2,:)./(K_o2 + state.u(2,:))));
glu_rate_prop = ((state.u(3,:)./(K_g + state.u(3,:))));

S(1,:) = alpha*state.u(1,:).*(O2_rate_prop.*glu_rate_prop - state.u(1,:));
S(2,:) = -v_o2*O2_rate_prop.*state.u(1,:);
S(3,:) = -v_g*glu_rate_prop.*state.u(1,:);
end


function u_init = init_func(locations)
%This function defines the initial conditions of the problem
global O2_0
global g_0
global xmax
global ymax
global zmax
global allcentres
global Nclusters
%

c_0 = 0.01;
sigma = 10;
sigmaz = sigma*((zmax*zmax)/(xmax*ymax));

N=3;
M = length(locations.x);
u_init = zeros(N,M);
for i=1:Nclusters
    %random start position of cell clusters
    centre = allcentres(i,:);
    %centre(3) = zmax/2;
    u1_x = sqrt(2*pi)*sigma*normpdf(locations.x,centre(1),sigma);
    u1_y = sqrt(2*pi)*sigma*normpdf(locations.y,centre(2),sigma);
    u1_z = sqrt(2*pi)*sigmaz*normpdf(locations.z,centre(3),sigmaz);
    u_init(1,:) = u_init(1,:) + c_0.*u1_x.*u1_y.*u1_z;
end

%uniform inital oxygen and glucose concentration
u_init(2,:) = O2_0;
u_init(3,:) = g_0;
end