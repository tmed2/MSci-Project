%This script implements the reaction-diffusion model of cancer using the
%finite element modelling methods in PDE toolbox.
%
%The system being modlled is that of a organ-on-a-chip and a cuboidal
%geometry has been chosen to model the system
%
%Strange instabilities resulting in negative concentrations; these are due
%to "steep" initial conditions
%

datetime('now')
clear variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Parameters
global D_c
global D_O2
global alpha
global K_mic
global v_max
global O2_0

D_c = 200; %um^2/day (can be up to 2000)
D_O2 = 86.4e4; %um^2/day
alpha = 2; %per day
K_mic = 0.12; %umol/gram
v_max = 82; %umol/gram/day
O2_0 = 3;
c_coeffs = [D_c;D_O2];

%Geometic parameters
global xmax
global ymax
global zmax
global tmax

xmax = 1000;
ymax = 1000;
zmax = 100;
tmax = 21;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic %start timer

%create a pde object, and define the geometry of the problem
model = createpde(2);
[x,y,z] = meshgrid(0:xmax:xmax,0:ymax:ymax,0:zmax:zmax); %distances in um
t = 0:1:tmax; %time in days
x = x(:);
y = y(:);
z = z(:);
K = convhull(x,y,z);
nodes = [x';y';z'];
elements = K';
geometryFromMesh(model,nodes,elements);

% figure;
pdegplot(model,'FaceLabels','on','FaceAlpha',0.5); %plots the geometry for checking
ax = gca;
ax.XLim = [0,xmax];
ax.YLim = [0,ymax];
ax.ZLim = [0,zmax];
xlabel('x / um');
ylabel('y / um');
zlabel('z / um');


%No flux of cells outside of the boundary
applyBoundaryCondition(model,'neumann',...
                             'Face',1:6,...
                             'g',0,...
                             'q',0);

applyBoundaryCondition(model,'dirichlet',...
                             'Face',[3,5,6],...
                             'h',[0,0;0,1],...
                             'r',[0,O2_0]);

specifyCoefficients(model,'m',0,...
                          'd',[1;1],...
                          'c',c_coeffs,...
                          'a',0,...
                          'f',@source_func);

setInitialConditions(model,@init_func);
generateMesh(model);
result = solvepde(model,t);

toc %stop timer

tic %start timer

%prints out images of plots at multiple times
[X,Y,Z] = meshgrid(0:10:xmax,0:10:ymax,0:5:zmax);
U1 = interpolateSolution(result,X,Y,Z,1,1:length(t));
U1 = squeeze(U1);

for j=1:length(t)
t_fig = t(j);
fig_title = 'Cancer Cell Concentration at ';
fig_title = fig_title + string(t_fig) + ' days';
U1_t = U1(:,j);
U1_t = reshape(U1_t,size(X));

h = figure('visible', 'off');
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
c = colorbar;
c.Label.String = 'w';
caxis([0 1]);

view(-45,45);
axis equal;
t_fig_str = strcat('w_t',num2str(t_fig));
print(h,'-bestfit',t_fig_str,'-dpdf');
print(h,t_fig_str,'-dpng');
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
U1_t = U2(:,j);
U1_t = reshape(U1_t,size(X));

h = figure('visible', 'off');
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
c = colorbar;
c.Label.String = 'O2';
caxis([0, 3*O2_0]);

view(-45,45);
axis equal;
t_fig_str = strcat('O2_t',num2str(t_fig));
print(h,'-bestfit',t_fig_str,'-dpdf');
print(h,t_fig_str,'-dpng');
close
end

toc %stop timer

%make a single plot
% [X,Y,Z] = meshgrid(0:5:400,0:5:400,0:2:100);
% V = interpolateSolution(result,X,Y,Z,1:length(t));
% V = V(:,50);
% V = reshape(V,size(X));
% figure;
% colormap jet;
% contourslice(X,Y,Z,V,[10,100,200,300,390],[10,100,200,300,390],[25,50,75]);
% ax1 = gca;
% ax1.XLim = [0,400];
% ax1.YLim = [0,400];
% ax1.ZLim = [0,100];
% xlabel('x / um');
% ylabel('y / um');
% zlabel('z / um');
% title('Contours of Dimentionless Cancer Cell Concentration at 10 days')
% c = colorbar;
% c.Label.String = 'u';
% view(-45,45);
% axis equal;

%animation needs work
% sol_anim = figure;
% V = interpolateSolution(result,X,Y,Z,1:length(t));
% xlabel('x / um');
% ylabel('y / um');
% zlabel('z / um');
% colormap jet;
% c = colorbar;
% c.Label.String = 'u';
% caxis([0,1]);
% view(0,90);
% axis tight manual
% ax = gca;
% ax.NextPlot = 'replaceChildren';
% ax.XLim = [0,400];
% ax.YLim = [0,400];
% ax.ZLim = [0,100];
% title('Evolution of Dimentionless Concentration over 100 days')
% loops = length(t);
% F(loops) = struct('cdata',[],'colormap',[]);
% for j = 1:loops
%     V_t = V(:,j);
%     V_t = reshape(V_t,size(X));
%     contourslice(X,Y,Z,V_t,[100,200,300,390],[100,200,300,390],[10,25,50,75,90]);
%     F(j) = getframe;
%     cla(ax)
% end
% movie(sol_anim,F,1);


function S = source_func(region,state)
%This defines the quadratic loss term for logistic proliferation
global alpha;
global K_mic;
global v_max
%

N=2;
nr = length(region.x);
S = zeros(N,nr);

O2_rate_prop = ((state.u(2,:)./(K_mic + state.u(2,:))));

S(1,:) = alpha*state.u(1,:).*(O2_rate_prop - state.u(1,:));
S(2,:) = -v_max*O2_rate_prop.*state.u(1,:);
end


function u_init = init_func(locations)
%This function defines the initial conditions of the problem, a gaussian
%has been chosen to approximate a point source for the cancer cells and a
%constant has been chosen for the inital O2 concentration
global O2_0;
global xmax
global ymax
global zmax
%

c_0 = 0.02;
sigma = 10;
sigmaz = sigma*((zmax*zmax)/(xmax*ymax));
Nclusters = 5;

M = length(locations.x);
u_init = zeros(2,M);
for i=1:Nclusters
    %random start position of cell clusters
    centre = (0.2 + 0.6*rand(1,3)).*[xmax,ymax,zmax];
    centre(3) = zmax/2;
    u1_x = sqrt(2*pi)*sigma*normpdf(locations.x,centre(1),sigma);
    u1_y = sqrt(2*pi)*sigma*normpdf(locations.y,centre(2),sigma);
    u1_z = sqrt(2*pi)*sigmaz*normpdf(locations.z,centre(3),sigmaz);
    u_init(1,:) = u_init(1,:) + c_0.*u1_x.*u1_y.*u1_z;
end

u_init(2,:) = O2_0; %uniform inital o2 concentration
end