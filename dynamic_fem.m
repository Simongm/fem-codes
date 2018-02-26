function dynamic_fem


% %% initialize movie:
% close all;figure
% myVideo = VideoWriter('dynamic.avi'); % create object
% myVideo.FrameRate = 30;  % frame rate Default 30
% open(myVideo);


%% define constants:
G=3e9; %shear modulus
nu=0.35; %poisson ratio
rho=3500; %density
E=2*G*(1+nu);
vs=(G/rho)^0.5;
vp=((2*G*(1-nu))/(rho*(1-2*nu)))^0.5;
lambda=(2*G*nu)/(1-2*nu);
% or lambda2=vp^2*rho-2*G;


%% construct the grid
nz=200;nx=200; %nb of elements
x=linspace(0,0.3,nx);
z=linspace(0,0.3,nz);
dx=abs(x(2)-x(1));
dz=abs(z(2)-z(1));
[xx,zz]=meshgrid(x,z);
dt=1/(vp*sqrt(1/dx^2+1/dz^2)); %timestep according to courant stability criterion
% nt=200; %nb of iterations
% t=0:nt;t=t*dt;
% tmax=dt*nt;
tmax=0.001;
t=linspace(0,tmax,round(tmax/dt));nt=length(t);
%% initializes the vectors
vx_old=zeros(size(xx));
vz_old=zeros(size(xx));
txx_old=zeros(size(xx));
tzz_old=zeros(size(xx));
txz_old=zeros(size(xx));
vx_new=vx_old;
vz_new=vz_old;
txx_new=txx_old;
tzz_new=tzz_old;
txz_new=txz_old;


%% source terms
% amplitude of gaussian pulse for point source
f0=50/(tmax); 
t0=2/f0;
amp=1*exp(-f0^2*(t-t0).^2); amp=amp-amp(1);


%% prestress to start dynamic rupt
size_patch=0.008; % size half patch in m
nb_pt_patch=round(size_patch/dx);
nxs=round(nx/2);nzs=round(nz/2);Dc=0.085e-5;
rup=zeros(1,nx);friction=zeros(1,nx);
prestress(1:nx)=1.7*10^6;prestress(nxs-nb_pt_patch:nxs+nb_pt_patch)=1.8001*10^6;peak(1:nz)=1.8e6;
u_interface_old=zeros(1,nx);
mu_d=0.8; %actually coresponds to tau_f/tau_p ratio
tau_f=mu_d*peak(1);
cc=jet(nt);
figure; 
% hold on;
stress_init=zeros(size(xx))+prestress(1,1);%stress_init(nzs,nxs-nb_pt_patch:nxs+nb_pt_patch)=1.8001*10^6;

for i=1:nt
    loading_rate=0.2e9;
    
%     txz_old=txz_old+loading_rate*dt;
    [vx_new,vz_new,txx_new,tzz_new,txz_new]=tnew_virieux(xx,zz,vx_old,vz_old,txx_old,tzz_old,txz_old,vx_new,vz_new,txx_new,tzz_new,txz_new,lambda,G,dt,rho);

    % src term for kinematic ruputre
%     vx_new(round(length(x)/2),round(length(z)/2))=amp(i);
%     tyz_new(100,:)=amp(i);


    %% or spontaneous rupture
    u_interface_new=u_interface_old+vx_new(nzs,:)*dt;
     for j=1:nx %go along interface
         if (rup(j)==0); %check elements that have not broken yet
         if (abs(txz_new(nzs,j)+prestress(j))>peak(j));   
           rup(j)=1;
         end;
         end;
     end
     for j=1:nx;
        if (rup(j)==1);
            friction(j)=1-u_interface_new(j)*(1-mu_d)/Dc; 
            if (friction(j)<mu_d) 
                friction(j)=mu_d; 
            end
            txz_new(nzs,j)=peak(j)*friction(j)-prestress(j);

%            if (vx_new(nzs,j)<0); %for a self-healing rupture
%                rup(j)=0;
%            end;
        end
     end
     txz_new(nzs,txz_new(nzs,:)<(peak(1)*mu_d-prestress(1,1)))=peak(1)*mu_d-prestress(1,1);
     txz_new(nzs,txz_new(nzs,:)>(peak(1)-prestress(1,1)))=peak(1)-prestress(1,1);
     u_interface_old=u_interface_new;
     
    %to increase the loading
    peak=peak-loading_rate*dt;     
    tau_f=tau_f-loading_rate*dt;
    mu_d=(tau_f/peak(1));
     %% plots
%for stress at interface-------------
    if rem(i*dt,1*dt)==0
%         plot(xx(nzs,:),txz_new(nzs,:)+(i-1)*loading_rate*dt,'color',cc(round(i*nt/500),:))
%         ylim([1e5,3e6])
        stress2plt=txz_new+stress_init+(i-1)*loading_rate*dt;
        surf(xx,zz,stress2plt,'edgecolor','none')
        title([num2str(t(i)),'s'])
%         view(2)
%         colorbar()
%     zlim([tau_f, tau_p])
        drawnow()
%         writeVideo(myVideo, getframe); 
%         hold on
    end
%----------------------------------

    
    
    %%update vectors
    vx_old=vx_new;
    vz_old=vz_new;
    txx_old=txx_new;
    tzz_old=tzz_new;
    txz_old=txz_new;
    
    
    
end

% close(myVideo)

function [vx_new,vz_new,txx_new,tzz_new,txz_new]=tnew_virieux(xx,zz,vx_old,vz_old,txx_old,tzz_old,txz_old,vx_new,vz_new,txx_new,tzz_new,txz_new,lambda,G,dt,rho)

%for velocities we use the old stress variables
dtxxdxx=diff(txx_old,1,2)./diff(xx,1,2);dtxxdxx=dtxxdxx(1:end-1,:);dtxzdzz=diff(txz_old,1,1)./diff(zz,1,1);dtxzdzz=dtxzdzz(:,2:end);
vx_new(1:end-1,2:end)=vx_old(1:end-1,2:end)+dt/rho*(dtxxdxx+dtxzdzz);

dtxzdxx=diff(txz_old,1,2)./diff(xx,1,2);dtxzdxx=dtxzdxx(2:end,:);dtzzdzz=diff(tzz_old,1,1)./diff(zz,1,1);dtzzdzz=dtzzdzz(:,1:end-1);
vz_new(2:end,1:end-1)=vz_old(2:end,1:end-1)*dt/rho*(dtxzdxx+dtzzdzz);


% %% sponge boundaries  
%  for j=1:nz;for i=1:nn;sx(j,i)=sx(j,i)*(ax1*i*dx+cx1);end;end;
%  for j=1:nz;for i=nx-nn:nx;sx(j,i)=sx(j,i)*(ax2*i*dx+cx2);end;end;
%  for j=nz-nn:nz;for i=1:nx;sz(j,i)=sz(j,i)*(az*j*dx+cz);end;end;
 
 
%for stresses we use the newly calculated velocities
dvxdxx= diff(vx_new,1,2)./diff(xx,1,2);dvxdxx=dvxdxx(1:end-1,:);dvzdzz=diff(vz_new,1,1)./diff(zz,1,1);dvzdzz=dvzdzz(:,1:end-1);
txx_new(1:end-1,1:end-1)=txx_old(1:end-1,1:end-1)+dt*((lambda+2*G)*dvxdxx+lambda*dvzdzz);
tzz_new(1:end-1,1:end-1)=tzz_old(1:end-1,1:end-1)+dt*(lambda*dvxdxx+(lambda+2*G)*dvzdzz);

dvxdzz=diff(vx_new,1,1)./diff(zz,1,1);dvxdzz=dvxdzz(:,2:end);dvzdxx=diff(vz_new,1,2)./diff(xx,1,2);dvzdxx=dvzdxx(2:end,:);
txz_new(2:end,2:end)=txz_old(2:end,2:end)+dt*G*(dvxdzz+dvzdxx);

%% impose BC

% vx_new(1,:)=0;vx_new(:,1)=0;vx_new(end,:)=0;vx_new(:,end)=0;
% vz_new(1,:)=0;vz_new(:,1)=0;vz_new(end,:)=0;vz_new(:,end)=0;

txx_new(1,:)=0;txx_new(:,1)=0;txx_new(end,:)=0;txx_new(:,end)=0;
tzz_new(1,:)=0;tzz_new(:,1)=0;tzz_new(end,:)=0;tzz_new(:,end)=0;
txz_new(1,:)=0;txz_new(:,1)=0;txz_new(end,:)=0;txz_new(:,end)=0;

