%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  THzsource_FDTD_2Dv1.m
% 
% 2D FDTD suimulation with progagting THz source.
% Fields are Ex,Ez, Hy , with source and geometry independent of y
%   [ making (Ex,Ez,Hy) and (Ez,Hx,Hz) independent fields ]
% Two sources in simulation, both progagating in the z direction, 
% at a vleovity that can be tuned <c or >c 
% the two sources are constrained to x=2mm, and x = 3mm
%
% MJC, 16/7/2014



% Clearing variables in memory and Matlab command screen
clear all;
%clc;

% Grid Dimension in x (xdim) and z (ydim) directions
xdim=1000;
zdim=800;

% Total no of time steps
time_tot=800;

% Courant stability factor
S=1/(2^0.5);

% Parameters of free space (permittivity and permeability and speed of
% light) are all not 1 and are given real values
epsilon0=(1/(36*pi))*1e-9;
mu0=4*pi*1e-7;
c=3e+8;

% Spatial grid step length (spatial grid step= 1 micron and can be changed)
delta=5e-6;
% Temporal grid step obtained using Courant condition
deltat=S*delta/c;

% Position of the source (center of the field domain)

%for lmn = 1:10;

%varname = ['xsource',num2str(lmn)];
%['xsource',num2str(lmn)],'i']=
xsource=2.5e-3;
    xsourcei = round(xsource/delta);
xsource2=2.55e-3;
    xsource2i = round(xsource2/delta);
xsource3=2.6e-3;
    xsource3i = round(xsource3/delta);
xsource4=2.65e-3;
    xsource4i = round(xsource4/delta);

    
zi = 1:zdim;
z = zi*delta;
xi = 1:xdim;
x = xi*delta;

% Initialization of field matrices
Ex=zeros(xdim,zdim);
Ez=zeros(xdim,zdim);
Hy=zeros(xdim,zdim);

% Initialization of permittivity and permeability matrices
epsilon=epsilon0*ones(xdim,zdim);
%epsilon(1:xdim/2,1:zdim)=1.0*epsilon0;
epsilon(250:750,1:400)=10.0*epsilon0;
%epsilon(580:xdim,1:zdim)=10.0*epsilon0;


% Left part of domain containing glass of n=1.5
mu=mu0*ones(xdim,zdim);

% Initializing electric and magnetic conductivity matrices
sigma=4e-4*ones(xdim,zdim)*0;
sigma_star=4e-4*ones(xdim,zdim)*0;

%Choice of nature of source
gaussian=1;
sine=0;
% The user can give a frequency of his choice for sinusoidal (if sine=1 above) waves in Hz 
frequency=1.5e+13;
%Choose any one as 1 and other as 0. 


figure(1)
% Update loop begins
for n=1:1:time_tot

    
    % Setting spatial update limits
    n1=1;
    n2=xdim-1;
    n11=1;
    n21=zdim-1;
        
    %Vector update instead of for-loop for Hy fields
    Hy(n1:n2,n11:n21) = ((mu(n1:n2,n11:n21)-0.5*deltat*sigma_star(n1:n2,n11:n21))./(mu(n1:n2,n11:n21)+0.5*deltat*sigma_star(n1:n2,n11:n21))).*Hy(n1:n2,n11:n21) ... 
        +(deltat/delta)./(mu(n1:n2,n11:n21)+0.5*deltat*sigma_star(n1:n2,n11:n21)).* ...  
               ( Ez(n1+1:n2+1,n11:n21)-Ez(n1:n2,n11:n21)    -Ex(n1:n2,n11+1:n21+1)+Ex(n1:n2,n11:n21) );
   
           
    
    %Vector update instead of for-loop for Ex and Ez fields
    Ez(n1+1:n2+1,n11+1:n21+1)=((epsilon(n1+1:n2+1,n11+1:n21+1)-0.5*deltat*sigma(n1+1:n2+1,n11+1:n21+1))./(epsilon(n1+1:n2+1,n11+1:n21+1)+0.5*deltat*sigma(n1+1:n2+1,n11+1:n21+1))).*Ez(n1+1:n2+1,n11+1:n21+1) ...
        +(deltat/delta)*(Hy(n1+1:n2+1,n11+1:n21+1)-Hy(n1:n2,n11+1:n21+1))./(epsilon(n1+1:n2+1,n11+1:n21+1)+0.5*deltat*sigma(n1+1:n2+1,n11+1:n21+1));
    
    Ex(n1+1:n2+1,n11+1:n21+1)=((epsilon(n1+1:n2+1,n11+1:n21+1)-0.5*deltat*sigma(n1+1:n2+1,n11+1:n21+1))./(epsilon(n1+1:n2+1,n11+1:n21+1)+0.5*deltat*sigma(n1+1:n2+1,n11+1:n21+1))).*Ex(n1+1:n2+1,n11+1:n21+1) ...
        -(deltat/delta)*(Hy(n1+1:n2+1,n11+1:n21+1)-Hy(n1+1:n2+1,n11:n21))./(epsilon(n1+1:n2+1,n11+1:n21+1)+0.5*deltat*sigma(n1+1:n2+1,n11+1:n21+1));
    
    
    % Smoothening ends of the line source
    Ez(1:xdim,1)=Ez(1:xdim,2);
    Ez(1:xdim,zdim)=Ez(1:xdim,zdim-1);
    
    % Source conditions
    %if sine
    if sine==1
        tstart=1;
        N_lambda=c/(frequency*delta);
        Ez(xsource,1:zdim)=sin(((2*pi*(c/(delta*N_lambda))*(n-tstart)*deltat)));
    end
    %if gaussian
    if gaussian==1
        %%%% source type
        %  1 =  single guassian
        %  2 =  guassian, opposing pair
        %  3 =  gaussian derivative, opposing pair
        %  4 =  double (time-spaced) gaussian ; opposing pair
        %  5 =  
        %  6 =  
        sourceType = 5;
        
        
        sigmat = 2000e-15;
        sigmaz = 4000e-6;
        theta = 30*pi/180;
        z0 = 1000e-6;
        
        ng = 2;
        betaS = 1/ng/sin(theta);
        if sourceType==1 % single guassian
            Ez(xsourcei,zi) = exp(-(n*deltat -600e-15- (z-z0)*sin(theta)/c*ng ).^2/sigmat^2).* ...
                exp(-(z-z0).^2*cos(theta)^2/sigmaz^2);
            
        elseif sourceType==2   %guassian, opposing pair
            Ez(xsourcei,zi) = exp(-(n*deltat -600e-15- (z-z0)*sin(theta)/c*ng ).^2/sigmat^2).* ...
                exp(-(z-z0).^2*cos(theta)^2/sigmaz^2) ;
          
        elseif sourceType==3  %gaussian derivative, opposing pair
            
            ttt = (n*deltat -600e-15- (z-z0)*sin(theta)/c*ng )*1e15;
            Ez(xsourcei,zi) = ttt.*exp(-(n*deltat -600e-15- (z-z0)*sin(theta)/c*ng ).^2/sigmat^2).* ...
                exp(-(z-z0).^2*cos(theta)^2/sigmaz^2) ;
            
         
        elseif sourceType==4  %double (time-spaced) gaussian ; opposing pair
            
            sigmat = 2000e-15;
            sigmaz = 2000e-6;
            theta = 30*pi/180;
            z0 = 1500e-6;
            
            Ez(xsourcei,zi) = exp(-(n*deltat -600e-15- (z-z0)*sin(theta)/c*ng ).^2/sigmat^2).* ...
                exp(-(z-z0).^2*cos(theta)^2/sigmaz^2) + ...
                exp(-(n*deltat -2000e-15- (z-z0)*sin(theta)/c*ng ).^2/sigmat^2).* ...
                exp(-(z-z0).^2*cos(theta)^2/sigmaz^2);
            
            
        elseif sourceType==5  %Optical rectification of asymetric guassian
            % assuming a laser pulse of form
            %     Alaser = [ erf(t/sRise) + 1 ]*exp(-t/sFall)
            % and optical rectificfayion giving THz field profile as time
            % derivatitave of laser envelope, we then have
            %    ETHz(t) = exp(-t/sFall)*[ 2/sRise/sqrt(2)exp(-t^2/sRise^2) - ...
            %                 1/sFall*erf(t/sRise) - 1/sFall ]
             
            t0 = 600e-15;
            t = n*deltat - t0 - (z-z0)*sin(theta)/c*ng ;
            sR = 0.5e-12;
            sF = 15e-12;
      
            
             Ez(xsourcei,zi) =  exp(-t/sF ).* ...
              ( 2/sR/sqrt(2).*exp(-t.^2/sR^2) - ...
                             1/sF*erf(t/sR) - 1/sF ).* ...
                             exp(-(z-z0).^2*cos(theta)^2/sigmaz^2);
            
             
            Ez(xsource2i,zi) = exp(-t/sF ).* ...
              ( 2/sR/sqrt(2).*exp(-(t-5e-13).^2/sR^2) - ...
                             1/sF*erf(t/sR) - 1/sF ).* ...
                         exp(-(z-z0).^2*cos(theta)^2/sigmaz^2);
            
            Ez(xsource3i,zi) = exp(-t/sF ).* ...
              ( 2/sR/sqrt(2).*exp(-(t-10e-13).^2/sR^2) - ...
                             1/sF*erf(t/sR) - 1/sF ).* ...
                         exp(-(z-z0).^2*cos(theta)^2/sigmaz^2);
            Ez(xsource4i,zi) = exp(-t/sF ).* ...
              ( 2/sR/sqrt(2).*exp(-(t-15e-13).^2/sR^2) - ...
                             1/sF*erf(t/sR) - 1/sF ).* ...
                         exp(-(z-z0).^2*cos(theta)^2/sigmaz^2);
            
        end;
        
    end
    if mod(n-1,42) == 0,
        
        zframe = ( 3e8*(n*deltat-5433e-15) + 1800*1e-6 )*1e6;
        
        
        
        if 1 ==1
            figure(1)
            pcolor((delta*(1e+6:1e+6:zdim*1e+6))',delta*(1e+6:1e+6:xdim*1e+6),Ez);
            shading interp;
            drawnow
            
        end
        
        if 1==2, %main plots
            

            
            %Movie type colour scaled image plot of Ez
            subplot(311)
            imagesc(delta*(1e+6:1e+6:xdim*1e+6),(delta*(1e+6:1e+6:zdim*1e+6))',Ez'); %,[-1,2]);
            colorbar; %('north');
            line(xsource*delta*[1 1]*1e6,[0,500])
            title(['\fontsize{14}E_z at time = ',num2str(round(n*deltat*1e+15)),' fs']);
            axis image
            %axis([200,400,0,500])
            %zframe = ( 3e8*(n*deltat-5433e-15) + 2400*1e-6 )*1e6
            axis([2000-100,3000+100,0+zframe,1500+zframe])
            xlabel('x (in um)','FontSize',14);
            ylabel('z (in um)','FontSize',14);
            %set(gca,'FontSize',11);
            grid on
            hold on
            plot(2500,zframe+600,'o')
            hold off
            
            subplot(312)
            imagesc(delta*(1e+6:1e+6:xdim*1e+6),(delta*(1e+6:1e+6:zdim*1e+6))',Ex'); %,[-1,1]);
            colorbar; %('north');
            line(xsource*delta*[1 1]*1e6,[0,500])
            axis image
            %axis([200,400,0,500])
            % zframe = ( 3e8*(n*deltat-5433e-15) + 1800*1e-6 )*1e6
            axis([2000-100,3000+100,0+zframe,1500+zframe])
            title(['\fontsize{14}E_x at time = ',num2str(round(n*deltat*1e+15)),' fs']);
            xlabel('x (in um)','FontSize',14);
            ylabel('z (in um)','FontSize',14);
            %set(gca,'FontSize',11);
            grid on
            hold on
            plot(2500,zframe+600,'o')
            hold off
            
            subplot(313)
            imagesc(delta*(1e+6:1e+6:xdim*1e+6),(delta*(1e+6:1e+6:zdim*1e+6))',c*mu0*Hy');
            colorbar; %('north');
            line(xsource*delta*[1 1]*1e6,[0,500])
            axis image
            % zframe = ( 3e8*(n*deltat-5433e-15) + 2400*1e-6 )*1e6
            axis([2000-100,3000+100,0+zframe,1500+zframe])
            title(['\fontsize{14}cB_y at time = ',num2str(round(n*deltat*1e+15)),' fs']);
            xlabel('x (in um)','FontSize',14);
            ylabel('z (in um)','FontSize',14);
            %set(gca,'FontSize',11);
            grid on
            hold on
            plot(2500,zframe+600,'o')
            hold off
            drawnow
            
        end
        %getframe;
        
        %         pagefill('p');
        %     print('-dpdf',['THzsource_FDTD_' num2str(n*deltat*1e15,'%3.0f') 'fs']);
        
        
        if 1==2 %addiational force plots
            %%%%%%%%%%%%%%%%
            
            Fx = Ex + c*mu0*Hy;
            Fz = Ez;
            
            xxi = xi(x>2000e-6 & x< 3000e-6);
            zzi = zi(z>500e-6 & z< 3000e-6);
            
            xxi2 = min(xxi):10:max(xxi);
            zzi2 = min(zzi):10:max(zzi);
            
            xx = x(xxi2)*1e6;
            zz = z(zzi2)*1e6;
            
            
            [XX,ZZ] = ndgrid(xx,zz);
            
            figure(2)
            subplot(211)
            imagesc(delta*(1e+6:1e+6:xdim*1e+6),(delta*(1e+6:1e+6:zdim*1e+6))',Hy');
            hold on
            quiver(XX',ZZ',Ex(xxi2,zzi2)',Ez(xxi2,zzi2)');
            axis image
            axis([1800,3200,[0,1500]+zframe])
            title(['\fontsize{14} {Ex,Ez} and H_y (colourmap) at time = ',num2str(round(n*deltat*1e+15)),' fs']);
            xlabel('x (in um)','FontSize',14);
            ylabel('z (in um)','FontSize',14);
            set(gca,'FontSize',14);
            grid on
            colorbar
            
            subplot(223)
            imagesc(delta*(1e+6:1e+6:xdim*1e+6),(delta*(1e+6:1e+6:zdim*1e+6))',Fz');
            hold on
            quiver(XX',ZZ',Fx(xxi2,zzi2)',Fz(xxi2,zzi2)');
            axis image
            axis([1800,3200,[-500,500]+zframe+500])
            title(['\fontsize{14} {Fx,Fz} and F_z (colourmap), \theta = ',...
                num2str(theta*180/pi,'%3.1f'),'^0,', ...
                ' \beta_S = ' num2str(betaS,'%3.2f')]);
            xlabel('x (in um)','FontSize',14);
            ylabel('z (in um)','FontSize',14);
            set(gca,'FontSize',14);
            grid on
            colorbar
            
            subplot(224)
            imagesc(delta*(1e+6:1e+6:xdim*1e+6),(delta*(1e+6:1e+6:zdim*1e+6))',Fx');
            hold on
            quiver(XX',ZZ',Fx(xxi2,zzi2)',Fz(xxi2,zzi2)');
            axis image
            axis([1800,3200,[-500,500]+zframe+500])
            title(['\fontsize{14} {Fx,Fz} and F_x (colourmap), \theta = ',...
                num2str(theta*180/pi,'%3.1f'),'^0,', ...
                ' \beta_S = ' num2str(betaS,'%3.2f')]);
            xlabel('x (in um)','FontSize',14);
            ylabel('z (in um)','FontSize',14);
            set(gca,'FontSize',14);
            grid on
            colorbar
            
            drawnow
        end;
        
        
    end;
    
    %store some spatial data for temporal evoluation plots
    % Ez(z,x=2.5mm,t)
  
    Ezt1(:,n) = Ez(420,:);
    Ezt2(:,n) = Ez(500,:);
    
    Hyt1(:,n) = Hy(420,:);
    Hyt2(:,n) = Hy(500,:);
    
    Ext1(:,n) = Ex(420,:);
    Ext2(:,n) = Ex(500,:);
  
end

