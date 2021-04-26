% 14/5/07

% Pedestrian lateral excitation
% No change in velocity at footstep - use velocity at end of previous step
% Adjust footstep width for stability

% Based on pedlatdif3 but changed to use Hof et al (2007) notation
% - 'absolute' position of CoM and CoP (relative to bridge)
% Average over long time period to cover different relative phases

warning off
clear all;
randn('state',0)
%figure;
% Pedestrian parameters (averages from Hof et al 2007)
%l = 0.9;		% Leg length (m)
r = 1.17;		% Effective pendulum length (m) = 1.34 * leg length [Hof et al 2005,2007]
%m = 77;			% Mass (kg)
g = 9.81;
m = 76.9;
phase=[]
p1=[]
l1={}
% Hof 'normal' walking (1m/s normalised)
%d = 0.086/2;	% Initial half step width (m)
%fp = 1/1.31/sqrt(l);			% Walking frequency (full lateral cycle) (Hz)
%bmin = 0.0161;	% Hof et al 2007 'margin' (m)
%bsd = 0.00378;

% Hof 'fast' walking (1.25m/s normalised)
%d = 0.088/2;	% Initial half step width (m)
%fp = 1/1.19/sqrt(l);			% Walking frequency (full lateral cycle) (Hz)
%bmin = 0.0181;	% Hof et al 2007 'margin' (m)
%bsd = 0.00403; 
%%gamma = 0.2;	% Hof fig 3 gradient - 1
fb_all = sort([0.125:.0125:3]);
Nfp=64
Fh=zeros(Nfp*max(size(fb_all)),1)
Fj=zeros(Nfp,max(size(fb_all)))
Nph=1
for toffs=1:1
for ip=1:Nfp
%for ip=1:Nfp
% Modified parameters to match inverted pendulum results for no bridge motion
fp = 0.65+0.6*(ip-1)/Nfp;			% Walking frequency (full lateral cycle) (Hz)
%fp = 0.7;			% Walking frequency (full lateral cycle) (Hz)
%fp=1.0;
d = 0.046;	% Initial half step width (m)
bmin = 0.0157;	%  'margin' (m)
bsd = 0.00403; 

% Modified parameters to match inverted pendulum results for no bridge motion, 0.8Hz
%fp = 0.8;			% Walking frequency (full lateral cycle) (Hz)
%d = 0.049;	% Initial half step width (m)
%bmin = 0.0141;	%  'margin' (m)


P = 0.5;	% Stance time / stride time.  0.5 = instantaneous transfer from one foot to the other

%K = 1700;	% Stabilising spring stiffness (Donelan et al. 2004) 

% Pedestrian derived parameters
%r = 1.17;		% Effective pendulum length (m) = 1.34 * leg length [Hof et al 2005,2007]
sgr = sqrt(g/r);
r=1.17;
%yfac = 2*(2*fp)^2/sgr^2
yfac = 0;
ydotfac = 0;
%ymax = d*(1-1/cosh(-sgr/4/fp));
%ch = cosh(sgr/2/fp);
%sh = sinh(sgr/2/fp);
%yfac = ch/(ch-1)-1
%ydotfac = sh/(ch-1)
tref = (1.5-P);

% Bridge parameters
% Frequency of vibration (Hz)
%fb = 1/2.01*fp;	%appox 0.4Hz
%%fb = 1/1.69*fp;	%appox 0.53Hz; % CSB Mode L2	
%fb = 1/1.09*fp;	%appox 0.74Hz; % CSB Mode L3
%%fb = 1/1.01*fp;	%appox 0.89Hz; % Brownjohn Mode LS1
%%fb = 1/0.87*fp;	%appox 1.03Hz; % Millennium Bridge Mode NL1
%fb = 1/0.51*fp;	%appox 1.6Hz;
%fb = 1/0.41*fp;	%appox 2Hz;
%fb = 1/0.27*fp;	%appox 3Hz;

%fb_all = fp./[2.01:-0.02:0.01];

toD=[8, 32, 56, 76]./4;

%fb_all = fp./[0.83:-0.02:0.71]; % range for 0.7Hz walking frequency
cmap=[0.1,0.65,0.7, 0.2,0.15,0.2, 1,0.2,0.25, 0.8,0.2,0.25, 0.1,0.0,1.0];
cmap=reshape(cmap,[3 5])';
fbs = fb_all;
F=[];
it=1;
for fi=1:max(size(fbs))
fb=fbs(fi)
tau=0;
	Nd=10;
	clear k2;
	clear k4;
	k2=zeros(Nd,1);
	Nk2=zeros(Nd,1);
	k4=zeros(Nd,1);
	X = 0.006;		% Amplitude (m)
	N = 2000;		% No of cylces
	n = 5000;	% No of time steps per footstep
	n0 = n;

	dt = 1/2/fp/n;
	tmax = 1/fb*N;
	t = (dt:dt:tmax*2)';



	%tmod = mod(t,1/fp);

	% Bridge derived parameters
	ob = fb*2*pi;	% bridge circular natural frequency 
	a = X/(1+(sgr/ob)^2);
	%tau = linspace(0,1/fp,181);
	x = X*sin(ob*(t-tau));
	xdot = X*ob*cos(ob*(t-tau));

	% Set parameters for initial step
	clear y ydot ydotdot
	t0(1) = 0;
	y0(1) = 0;
	u(1) = d;
	ydot0(1) = d*sgr*tanh(1/4/fp*sgr);
	ydotdot0(1) = sgr*sgr*(u(1));
	ydot=ydot0;
	foot(1) = -1; % +ve for right foot, -ve for left foot
	imin = 1;
	t0 =  [0];
	k=[];

	k3=[];

	for j=1:N-1
	   %imax = max(find(t<=t0(j)+1/2/fp));
	   n;
	   A(j) = y0(j) - u(j) + a*sin(ob*(t0(j)-tau));
	   B(j) = ydot0(j)/sgr + (ob/sgr)*a*cos(ob*(t0(j)-tau));
	   imax = imin-1+n;
	   S=size(t);
	   imax=min(imax,S(1));
	   if imax<=imin
		break
	   end
		y(imin:imax,:) = u(j) + A(j)*cosh(sgr*(t(imin:imax)-t0(j))) + B(j)*sinh(sgr*(t(imin:imax)-t0(j))) - a*sin(ob*(t(imin:imax)-tau));
		   ydot(imin:imax,:) = A(j)*sgr*sinh(sgr*(t(imin:imax)-t0(j))) + B(j)*sgr*cosh(sgr*(t(imin:imax)-t0(j))) - a*ob*cos(ob*(t(imin:imax)-tau));
		   ydotdot(imin:imax,:) = A(j)*sgr^2*cosh(sgr*(t(imin:imax)-t0(j))) + B(j)*sgr^2*sinh(sgr*(t(imin:imax)-t0(j))) + a*ob^2*sin(ob*(t(imin:imax)-tau));
		   	   u(j+1)=y(imax,1)+ydot(imax,1)/sgr+foot(j)*bmin;
	   foot(j+1)=-foot(j);

	   
	   i0 = round(imin-1+(imax-imin+1)*tref);
	   %T=t0:dt:(t0+n*dt);
	   y0(j+1) = y(imax);
	   ydot0(j+1) = ydot(imax);
	   ydotdot0(j+1) = ydotdot(imax);
	   
	   
	   % absolute
	   %z = y(i0)+x(i0);
	   %zdot = ydot(i0)+xdot(i0);
	   %x0 = x(i0);
	   
	   % relative
	   z = y(i0);
	   zdot = ydot(i0);
	   x0 = 0;
	   
	   % half way between absolute and relative (for yfac = 0)
	   %z = y(i0);
	   %zdot = ydot(i0)+xdot(i0)/2;
	   %x0 = 0;

		   %u(j+1) = y0(j+1) + ydot0(j+1)/sgr + foot(j+1)*bmin;  % Hof et al 2007 (increasing bmin has little effect)
		   %u(j+1) = y0(j+1)*(1+yfac) + ydot0(j+1)/sgr + foot(j+1)*bmin;  % Hof et al 2007,  modified - gives positive k
		%u(j+1) = y0(j+1) + ydot0(j+1)/sgr/tanh(1/4/fp*sgr);  % proportional to velocity, target step duration
	   %u(j+1) = (y0(j+1) + ydot0(j+1)/sgr)*1.15 + foot(j+1)*bmin;  % Hof et al 2007, including gradient of CoP vs XcoM
	   %u(j+1) = y0(j+1)*(1+yfac) + ydot0(j+1)/sgr*ydotfac; % Target zero displacement at end of step
		%u(j+1) = (y0(j+1) + foot(j+1)*bmin)*2;  % based on displacement only - unstable, at least for these parameters
	   
	      %v0(j+1) = (y(imax)-y(imax-1)+x(imax)-x(imax-1))/dt; % absolute velocity    
	      %d(j+1) = abs(v0(j+1))/sgr/tanh(1/4/fp*sgr); % proportional to velocity, target step duration
	      %d(j+1) = d(1)+(abs(v0(j+1))-v0(1))/sgr/tanh(1/4/fp*sgr)*c1; % linear function of velocity
	      %d(j+1) = d(1)+foot*y(imax)*c2; % linear function of displacement (c2 = 1 = width between feet constant)
	      %d(j+1) = d(1) + (abs(v0(j+1))-v0(1))/sgr/tanh(1/4/fp*sgr)*c1 + foot*y(imax)*c2; % combined feedback
	      %d(j+1) = -foot*v0(j+1)/sgr + bmin;  % Hof et al 2007      
	      %d(j+1) = -foot*v0(j+1)/sgr*0.8 + 0.0238;
	      %d(j+1) = -foot*(v0(j+1)/sgr*(1+gamma) + (z1(imax))*gamma) + bmin; %+ x(imax)*gamma ???
		
   	    fs=(r*fp./1.34./1.35).^2*fp;
		    n=ceil(n0*(1+((0.046)^2-(y(imax,1)-u(j+1))^2)./(4*fs*fs)));
		
		    if n<=1
		       n=2;
		    end
	    	    t0(j+1)=t0(j)+n*dt;
		    
		    H=-m*(-sin(ob*(t(imin:imax)))*X*ob^2+ydotdot(imin:imax));
		    k(j) = trapz(2.0*H.*(cos(ob*(t(imin:imax)-tau))).*dt)./(t0(j+1)-t0(j));
		    k3(j)=mean(abs(u(j)-y(imin:imax)));
		    imin = imax+1;
		end

		y0 = y0(1:end-1);
		ydot0 = ydot0(1:end-1);
		foot = foot(1:end-1);
		u = u(1:end-1);

		XcoM = y+ydot/sgr;
		XcoM0 = y0+ydot0/sgr;

		u2 = reshape([u;u],1,length(u)*2);
		t2 = reshape([t0;t0],1,length(t0)*2);
		t2 = [t2(2:end),t0(end)*2-t0(end-1)];

		t = t(1:length(y));
		x = x(1:length(y));

		% Find component of force in phase with bridge velocity
		start = 1;
		stop = length(t0);


		% Component in phase with bridge displacement
		phase=(mod(t0-tau,1/fb))*fb;
		k=k';
		for j=start:max(1,stop-1)
		    k2(max(1,ceil(phase(j)*(Nd))))=k2(max(1,ceil(phase(j)*(Nd))))+k(j);
		    Nk2(max(1,ceil(phase(j)*(Nd))))=Nk2(max(1,ceil(phase(j)*(Nd))))+1;
		    k4(max(1,ceil(phase(j)*(Nd))))=k3(j);
		end

		if k2(1)>9998
		    k2(1)=k2(end);
		else
		    k2(end)=k2(1);
		end
		P=0:(1/Nd):1.0;
	k2./Nk2;
	F(fi)=trapz(P(Nk2>0),k2(Nk2>0)./Nk2(Nk2>0),1);
	Fj(ip,fi)=1./mean(diff(t0(start:2:stop)));
	clear y; 
	clear ydot;
	clear y0;
	clear ydot0;
	%clear t0;
	clear phase;
	clear ydotdot;
	clear u; 
	clear x; 
	clear t;
	clear u2;
	clear t2;
	clear A;
	clear B;
        F(fi);
	Fh((ip-1)*max(size(fb_all))+fi,toffs)=F(fi);
	if(sum(toD==fi))
	it=it+1;
	%area(p(K>0), K(K>0), 'FaceColor', [1,0,0], 'FaceAlpha', 0.3, 'EdgeAlpha',0);
	%ax2=subplot(3,1,2);
	%plot(P(k4<9998.0),k4(k4<9998.0), 'Color', lcol, 'LineWidth', 3.5); hold on;
	%plot(xlim(), [0.04 0.04], 'r--', 'LineWidth',3.5); hold on;
	%xlabel('\tau f_b','FontSize',55)
	%ylabel('H','FontSize',55)
	end
end
end
end

[XX,YY]=meshgrid(fbs, linspace(0.6,60,1.2))
csvwrite('XX.csv', XX)
csvwrite('YY.csv', YY)
csvwrite('Fj.csv', Fj)
csvwrite('Fh.csv', Fh)
save('fratio.mat')
    

