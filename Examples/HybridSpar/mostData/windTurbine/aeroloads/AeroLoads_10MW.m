%% INIT
clearvars -except WindTurbine_type
clc
%% SETTINGS 
if 1
%% Plot
plotflag=1;
loadplot_delta_wind=8;
%% Look up table settings
wind_vec=-2:0.5:25;
o_discr=17;                                % rotor speed discretisation points
theta_discr=15;                            % blade pitch discretisation points
o_A =3.6*pi/30;                            % amplitude of rotor speed [rad/s]
theta_A = 12*pi/180;                       % amplitude of blade pitch [rad]  
window_outlier=7;
%% BEM Settings
BEM_data.rho_air=1.225;
BEM_data.root_tol_rel=1/200;
BEM_data.maxit=20;
BEM_data.func_tol2=1e-5;
BEM_data.root_tol_rel2=1/200;   
BEM_data.maxit2=200;
BEM_data.eps=1e-6;
end
%% DATA 
if 1
%% Load files
cd ..
addpath('turbine_properties') 
addpath('control')
cd aeroloads
load("Properties_10MW.mat");   % Load windturbine properties
load("Bladedata_10MW.mat");    % Load blade data
load('SteadyStates_10MW.mat'); % Load steady states data 
%% Look up table 
o_discr=o_discr+double(mod(o_discr,2)==0);
theta_discr=theta_discr+double(mod(theta_discr,2)==0);
omegaerr=linspace(-o_A,o_A,o_discr);
bladepitcherr=linspace(-theta_A,theta_A,theta_discr);

wind_vec=unique([0 wind_vec]);
%% BEM data
BEM_data.RhubRblade=[bladedata.radius(1) bladedata.radius(end)];
BEM_data.twist=0.5*(bladedata.twist(1:end-1)+bladedata.twist(2:end));
BEM_data.chord=0.5*(bladedata.chord(1:end-1)+bladedata.chord(2:end));
BEM_data.BlCrvAng=pi/180*0.5*(bladedata.BlCrvAng(1:end-1)+bladedata.BlCrvAng(2:end));
BEM_data.BlSwpAC=0.5*(bladedata.BlSwpAC(1:end-1)+bladedata.BlSwpAC(2:end));
BEM_data.BlCrvAC=0.5*(bladedata.BlCrvAC(1:end-1)+bladedata.BlCrvAC(2:end));
BEM_data.r=0.5*(bladedata.radius(1:end-1)+bladedata.radius(2:end));
BEM_data.r_int=bladedata.radius(2:end)'-bladedata.radius(1:end-1)';
BEM_data.airfoil=bladedata.airfoil;
BEM_data.airfoil_index=bladedata.airfoil_index(1:end-1);
BEM_data.precone=WTcomponents.hub.precone;

end
%% CALC ROSCO AEROLOADS LOOK UP TABLE 
if 1
%% Compute F_aero
SS=SteadyStates.ROSCO.SS;
omega_vec=max(0,interp1([0 SS.WINDSPEED],[0 SS.ROTSPD],wind_vec,"linear","extrap"));
bladepitch_vec=max(0,interp1([0 SS.WINDSPEED],[0 SS.BLADEPITCH],wind_vec,"linear","extrap"));

F_aero=zeros(length(wind_vec),o_discr,theta_discr,6);

for i=1:length(wind_vec)
    for ii=1:o_discr
        for iii=1:theta_discr
            [F_aero(i,ii,iii,:)]=BEM(wind_vec(i),omega_vec(i)+omegaerr(ii),bladepitch_vec(i)+bladepitcherr(iii),BEM_data);
        end
    end
end

%% Remove outliers
F_aero_no_out=0*F_aero;
for i=1:6
    for j=1:length(wind_vec)

        Z=squeeze(F_aero(j,:,:,i));
        outliers_i = isoutlier(Z,'movmedian',window_outlier,1);
        outliers_j = isoutlier(Z,'movmedian',window_outlier,2);
        [outliers_i_idx,outliers_j_idx]=find(outliers_i & outliers_j);
        [X,Y]=meshgrid(bladepitcherr,omegaerr);
        X2=X;Y2=Y;Z2=Z;
        X2(unique(outliers_i_idx),:)=[];
        X2(:,unique(outliers_j_idx))=[];
        Y2(unique(outliers_i_idx),:)=[];
        Y2(:,unique(outliers_j_idx))=[];
        Z2(unique(outliers_i_idx),:)=[];
        Z2(:,unique(outliers_j_idx))=[];
        try
            Z_interp= interp2(X2, Y2, Z2, X, Y, 'spline');
            Z(outliers_i & outliers_j)=Z_interp(outliers_i & outliers_j);
            F_aero_no_out(j,:,:,i)=Z;
        catch
            F_aero_no_out(j,:,:,i)=F_aero(j,:,:,i);
        end
    end
end
% Remove outliers 2
% F_aero_no_out2=0*F_aero;
% for i=1:6
%     for j=1:length(wind_vec)
% 
%         Z=squeeze(F_aero(j,:,:,i));
%         outliers_i = isoutlier(Z,'movmedian',window_outlier,1);
%         outliers_j = isoutlier(Z,'movmedian',window_outlier,2);
%         Z_nan=Z;
%         Z_nan(outliers_i & outliers_j)=nan;
%         F_aero_no_out2(j,:,:,i)=fillmissing2(Z_nan,"linear");
% 
%     end
% end

%% Create aeroloads struct
aeroloads.ROSCO.SS=SS;
aeroloads.ROSCO.wind_vec=wind_vec;
aeroloads.ROSCO.omega_vec=omega_vec;
aeroloads.ROSCO.bladepitch_vec=bladepitch_vec;
aeroloads.ROSCO.omegaerr=omegaerr;
aeroloads.ROSCO.bladepitcherr=bladepitcherr;
aeroloads.ROSCO.F_aero=F_aero_no_out;

%% Plots
if plotflag
%% Breakpoints lookup-table
figure()
subplot(2,1,1)
scatter(wind_vec, omega_vec, 20, 'r', 'filled');
hold on
errorbar(wind_vec, omega_vec, omegaerr(end)-omegaerr(1), 'b', 'LineWidth', 1); 
xlabel('Wind Speed (m/s)');
ylabel('Rotor Speed (rad/s)');
grid on;
legend('Steady States Value','Deltas');
subplot(2,1,2)
scatter(wind_vec, bladepitch_vec*180/pi, 20, 'r', 'filled'); 
hold on
errorbar(wind_vec, bladepitch_vec*180/pi, (bladepitcherr(end)-bladepitcherr(1))*180/pi, 'b', 'LineWidth', 1);
xlabel('Wind Speed (m/s)');
ylabel('Bladepitch (deg)');
grid on;
legend('Steady States Value','Deltas');
sgtitle('ROSCO: Look-up table Breakpoints');

%% Loads
loadNames={'F_x','F_y','F_z','T_x','T_y','T_z'};
loadUdm={'N','N','N','Nm','Nm','Nm'};
for i=1:6
    figure()
    for j=1:loadplot_delta_wind:length(wind_vec)
        colormap(winter);
        s=surf(omega_vec(j)+omegaerr,(bladepitch_vec(j)+bladepitcherr)*180/pi,squeeze(F_aero_no_out(j,:,:,i))');
        s.FaceAlpha = 0.3;
        shading interp;
        
        hold on
        plot3(omega_vec(j),bladepitch_vec(j)*180/pi,F_aero(j,ceil(length(omegaerr)/2),ceil(length(bladepitcherr)/2),i),'.r','markersize',25)
        text(omega_vec(j),bladepitch_vec(j)*180/pi,F_aero(j,ceil(length(omegaerr)/2),ceil(length(bladepitcherr)/2),i),...
            [num2str(wind_vec(j),2) 'm/s'],"FontSize",12)

        title(['ROSCO Blade root loads: ' loadNames{i}])
        grid on
        xlabel('Rotor Speed (rad/s)')
        ylabel('Blade pitch (deg)')
        zlabel([loadNames{i} '(' loadUdm{i} ')'])
    end
end
sf=length(findobj(allchild(0), 'flat', 'Type', 'figure'));
spreadfigures(flip(sf:-1:sf-5))

end
end
%% CALC BASELINE AEROLOADS LOOK UP TABLE 
if 1
%% Compute F_aero
SS=SteadyStates.Baseline.SS;
omega_vec=max(0,interp1([0 SS.WINDSPEED],[0 SS.ROTSPD],wind_vec,"linear","extrap"));
bladepitch_vec=max(0,interp1([0 SS.WINDSPEED],[0 SS.BLADEPITCH],wind_vec,"linear","extrap"));

F_aero=zeros(length(wind_vec),o_discr,theta_discr,6);

for i=1:length(wind_vec)
    for ii=1:o_discr
        for iii=1:theta_discr
            [F_aero(i,ii,iii,:)]=BEM(wind_vec(i),omega_vec(i)+omegaerr(ii),bladepitch_vec(i)+bladepitcherr(iii),BEM_data);
        end
    end
end

%% Remove outliers
F_aero_no_out=0*F_aero;
for i=1:6
    for j=1:length(wind_vec)

        Z=squeeze(F_aero(j,:,:,i));
        outliers_i = isoutlier(Z,'movmedian',window_outlier,1);
        outliers_j = isoutlier(Z,'movmedian',window_outlier,2);
        [outliers_i_idx,outliers_j_idx]=find(outliers_i & outliers_j);
        [X,Y]=meshgrid(bladepitcherr,omegaerr);
        X2=X;Y2=Y;Z2=Z;
        X2(unique(outliers_i_idx),:)=[];
        X2(:,unique(outliers_j_idx))=[];
        Y2(unique(outliers_i_idx),:)=[];
        Y2(:,unique(outliers_j_idx))=[];
        Z2(unique(outliers_i_idx),:)=[];
        Z2(:,unique(outliers_j_idx))=[];
        try
            Z_interp= interp2(X2, Y2, Z2, X, Y, 'spline');
            Z(outliers_i & outliers_j)=Z_interp(outliers_i & outliers_j);
            F_aero_no_out(j,:,:,i)=Z;
        catch
            F_aero_no_out(j,:,:,i)=F_aero(j,:,:,i);
        end
    end
end
% Remove outliers 2
% F_aero_no_out2=0*F_aero;
% for i=1:6
%     for j=1:length(wind_vec)
% 
%         Z=squeeze(F_aero(j,:,:,i));
%         outliers_i = isoutlier(Z,'movmedian',window_outlier,1);
%         outliers_j = isoutlier(Z,'movmedian',window_outlier,2);
%         Z_nan=Z;
%         Z_nan(outliers_i & outliers_j)=nan;
%         F_aero_no_out2(j,:,:,i)=fillmissing2(Z_nan,"linear");
% 
%     end
% end



%% Create aeroloads struct
aeroloads.Baseline.SS=SS;
aeroloads.Baseline.wind_vec=wind_vec;
aeroloads.Baseline.omega_vec=omega_vec;
aeroloads.Baseline.bladepitch_vec=bladepitch_vec;
aeroloads.Baseline.omegaerr=omegaerr;
aeroloads.Baseline.bladepitcherr=bladepitcherr;
aeroloads.Baseline.F_aero=F_aero_no_out;

%% Plots
if plotflag
%% Breakpoints lookup-table
figure()
subplot(2,1,1)
scatter(wind_vec, omega_vec, 20, 'r', 'filled');
hold on
errorbar(wind_vec, omega_vec, omegaerr(end)-omegaerr(1), 'b', 'LineWidth', 1); 
xlabel('Wind Speed (m/s)');
ylabel('Rotor Speed (rad/s)');
grid on;
legend('Steady States Value','Deltas');
subplot(2,1,2)
scatter(wind_vec, bladepitch_vec*180/pi, 20, 'r', 'filled'); 
hold on
errorbar(wind_vec, bladepitch_vec*180/pi, (bladepitcherr(end)-bladepitcherr(1))*180/pi, 'b', 'LineWidth', 1);
xlabel('Wind Speed (m/s)');
ylabel('Bladepitch (deg)');
grid on;
legend('Steady States Value','Deltas');
sgtitle('Baseline: Look-up table Breakpoints');
%% Loads
loadNames={'F_x','F_y','F_z','T_x','T_y','T_z'};
loadUdm={'N','N','N','Nm','Nm','Nm'};
for i=1:6
    figure()
    for j=1:loadplot_delta_wind:length(wind_vec)
        colormap(winter);
        s=surf(omega_vec(j)+omegaerr,(bladepitch_vec(j)+bladepitcherr)*180/pi,squeeze(F_aero_no_out(j,:,:,i))');
        s.FaceAlpha = 0.2;
        shading interp;
        
        hold on
        plot3(omega_vec(j),bladepitch_vec(j)*180/pi,F_aero(j,ceil(length(omegaerr)/2),ceil(length(bladepitcherr)/2),i),'.r','markersize',25)
        text(omega_vec(j),bladepitch_vec(j)*180/pi,F_aero(j,ceil(length(omegaerr)/2),ceil(length(bladepitcherr)/2),i),...
            [num2str(wind_vec(j),2) 'm/s'],"FontSize",12)

        title(['Baseline Blade root loads: ' loadNames{i}])
        grid on
        xlabel('Rotor Speed (rad/s)')
        ylabel('Blade pitch (deg)')
        zlabel([loadNames{i} '(' loadUdm{i} ')'])
    end
end
sf=length(findobj(allchild(0), 'flat', 'Type', 'figure'));
spreadfigures(flip(sf:-1:sf-5))


end
end
%% SAVE 
save('Aeroloads_10MW','aeroloads')
%% FUNCTIONS
function F_aero = BEM(v_wind_rel,omega,bladepitch,data)

pntm=[0 0 0].*data.r;
Cn=0;   Ct=0;    CM=0;
a=0;
at=0;

for node=1:length(data.r)

    VXoverVY=v_wind_rel/(omega*data.r(node)*cosd(data.precone));
    s=3*data.chord(node)/2/pi/data.r(node);
    i=1;
    eps=1;

    while eps>data.root_tol_rel  &&  i<=data.maxit
        phi_i=atan((1-a)/(1+at)*VXoverVY);
        alfa=rad2deg(phi_i-bladepitch)-data.twist(node);

        CL=interp1(data.airfoil(:,1,data.airfoil_index(node)),data.airfoil(:,2,data.airfoil_index(node)),alfa);
        CD=interp1(data.airfoil(:,1,data.airfoil_index(node)),data.airfoil(:,3,data.airfoil_index(node)),alfa);
        CM=interp1(data.airfoil(:,1,data.airfoil_index(node)),data.airfoil(:,4,data.airfoil_index(node)),alfa);
        Cn=CL*cos(phi_i)+CD*sin(phi_i);
        Ct=CL*sin(phi_i)-CD*cos(phi_i);

        Ftip=2/pi*acos(exp(-1.5*(data.RhubRblade(2)-data.r(node))/data.r(node)/abs(sin(phi_i+1e-5)))); % loss tip
        Fhub=2/pi*acos(exp(-1.5*(data.r(node)-data.RhubRblade(1))/data.r(node)/abs(sin(phi_i+1e-5)))); % loss hub
        F=Ftip*Fhub;

        k=s*Cn/4/F/sin(phi_i)^2;
        kt=s*Ct/4/F/sin(phi_i)/cos(phi_i);
        if (phi_i>0 && k<=2/3)
            a=k/(1+k);

        elseif (phi_i>0 && k>2/3)
            g1=2*F*k-(10/9-F);
            g2=2*F*k-F*(4/3-F);
            g3=2*F*k-(25/9-2*F);
            a=(g1-sqrt(g2))/g3;

        elseif (phi_i<0 && k>1)
            a=k/(k-1);

        elseif (k<=1 && phi_i<0)
            a=0;
        end

        at=kt/(1-kt);

        phi_f=atan((1-a)/(1+at)*VXoverVY);
        i=i+1;
        eps=abs((phi_i-phi_f)/phi_i);


    end


    if i>data.maxit || a==0

        if R(data.eps,v_wind_rel,bladepitch,omega,s,node,data)*R(pi/2-data.eps,v_wind_rel,bladepitch,omega,s,node,data)<0
            phi=Root_find(@(phi)R(phi,v_wind_rel,bladepitch,omega,s,node,data),[data.eps pi/2],data);

        elseif Rpb(-pi/4,v_wind_rel,bladepitch,omega,s,node,data)*Rpb(-data.eps,v_wind_rel,bladepitch,omega,s,node,data)<0
            phi=Root_find(@(phi)Rpb(phi,v_wind_rel,bladepitch,omega,s,node,data),[-pi/4 -data.eps],data);

        else
            phi=Root_find(@(phi)R(phi,v_wind_rel,bladepitch,omega,s,node,data),[pi/2+data.eps pi-data.eps],data);

        end

        alfa=rad2deg(phi-bladepitch)-data.twist(node);
        CL=interp1(data.airfoil(:,1,data.airfoil_index(node)),data.airfoil(:,2,data.airfoil_index(node)),alfa);
        CD=interp1(data.airfoil(:,1,data.airfoil_index(node)),data.airfoil(:,3,data.airfoil_index(node)),alfa);
        CM=interp1(data.airfoil(:,1,data.airfoil_index(node)),data.airfoil(:,4,data.airfoil_index(node)),alfa);
        Cn=CL*cos(phi)+CD*sin(phi);
        Ct=CL*sin(phi)-CD*cos(phi);
        Ftip=2/pi*acos(exp(-1.5*(data.RhubRblade(2)-data.r(node))/data.r(node)/abs(sin(phi+1e-5))));
        Fhub=2/pi*acos(exp(-1.5*(data.r(node)-data.RhubRblade(1))/data.r(node)/abs(sin(phi+1e-5))));
        F=Ftip*Fhub;

        kt=s*Ct/4/F/sin(phi)/cos(phi);
        at=kt/(1-kt);
        Vrel2=((omega*data.r(node))*(1+at)/cos(phi))^2;

    else

        Vrel2=(v_wind_rel*(1-a))^2+((omega*data.r(node))*(1+at))^2;

    end
    pntm(node,:)=0.5*data.rho_air*data.chord(node)*Vrel2*[Cn,Ct,CM*data.chord(node)];


    if anynan([a at])
        a=0; at=0;
    end

end

% Remove NaN

if anynan(pntm)

    temp_y=[[0 0 0];pntm;[0 0 0]];
    temp_x=(1:size(temp_y,1))';
    temp_y=[interp1(temp_x(~isnan(temp_y(:,1))),temp_y(~isnan(temp_y(:,1)),1),temp_x,'linear'),...
            interp1(temp_x(~isnan(temp_y(:,2))),temp_y(~isnan(temp_y(:,2)),2),temp_x,'linear'),...
            interp1(temp_x(~isnan(temp_y(:,3))),temp_y(~isnan(temp_y(:,3)),3),temp_x,'linear')];
    pntm=temp_y(2:end-1,:);

end


F_aero=[ data.r_int*(pntm(:,1).*cos(data.BlCrvAng));
        -data.r_int*pntm(:,2);
         data.r_int*(pntm(:,1).*sin(-data.BlCrvAng));
         data.r_int*(pntm(:,2).*(data.r-data.RhubRblade(1))*cosd(data.precone)+pntm(:,3).*sin(data.BlCrvAng)+pntm(:,1).*sin(-data.BlCrvAng).*data.BlSwpAC);
         data.r_int*(pntm(:,1).*(data.r-data.RhubRblade(1)));
         data.r_int*(pntm(:,3).*cos(data.BlCrvAng)-pntm(:,1).*cos(data.BlCrvAng).*data.BlSwpAC-pntm(:,2).*data.BlCrvAC)];




end

function [root] = Root_find(Fun,Int,data)

a = Int(1,1);
b = Int(1,2);
fa = Fun(a);
fb = Fun(b);
it = 0;
s=b;

while abs((b-a)/a)>data.root_tol_rel2  &&  it<data.maxit2  &&  abs(fb)>data.func_tol2

    it = it + 1;
    s = (a + b)/2;
    fs = Fun(s);

    if fa*fs < 0
        b = s;
        fb = fs;
    else
        a = s;
        fa = fs;
    end

end

root = s;

end

function resR=R(phi,U_inf,blpitch,omega,solidity,node,data)


alfa=rad2deg(phi-blpitch)-data.twist(node);

CL=interp1(data.airfoil(:,1,data.airfoil_index(node)),data.airfoil(:,2,data.airfoil_index(node)),alfa);
CD=interp1(data.airfoil(:,1,data.airfoil_index(node)),data.airfoil(:,3,data.airfoil_index(node)),alfa);

Cn=CL*cos(phi)+CD*sin(phi);
Ct=CL*sin(phi)-CD*cos(phi);

Ftip=2/pi*acos(exp(-1.5*(data.RhubRblade(2)-data.r(node))/data.r(node)/abs(sin(phi+1e-5)))); 
Fhub=2/pi*acos(exp(-1.5*(data.r(node)-data.RhubRblade(1))/data.r(node)/abs(sin(phi+1e-5)))); 
F=Ftip*Fhub;

k=solidity*Cn/4/F/sin(phi)^2;
kt=solidity*Ct/4/F/sin(phi)/cos(phi);

if  k<=2/3                
    a=k/(1+k);                         

else
    g1=2*F*k-(10/9-F);
    g2=2*F*k-F*(4/3-F);
    g3=2*F*k-(25/9-2*F)+1e-6;
    a=(g1-sqrt(g2))/g3;                

end

at=kt/(1-kt);                          

resR=sin(phi)/(1-a)-cos(phi)/(omega*data.r(node)/U_inf*(1+at));


end

function resRpb=Rpb(phi,U_inf,blpitch,omega,solidity,node,data)

alfa=rad2deg(phi-blpitch)-data.twist(node);

CL=interp1(data.airfoil(:,1,data.airfoil_index(node)),data.airfoil(:,2,data.airfoil_index(node)),alfa);
CD=interp1(data.airfoil(:,1,data.airfoil_index(node)),data.airfoil(:,3,data.airfoil_index(node)),alfa);

Cn=CL*cos(phi)+CD*sin(phi);
Ct=CL*sin(phi)-CD*cos(phi);

Ftip=2/pi*acos(exp(-1.5*(data.RhubRblade(2)-data.r(node))/data.r(node)/abs(sin(phi+1e-5)))); 
Fhub=2/pi*acos(exp(-1.5*(data.r(node)-data.RhubRblade(1))/data.r(node)/abs(sin(phi+1e-5)))); 
F=Ftip*Fhub;

k=solidity*Cn/4/F/sin(phi)^2;
kt=solidity*Ct/4/F/sin(phi)/cos(phi);
                       
resRpb=sin(phi)*(1-k)-cos(phi)/(omega*data.r(node)/U_inf)*(1-kt);


end


