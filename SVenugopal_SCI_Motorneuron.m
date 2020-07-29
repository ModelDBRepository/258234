%% % A two-compartment motoneuron model (MN)
% Author: Sharmila Venugopal (vsharmila@ucla.edu; svenugopal10@gmail.com)
% The motor neuron model reproduces dendritic persistent Ca2+-induced sustained MN discharge
% following chronic spinal injury
% If you use/benefit from the model, please cite: Venugopal S, et al, J Neurophysiology, 2010
%% 
clc;
clear all;
close all;
%% Functions calls for setting initial conditions, parameter values and generating a ramp current
%setParams();
y0 = setInit();
%% Simulation time and time steps
tmax=10000;
tstep=0.1;
tspan=0:tstep:tmax;
%% Solve using a stiff method: ode23s (DO NOT USE ode45)
[t,y]=ode23s(@rates,tspan,y0);
Vs = y(:,1);
Vd = y(:,7);
makeplots(t,Vs,Vd);
%% BEGIN FUNCTIONS
%% setInit(): sets and returns an initial value vector
function y0 = setInit()
    % soma initial conditions
    vs=-57.34;
    snah=0.5829;
    sn=0.1239;
    scam=0.004199;
    scah=0.9219;
    sca=0.0001406;
    % dendrite initial conditions
    vd=-56.64;
    dcas=0.08493;
    dca=0.01724;
    mnap=0.1; 
    hnap=0.9;
    
    y0 = [vs; snah; sn; scam; scah; sca; vd; dcas; dca; mnap; hnap];
end
%%     %% Voltage-dependent functions for ion channel gating
    %% % Somatic channel gating functions
    function ss_snaminf = snaminf(vs)
        snamth=-35.0; snamslp=7.8;
        ss_snaminf=1/(1+exp(-(vs-snamth)/snamslp));
    end

    function ss_snahinf = snahinf(vs)
        snahth=-55.0; snahslp=7.0;
        ss_snahinf=1.0/(1.0+exp((vs-snahth)/snahslp));
    end

    function ss_snahtau = snahtau(vs)
        snahv=-50.0; snaha=15.0; snahb=16.0; snahc=30.0;
        ss_snahtau=snahc/(exp((vs-snahv)/snaha)+exp(-(vs-snahv)/snahb));
    end

    function ss_sninf = sninf(vs)
        nth=-28.0; nslp=12.0;
        ss_sninf=1.0/(1.0+exp(-(vs-nth)/nslp));
    end

    function ss_sntau = sntau(vs)
        nv=-40.0; na=40.0; nb=50.0; nc=7.0;
        ss_sntau=nc/(exp((vs-nv)/na)+exp(-(vs-nv)/nb));
    end

    function ss_scaminf = scaminf(vs)
        camth=-30.0; camslp=5.0;
        ss_scaminf=1.0/(1.0+exp(-(vs-camth)/camslp));
    end

    function ss_scahinf = scahinf(vs)
        cahth=-45.0; cahslp=5.0;
        ss_scahinf=1.0/(1.0+exp((vs-cahth)/cahslp));
    end
    %% Dendritic channel gating functions
    function ss_dcasinf = dcasinf(vd)
        dcasth=-39.0; dcasslp=7.0;
        ss_dcasinf=1.0/(1.0+exp(-(vd-dcasth)/dcasslp));
    end

    function ss_minfnap = minfnap(vd)
        thetamnap=-48; Kmnap=3;
        ss_minfnap=1/(1+exp(-(vd-thetamnap)/Kmnap));
    end

    function ss_hinfnap = hinfnap(vd)
        thetahnap=-35; Khnap=6;
        ss_hinfnap=1/(1+exp((vd-thetahnap)/Khnap));
    end
%% rates(t,y): incorporates the rate equations of the model
function dydt = rates(t,y)

    %% % Model parameters
    vna=115.0; vk=-20.0; vl=0.0; vca=140.0; vrest=-60.0; cm=1.0; gl=0.51;
    p=1; kd=0.2; f=0.01; kca=2.0; alpha=0.009;
    % Soma parameters
    sgna=80.0; sgk=100.0; sgca=14.0; camtau=4.0; cahtau=40.0; sgkca=6;
    % Dendrite parameters
    dcastau=40.0; taunap=40; naptau=1000; 
    % Soma-dendrite coupling parameters
    gc=0.1; parea=0.1; perkca=1.0; 
    % Parameters sensitive to SCI
    dgkca=1.0; gnap=0.1; dgcas=0.25; alpha2=0.009; pernap=1;    
    %% Assign state variables    
    vs = y(1); snah = y(2); sn = y(3); scam = y(4); scah = y(5); sca = y(6); vd = y(7); dcas = y(8); dca = y(9); mnap = y(10); hnap = y(11);
    
    %% Get instantaneous value for the injected ramp current 
    Iramp = ramp(t);
    %% 
    % Somatic variables
    dydt(1)=(Iramp-sgna*snaminf(vs)^3*snah*(vs-(vna+vrest))-sgca*scam^2*scah*(vs-(vca+vrest))-sgk*sn^4*(vs-(vk+vrest))-perkca*sgkca*(sca^p/(kd^p+sca^p))*(vs-(vk+vrest))-gl*(vs-(vl+vrest))+gc*(vd-vs)/parea)/cm;
    dydt(2)=(snahinf(vs)-snah)/snahtau(vs);
    dydt(3)=(sninf(vs)-sn)/sntau(vs);
    dydt(4)=(scaminf(vs)-scam)/camtau;
    dydt(5)=(scahinf(vs)-scah)/cahtau;
    dydt(6)=-f*alpha*sgca*scam^2*scah*(vs-(vca+vrest))-f*kca*sca;
    
    % Dendritic variables
    dydt(7)=(-dgcas*dcas*(vd-(vca+vrest))-pernap*gnap*mnap*hnap*(vd-(vna+vrest))-perkca*dgkca*(dca^p/(kd^p+dca^p))*(vd-(vk+vrest))-gl*(vd-(vl+vrest))+gc*(vs-vd)/(1.0-parea))/cm;
    dydt(8)=(dcasinf(vd)-dcas)/dcastau;
    dydt(9)=-f*alpha2*dgcas*dcas*(vd-(vca+vrest))-f*kca*dca;
    dydt(10)=(minfnap(vd)-mnap)/taunap;
    dydt(11)=(hinfnap(vd)-hnap)/naptau;

    dydt = dydt';
end
%% ramp(t): Current ramp generator
function Iramp = ramp(t)
    offsetr=0;
    scaler=0.005;
    switchr=4000;

    if t<=4000
        Iramp = scaler*t+offsetr;    
    else
        Iramp = -scaler*t + (2*scaler*switchr);       
    end
end
%% makeplots(t,Vs,Vd): Function to draw output plots
function makeplots(t,Vs,Vd)
    XMIN = 0;
    XMAX = 10000;
    YMIN = -70;
    YMAX = 50;
    
    figure
    plot(t,Vs,'k','Linewidth',0.75);
    hold on
    plot(t,Vd,'r','Linewidth',1);
    axis([XMIN XMAX YMIN YMAX]);
    xlabel('Time (ms)');
    ylabel('V (mV)');
    legend('V_{soma}', 'V_{dendrite}');
end