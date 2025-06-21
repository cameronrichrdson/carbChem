clear all; close all;

% Input 
%==========================================================================
prompt = "Input Temp (C):";
T = input(prompt);         % Temp in Celsius (Lab ~21C most of the time)

prompt = "Input Salinity (g/kg):";
S = input(prompt);   
      
L = 1;       % Liters of synthetic seawater to be prepared

% Molality (mol/Kg art.SW) of constituent ions in the ASW 
    %  The Molality (mol/kg) value for ions from 
    %  "Guide to best practices for ocean CO2 measurements"
m_Cl = (0.54922);
m_SO4 = (0.02824);
m_Na = (0.46911);
m_Mg = (0.05283);
m_Ca = (0.01036);
m_K = (0.01021);

m_B = (0.00042); %assuming all boron from one source
%m_HCO3 = (0.00177); %assuming all alk from HCO3
m_CO32 = 0.0011; 

% Concentration of MgCl2 & CaCl2 solutions (mol/L) to be used
% (arbitarily set to 1 mol/L for now, to simplify calculation)
MgCl2_conc = 0.92;
CaCl2_conc = 0.984; 

%==========================================================================



% Equation of State for Seawater Density
% Seawater Density as a function of T & S at the normal atmospheric
% pressure (i.e., P = 0). 
% Millero & Poisson ('81) International one-atmosphere equation of state of
% seawater.

% Density of pure water
rhow = 999.842594 + 6.793952e-2*T -9.095290e-3*T^2 ...
            + 1.001685e-4*T^3 -1.120083e-6*T^4 + 6.536332e-9*T^5;

% Density of seawater at 1 atm, P=0
A =   8.24493e-1 - 4.0899e-3*T + 7.6438e-5*T^2 - 8.2467e-7*T^3 ...
    + 5.3875e-9*T^4;
B = -5.72466e-3 + 1.0227e-4*T - 1.6546e-6*T^2; 
C = 4.8314e-4;   

density = rhow + A*S + B*S^(3/2) + C*S^2; % unit -> kg/m^3  
                                                   % 1 m^3 = 1000 L
Density = density/1000;                   % unit -> kg/L  


% Moles of constituent ions
%    (mol/kg)*(kg/L)*(L) = mol
Cl = (m_Cl)*(Density)*L;
SO4 = (m_SO4)*(Density)*L;
Na = (m_Na)*(Density)*L;
Mg = (m_Mg)*(Density)*L;
Ca = (m_Ca)*(Density)*L;
K = (m_K)*(Density)*L;

B = (m_B)*(Density)*L;
%HCO3 = (m_HCO3)*(Density)*L;
CO32 = (m_CO32)*(Density)*L;


% Molecular mass of salts
Na2SO4_mw = 142.04;
KCl_mw = 74.55;
NaCl_mw = 58.44;

%alk 
BH3O3_mw = 61.83;
%NaHCO3_mw = 84.01;
Na2CO3_mw = 105.99;




% BH3O3 to be added (g)
BH3O3 = B * BH3O3_mw;
    % moles of B coming from BH3O3
    b_oh = B * 3

    
% Na2CO3 to be added (g)
Na2CO3 = CO32 * Na2CO3_mw;
    % moles of Na+ coming from Na2CO3
       sod_co32 = CO32 * 2;

% NaHCO3 to be added (g)
%NaHCO3 = HCO3 * NaHCO3_mw;
    % moles of Na+ coming from NaHCO3
      % sod_hco3 = HCO3;

% Na2SO4 to be added (g)
Na2SO4 = SO4 * Na2SO4_mw;
    % moles of Na+ coming from Na2SO4
      sod_so4 = (SO4)*2; 
    
% KCl to be added (g)
KCl = K * KCl_mw;
     % moles of Cl- coming from KCl
       chlo_k = K;       
    
% MgCl solution to be added (mL)
MgCl2 = (Mg / MgCl2_conc)*1000;
    % moles of Cl- coming from MgCl2
      chlo_mg = (Mg)*2;  

% CaCl2 solution needed (mL)
CaCl2 = (Ca / CaCl2_conc)*1000;    
    % moles of Cl- coming from CaCl2
      chlo_ca = (Ca)*2;   
    


% Na+ to be added from NaCl
Na_rest = Na - (sod_so4+sod_co32);
% Cl- to be added from NaCl
Cl_rest = Cl - (chlo_k) - (chlo_mg) - (chlo_ca);
%         !!!!!! Check  !!!!!!  
   Na_rest - Cl_rest
       
% NaCl to be added (g)
NaCl = Na_rest * NaCl_mw; 


%initial alk estimate
alk =  -10^6*(Na_rest - Cl_rest);    

%initial DIC
K0 = exp((9345.17 / T) - 60.2409 + ...
     23.3585 * log10(T / 100) + ...
     S * (0.023517 - 0.00023656 * T + 0.0047036 * (T / 100)^2));

KH = 10^-1.46
pCO2 = 420*10^-6
DIC = 10^6*m_CO32 + (KH*pCO2)



disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~');   
fprintf('Total Volume to be prepared (L) = %f\n', L);
fprintf('Temperature (C) = %f\n', T);
fprintf('Salinity (g/kg) = %f\n', S);
disp('                                                        ');
fprintf('Na2SO4 to be added (g) = %f\n', Na2SO4); 
fprintf('BH3O3 to be added (g) = %f\n', BH3O3); 
%fprintf('NaHCO3 to be added (g) = %f\n', NaHCO3); 
fprintf('Na2CO32 to be added (g) = %f\n', Na2CO3);
fprintf('KCl to be added (g) = %f\n', KCl); 
fprintf('NaCl to be added (g) = %f\n', NaCl); 
fprintf('MgCl2 solution to be added (mL) = %f\n', MgCl2); 
fprintf('CaCl2 solution to be added (mL) = %f\n', CaCl2); 
fprintf('TA (uEq/kg) %f\n', alk')
fprintf('DIC estimate (uM) %f\n', DIC')
fprintf('Total boron estimate (uM) %f\n', 10^6*B')
%fprintf('Sum of salts (g) = %f\n', Totsalts);
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~');   
