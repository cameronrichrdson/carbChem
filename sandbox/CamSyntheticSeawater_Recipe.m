clear variables

% Input 
%==========================================================================
T = 21;      % Temp in Celsius (Lab ~21C most of the time)
S = 35;      % Salinity
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


% complex sw 
% m_Cl = (0.54586);
% m_SO4 = (0.02824);
% m_Na = (0.46906);
% m_Mg = (0.05282);
% m_Ca = (0.01028);
% m_K = (0.01021);
m_CO3 = (0.00026);
m_HCO3 = (0.00177);


% anions = Cl_rest + m_HCO3 + 2*(m_CO3) + 2*(m_SO4);
% cations = Na_rest + m_K + 2*(m_Ca) + 2*(m_Mg);
% 
% 
% cations - anions
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
CO3 = (m_CO3)*(Density)*L;
HCO3 = (m_HCO3)*(Density)*L;

% Molecular mass of salts
Na2SO4_mw = 142.04;
KCl_mw = 74.55;
NaCl_mw = 58.44;
MgCl2_6H2O_mw = 203.3;
CaCl2_2H2O_mw = 147.01;

NaHCO3_mw = 84.00;
Na2CO3_mw = 106.00;



% Na2SO4 to be added (g)
Na2SO4 = SO4 * Na2SO4_mw;
    % moles of Na+ coming from Na2SO4
      sod_so4 = (SO4)*2; 

% Na2CO3 to be added (g)
Na2CO3 = CO3 * Na2CO3_mw;
    % moles of Na+ coming from Na2CO3
      sod_CO32 = (CO3)*2; 

% NaHCO3 to be added (g)
NaHCO3 = HCO3 * NaHCO3_mw;
    % moles of Na+ coming from NaHCO3
      sod_HCO3 = (HCO3);       
    
% KCl to be added (g)
KCl = K * KCl_mw;
     % moles of Cl- coming from KCl
       chlo_k = K;       
    
% MgCl2.6H2O to be added (g)
MgCl2_6H2O = Mg  * MgCl2_6H2O_mw;
    % moles of Cl- coming from MgCl2.6H2O
      chlo_mg = (Mg)*2;  

% CaCl2.2H2O to be added (g)
CaCl2_2H2O = Ca * CaCl2_2H2O_mw; 
    % moles of Cl- coming from CaCl2.2H2O
      chlo_ca = (Ca)*2;   
    
% Na+ to be added from NaCl
Na_rest = Na - (sod_so4) - (sod_HCO3) - (sod_CO32);


% Cl- to be added from NaCl
Cl_rest = Cl - (chlo_k) - (chlo_mg) - (chlo_ca);

% NaCl to be added (g)
NaCl = Na_rest * NaCl_mw;
    

% Total H2O content
% incorporates mole fraction of water in hydrated salts [];'

H2O = (MgCl2_6H2O*0.53)+(CaCl2_2H2O*0.25);

% Total salt content    
Totsalts = (NaCl+Na2CO3+NaHCO3+Na2SO4+KCl+MgCl2_6H2O+CaCl2_2H2O)-(H2O);


disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~');   
fprintf('Total Volume to be prepared (L) = %f\n', L); 
disp('                                                        ');
fprintf('Na2SO4 to be added (g) = %f\n', Na2SO4); 
fprintf('KCl to be added (g) = %f\n', KCl); 
fprintf('NaCl to be added (g) = %f\n', NaCl);
fprintf('Na2CO3 to be added (g) = %f\n', Na2CO3); 
fprintf('NaHCO3 to be added (g) = %f\n', NaHCO3); 
fprintf('MgCl2_6H2O to be added (g) = %f\n', MgCl2_6H2O); 
fprintf('CaCl2_2H2O to be added (g) = %f\n', CaCl2_2H2O); 
disp('                                                        ');
fprintf('Sum of salts (g) = %f\n', Totsalts);
fprintf('Salinity (g/kgsw) = %f\n', Totsalts/1.025);
fprintf('Na_rest - Cl_rest (g) = %f\n', Na_rest - Cl_rest);
disp('                                                        ');
disp('Add these ingridients & bring the solution to the target volume');
disp('with DI-H2O in a Volumetric Flask');
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~');   
