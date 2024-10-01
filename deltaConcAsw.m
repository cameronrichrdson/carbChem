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
m_B = (0.00042); %assuming all boron from one source
m_HCO3 = (0.00177+0.00026); %assuming all alk from HCO3

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


% Molecular mass of ions
Na_mw	= 22.989769;
Mg_mw = 24.305;
Ca_mw = 40.078;
K_mw = 39.0983;
Cl_mw = 35.453;
SO42_mw = 96.06;

% Molecular mass of salts
Na2SO4_mw = 142.04;
KCl_mw = 74.55;
NaCl_mw = 58.44;
MgCl2_mw = 95.20832;
CaCl2_mw = 110.98;
%==========================================================================

% Measured Salts (g)
Na2SO4 = 4.1061;
KCl = 0.7819;
NaCl = 24.5727;
Na2CO3 = 0.1193;



% Measured solution (mL)
MgCl2 = 58.8;
CaCl2 = 10.78;

% Concentration of MgCl2 & CaCl2 solutions (mol/L) to be used
% (arbitarily set to 1 mol/L for now, to simplify calculation)
MgCl2_molarity = 0.92;
CaCl2_molarity = 0.984; 


mol_MgCl2 = MgCl2 * (MgCl2_molarity / 1000);
mol_CaCl2 = CaCl2 * (CaCl2_molarity / 1000);

% Ion Concentration
Mg_conc = ((mol_MgCl2))/(Density*L);
K_conc = ((KCl/KCl_mw))/(Density*L);
SO4_conc = ((Na2SO4/Na2SO4_mw))/(Density*L);
Ca_conc = ((mol_CaCl2))/(Density*L);
Na_conc = (((NaCl/NaCl_mw)+(2*Na2SO4/Na2SO4_mw)))/(Density*L);
Cl_conc = (((NaCl/NaCl_mw)+(2*mol_MgCl2)+(2*mol_CaCl2)+(KCl/KCl_mw)))/(Density*L);



CO32_conc = (Na2CO3/Na2CO3_mw)/Density*L

disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~');   
fprintf('Difference in ideal molality of Asw and prepared molality (mol/kg)'); 
disp('                                                        ');
fprintf('Where delta [n] = [Asw]ideal - [Asw]measured '); 
disp('                                                        ');
disp('                                                        ');
fprintf('delta [Mg2+] (mol/kg) = %f\n', m_Mg - Mg_conc);
fprintf('percent error = %f\n', abs(100-(m_Mg/Mg_conc)*100));
disp('                                                        ');
fprintf('delta [K+] (mol/kg) = %f\n', m_K - K_conc);
fprintf('percent error = %f\n', abs(100-(m_K/K_conc)*100));
disp('                                                        ');
fprintf('delta [SO42-] (mol/kg) = %f\n', m_SO4 - SO4_conc);
fprintf('percent error = %f\n', abs(100-(m_SO4/SO4_conc)*100));
disp('                                                        ');
fprintf('delta [Ca2+] (mol/kg) = %f\n', m_Ca - Ca_conc);
fprintf('percent error = %f\n', abs(100-(m_Ca/Ca_conc)*100));
disp('                                                        ');
fprintf('delta [Na+] (mol/kg) = %f\n', m_Na - Na_conc);
fprintf('percent error = %f\n', abs(100-(m_Na/Na_conc)*100));
disp('                                                        ');
fprintf('delta [Cl-] (mol/kg) = %f\n', m_Cl - Cl_conc);
fprintf('percent error = %f\n', abs(100-(m_Cl/Cl_conc)*100));
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~');   