prompt = "Input Temp (C):";
T = input(prompt);         % Temp in Celsius 

prompt = "Input Salinity (g/kg):";
S = input(prompt);   
      
prompt = "Volume to be prepared (L):";
L = input(prompt);       % Liters of synthetic seawater to be prepared


%Molality (mol/Kg art.SW) of constituent ions in the ASW 
    %  The relative concentration value for ions from 
    %  "Guide to best practices for ocean CO2 measurements"
chloride_rconc = 0.99889;
sulfate_rconc = 0.1400;
sodium_rconc = 0.55661;
magnesium_rconc = 0.06626;
calcium_rconc = 0.02127;
potassium_rconc = 0.0206;
boron_rconc = 0.000232;


m_Cl = ((chloride_rconc/35.45)*(S/0.180655))/10;
m_SO4 = ((sulfate_rconc/96.056)*(S/0.180655))/10;
m_Na = ((sodium_rconc/22.98976928)*(S/0.180655))/10;
m_Mg = ((magnesium_rconc/24.305)*(S/0.180655))/10;
m_Ca = ((calcium_rconc/40.078)*(S/0.180655))/10;
m_K = ((potassium_rconc/39.0983)*(S/0.180655))/10;
m_B = ((boron_rconc/10.81)*(S/0.180655))/10;

m_CO32 = 0.001; 

% Concentration of MgCl2 & CaCl2 solutions (mol/L) to be used
MgCl2_conc = 1.703211;
CaCl2_conc = 0.952078; 

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
    b_oh = B * 3;

    
% Na2CO3 to be added (g)
 Na2CO3 = CO32 * Na2CO3_mw;
%     % moles of Na+ coming from Na2CO3
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
(Na_rest - Cl_rest);
       
% NaCl to be added (g)
NaCl = Na_rest * NaCl_mw; 



%initial alk estimate
alk =  (Na_rest - Cl_rest); 

alk

%initial DIC
KH = 10^-1.46;
pCO2 = 420*10^-6;
DIC = 10^6*m_CO32 + (KH*pCO2);



disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~');   
fprintf('For Non standard seawater of n')
fprintf('Salinity (g/kg) = %f\n', S);
fprintf('Temperature (C) = %f\n', T);
fprintf('Total Volume to be prepared (L) = %f\n', L);
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


% Measured Salts (g)
Na2SO4 = 4.110131325306031;
BH3O3 = 0.02634529307909917;
Na2CO3 = 0.10861463089934326;
NaCl =  29.33378393227003;
KCl = 0.779827399788203;

% Measured solution (mL)
MgCl2 =15.889083744052027;
CaCl2 = 11.067012893336422;

% Concentration of MgCl2 & CaCl2 solutions (mol/L) to be used
MgCl2_molarity = MgCl2_conc;
CaCl2_molarity = CaCl2_conc; 


mol_MgCl2 = MgCl2 * (MgCl2_molarity / 1000);
mol_CaCl2 = CaCl2 * (CaCl2_molarity / 1000);

% Ion Concentration
Mg_conc = ((mol_MgCl2))/(Density*L);
K_conc = ((KCl/KCl_mw))/(Density*L);
SO4_conc = ((Na2SO4/Na2SO4_mw))/(Density*L);
Ca_conc = ((mol_CaCl2))/(Density*L);
Na_conc = (((NaCl/NaCl_mw)+(2*Na2SO4/Na2SO4_mw)))/(Density*L);
Cl_conc = (((NaCl/NaCl_mw)+(2*mol_MgCl2)+(2*mol_CaCl2)+(KCl/KCl_mw)))/(Density*L);
Boron_conc = ((BH3O3/BH3O3_mw))/(Density*L);
CO32_conc = (Na2CO3/Na2CO3_mw)/Density*L;

salinity = (35.453 * Cl_conc) *1.80655;

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
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~');   
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~');   
fprintf('Solution composition: \n');
fprintf('Temperature (C) = %f\n', T);
fprintf('Salinity (g/kg) = %f\n', salinity);
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~');
fprintf('[Na+] (mol/kg H2O) = %f\n', Na_conc/(1-(S/1000)));
fprintf('[K] (mol/kg H2O) = %f\n', K_conc/(1-(S/1000)));
fprintf('[Mg2+] (mol/kg H2O) = %f\n', Mg_conc/(1-(S/1000)));
fprintf('[Ca2+] (mol/kg H2O) = %f\n', Ca_conc/(1-(S/1000)));
fprintf('[Cl-] (mol/kg H2O) = %f\n', Cl_conc/(1-(S/1000)));
fprintf('[SO42-] (mol/kg H2O) = %f\n', SO4_conc/(1-(S/1000)));
fprintf('[CO32-] (mol/kg H2O) = %f\n', CO32_conc/(1-(S/1000)));
fprintf('[Total Boron] (mol/kg H2O) = %f\n', Boron_conc/(1-(S/1000)));
fprintf('[TA] (uEq/kg H2O) = %f\n', -10^6*((Na_conc+2*Mg_conc+2*Ca_conc+K_conc)-(2*CO32_conc+Cl_conc+2*SO4_conc)))/(1-(S/1000));



