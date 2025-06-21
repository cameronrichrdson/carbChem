% User inputs
T = input('Input Temp (C): ');
S = input('Input Salinity (g/kg): ');
L = input('Volume to be prepared (L): ');
Mg_multiple = input('Concentration multiplier for Mg: ');
Ca_multiple = input('Concentration multiplier for Ca: ');

% Relative concentration values for ions
chloride_rconc = 0.99889;
sulfate_rconc = 0.1400;
sodium_rconc = 0.55661;
magnesium_rconc = 0.06626;
calcium_rconc = 0.02127;
potassium_rconc = 0.0206;
boron_rconc = 0.000232;

% Molality of constituent ions
m_Cl = (chloride_rconc / 35.45) * (S / 0.180655) / 10;
m_SO4 = (sulfate_rconc / 96.056) * (S / 0.180655) / 10;
m_Na = (sodium_rconc / 22.98976928) * (S / 0.180655) / 10;
m_Mg = (magnesium_rconc / 24.305) * (S / 0.180655) / 10;
m_Ca = (calcium_rconc / 40.078) * (S / 0.180655) / 10;
m_K = (potassium_rconc / 39.0983) * (S / 0.180655) / 10;
m_B = (boron_rconc / 10.81) * (S / 0.180655) / 10;
m_CO32 = 0.001;

% Concentration of solutions (mol/L)
MgCl2_conc = 1.703211;
CaCl2_conc = 0.952078;

% Seawater density calculation
rhow = 999.842594 + 6.793952e-2 * T - 9.095290e-3 * T^2 + 1.001685e-4 * T^3 - 1.120083e-6 * T^4 + 6.536332e-9 * T^5;
A = 8.24493e-1 - 4.0899e-3 * T + 7.6438e-5 * T^2 - 8.2467e-7 * T^3 + 5.3875e-9 * T^4;
B = -5.72466e-3 + 1.0227e-4 * T - 1.6546e-6 * T^2;
C = 4.8314e-4;

% Pure water density at given temperature (kg/m^3)
rho_w = 999.842594 + 6.793952e-2 * T - 9.095290e-3 * T^2 + ...
        1.001685e-4 * T^3 - 1.120083e-6 * T^4 + 6.536332e-9 * T^5;

% Molecular weights (g/mol)
MW_Cl = 35.45;
MW_SO4 = 96.056;
MW_Na = 22.98976928;
MW_Mg = 24.305;
MW_Ca = 40.078;
MW_K = 39.0983;
MW_B = 10.81;
MW_CO32 = 60.01;

% Mass of added solutes (g)
mass_Cl = m_Cl * MW_Cl * L;
mass_SO4 = m_SO4 * MW_SO4 * L;
mass_Na = m_Na * MW_Na * L;
mass_Mg = m_Mg * MW_Mg * L * Mg_multiple;
mass_Ca = m_Ca * MW_Ca * L * Ca_multiple;
mass_K = m_K * MW_K * L;
mass_B = m_B * MW_B * L;
mass_CO32 = m_CO32 * MW_CO32 * L;

% Total added mass (g)
mass_solutes = mass_Cl + mass_SO4 + mass_Na + mass_Mg + mass_Ca + mass_K + mass_B + mass_CO32;

% New density approximation (kg/L)
Density = (rho_w * L + mass_solutes / 1000) / L;



density = rhow + A * S + B * S^(3/2) + C * S^2;
Density_old = density / 1000; % Convert to kg/L

% Moles of constituent ions
Cl = m_Cl * Density * L;
SO4 = m_SO4 * Density * L;
Na = (m_Na * Density * L) - ((m_Mg * Density * L));
Mg = (m_Mg * Density * L) * Mg_multiple;
Ca = (m_Ca * Density * L) * Ca_multiple;
K = m_K * Density * L;
B = m_B * Density * L;
CO32 = m_CO32 * Density * L;

% Molecular masses
Na2SO4_mw = 142.04;
KCl_mw = 74.55;
NaCl_mw = 58.44;
BH3O3_mw = 61.83;
Na2CO3_mw = 105.99;

% Masses to be added (g)
BH3O3 = B * BH3O3_mw;
Na2CO3 = CO32 * Na2CO3_mw;
Na2SO4 = SO4 * Na2SO4_mw;
KCl = K * KCl_mw;

% Solution volumes (mL)
MgCl2 = (Mg / MgCl2_conc) * 1000;
CaCl2 = (Ca / CaCl2_conc) * 1000;

% Adjust NaCl to maintain ionic strength
if Mg_multiple > 1
    Mg_orig = (m_Mg * Density * L);
    Mg_added = Mg_orig * (Mg_multiple - 1);
    NaCl_removed = Mg_added * NaCl_mw; % Remove this much NaCl to compensate
else
    NaCl_removed = 0;
end

% Remaining Na and Cl to balance
Na_rest = Na - (2 * SO4 + 2 * CO32);
Cl_rest = Cl - (K + 2 * Mg + 2 * Ca);
NaCl = Na_rest * NaCl_mw - NaCl_removed; % Adjusted NaCl

% Initial alkalinity and DIC
KH = 10^(-1.46);
pCO2 = 420e-6;
alk = Na_rest - Cl_rest;
DIC = 10^6 * m_CO32 + (KH * pCO2);

I = 1/2 * (Cl + (4 * SO4) + Na + (4 * Mg) + (4 * Ca) + K);

% Display results
fprintf('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n');
fprintf('Ionic strength (g/kg) = %.5f\n', I);
fprintf('Temperature (C) = %.2f\n', T);
fprintf('Total Volume to be prepared (L) = %.2f\n', L);
fprintf('\n');
fprintf('Na2SO4 to be added (g) = %.5f\n', Na2SO4);
fprintf('BH3O3 to be added (g) = %.5f\n', BH3O3);
fprintf('Na2CO32 to be added (g) = %.5f\n', Na2CO3);
fprintf('KCl to be added (g) = %.5f\n', KCl);
fprintf('NaCl to be added (g) = %.5f\n', NaCl);
fprintf('MgCl2 solution to be added (mL) = %.5f\n', MgCl2);
fprintf('CaCl2 solution to be added (mL) = %.5f\n', CaCl2);
fprintf('TA (uEq/kg) = %.5f\n', alk);
fprintf('DIC estimate (uM) = %.5f\n', DIC);
fprintf('Total boron estimate (uM) = %.5f\n', 10^6 * B);
fprintf('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n');