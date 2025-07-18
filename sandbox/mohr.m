
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
agno3_conc = 0.09746; %(mol/L)
agno3_aliq = 0.0025; %(L) amount of agno3 added for each titration
cacl2_conc = 2.3250368913; %*******estimate of cacl2 conc, determined gravimetrically
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


mol_agno3 = (agno3_conc * agno3_aliq); 
mol_cacl2 = mol_agno3/2; %mol_cacl2

dilutemol_cacl2 = mol_cacl2/ 0.13; %amount of cacl2.. mol cacl2 divided by arbitrary dilution factor

L_cacl2 = (dilutemol_cacl2 * 0.1)/cacl2_conc; %amount of stock cacl2 to add to each flask with 100mL DI


disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~');   
fprintf('Mohr Titration recipe \n');
fprintf('Expected aliquot of AgNO3 needed for determination = %f', agno3_aliq*1000);
fprintf('\n');
fprintf('For [AgNO3] (mol/L) = %f\n', agno3_conc);
fprintf('Prepare 5 erlenmeyer flasks each of composition: \n');
fprintf('\n');
fprintf('Stock solution (uL) = %f\n',L_cacl2*10^6);
fprintf('Amount of indicator (5%% K2CrO4) to add (mL) = %f\n',1);
fprintf('DI water (mL) = %f\n',0.100*1000);
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~');


agno3_aliqtrue = [1.8,2.1,1.9,2.0,1.7,1.9,1.95,1.9]; %enter vol of agno3 used for each titration
avg_aliq = mean(agno3_aliqtrue)*10^-3; %avg aliquot in mL

true_Lcacl2 = 40.3*10^-6; %actual amnt of stock dispensed uL
mol_cl = (agno3_conc * avg_aliq); 
mol_ca = (mol_cl/(true_Lcacl2))/2;

fprintf('true molarity of CaCl2 solution (mol/L) = %f\n', mol_ca)
