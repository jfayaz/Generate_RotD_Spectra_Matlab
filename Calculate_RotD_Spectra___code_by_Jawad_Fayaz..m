clear all; clc; close all; fclose all; direc = pwd;
%% ------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% author : JAWAD FAYAZ (email: jfayaz@uci.edu), Sarah Azar, Nicholas Luco
% visit: (https://jfayaz.github.io)

% ------------------------------ Instructions ------------------------------------- 
% This code develops the RotD50 Sa and RotD100 Sa Spectra of the Bi-Directional 
% Ground Motion records provided in the 'GMs.mat' file which must be in the current folder. 
% GMs.mat must contain 3 variables
%      'acc1'              --> n x 1 Cell Structure containing GM history in Direction 1 
%      'acc2'              --> n x 1 Cell Structure containing GM history in Direction 2
%      'dt'                --> n x 1 array containing dts for the GM histories
%      
% INPUT:
% This codes provides the option to have 3 different regions of developing the Spectra of ground motions with different period intervals (discretizations)
% The following inputs within the code are required:
% 
%     'Int_T_Reg_1'        --> Period Interval for the first region of the Spectrum 
%     'End_T_Reg_1'        --> Last Period of the first region of the Spectrum (where to end the first region)
%     'Int_T_Reg_2'        --> Period Interval for the second region of the Spectrum 
%     'End_T_Reg_2'        --> Last Period of the second region of the Spectrum (where to end the second region)
%     'Int_T_Reg_3'        --> Period Interval for the third region of the Spectrum 
%     'End_T_Reg_3'        --> Last Period of the third region of the Spectrum (where to end the third region)
%     'Damp'               --> Value of Damping for SDOF system
%     'Plot_Spectra'       --> whether to plot the generated Spectra of the ground motions (options: 'Yes', 'No')    
% 
% 
% OUTPUT:
% The output will be provided in cell structure 'SPECTRA.mat' file which contain RotD50 and RotD100 spectra of the GMs 
% in the same order as given in 'acc1' and 'acc2' variables. The cell structures are developed as:
%     'Periods (secs)' 'RotD50 Sa (g)' 'RotD100 Sa (g)' 
% 
%%-------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ======================== USER INPUTS =============================== %%

% For periods 0 to 'End_T_Reg_1' in an interval of 'Int_T_Reg_1' 
% (Note interval must not be less than 2 decimal points, which means min 'Int_T_Reg_1' = 0.01
Int_T_Reg_1   = 0.05;
End_T_Reg_1   = 1;

% For periods ['End_T_Reg_1'+'Int_T_Reg_2'] to 'End_T_Reg_2' in an interval of 'Int_T_Reg_2'
% (Note interval must not be less than 2 decimal points, which means min 'Int_T_Reg_2' = 0.01
Int_T_Reg_2   = 0.1;
End_T_Reg_2   = 2;

% For periods ['End_T_Reg_2'+'Int_T_Reg_3'] to 'End_T_Reg_3' in an interval of 'Int_T_Reg_3'
% (Note interval must not be less than 2 decimal points, which means min 'Int_T_Reg_3' = 0.01
Int_T_Reg_3   = 0.25;
End_T_Reg_3   = 5;

% Damping value for SDOF
Damp          = 0.05; 

% Plot Spectra  (options: 'Yes' or 'No')
Plot_Spectra  = 'Yes';

%%%%%%================= END OF USER INPUT ========================%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% =============== Generation of Spectra ================
T1 = [Int_T_Reg_1:Int_T_Reg_1:End_T_Reg_1];
T2 = [End_T_Reg_1+Int_T_Reg_2:Int_T_Reg_2:End_T_Reg_2];
T3 = [End_T_Reg_2+Int_T_Reg_3:Int_T_Reg_3:End_T_Reg_3];
Tn = [T1 T2 T3]';  

load(['./GMs.mat'])

for gm = 1:length(acc1)
    dtX = dt(gm,1);
    ugx = acc1{gm,1};
    dtY = dt(gm,1);
    ugy = acc2{gm,1};  
    fprintf('\n\n\nGenerating Spectrum for GM %d\n',gm)
    
    for i = 1:length(Tn)
        T = Tn(i);
        fprintf('\tGenerating Spectrum for GM %d for Period %.2f \n',gm,T)
        Sa50        = RotDxx_Calculations(ugx, ugy, dtX, dtY, 50, T, Damp);
        RotD50(i,1) = Sa50;
        Sa100       = RotDxx_Calculations(ugx, ugy, dtX, dtY, 100, T, Damp);
        RotD100(i,1) = Sa100;
    end
    
    SPECTRA{gm,1} = [Tn, RotD50, RotD100];
    a = num2cell(SPECTRA{gm,1});
    b = cell2struct(a,{'Period','RotD50','RotD100'},2);
    SPECTRA{gm,1} = b;
    
    if strcmpi(Plot_Spectra,'Yes') == 1
        figure(1)
        plot (Tn,RotD50,'o-','linewidth',2,'DisplayName',['GM',num2str(gm)])
        hold on        
        xlabel ('Period (sec)','fontsize',16,'fontWeight','bold')
        ylabel ('RotD50 Sa (g)','fontsize',16,'fontWeight','bold')
        title('RotD50 Spectra','fontsize',16,'fontWeight','bold')
        set(gca,'fontsize',14,'FontName', 'Times New Roman','LineWidth', 1.5)
        grid on; box off
        legend()

        figure(2)
        plot (Tn,RotD100,'o-','linewidth',2,'DisplayName',['GM',num2str(gm)])
        hold on        
        xlabel ('Period (sec)','fontsize',16,'fontWeight','bold')
        ylabel ('RotD100 Sa (g)','fontsize',16,'fontWeight','bold')
        title('RotD100 Spectra','fontsize',16,'fontWeight','bold')
        set(gca,'fontsize',14,'FontName', 'Times New Roman','LineWidth', 1.5)
        grid on; box off
        legend()
    end   
end

save('GM_SPECTRA.mat','SPECTRA')