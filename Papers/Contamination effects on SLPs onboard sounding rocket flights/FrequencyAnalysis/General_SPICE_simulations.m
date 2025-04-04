%----------------------------------------------------------------------%
%{
Automated Spice Script

Takes user inputs for plasma parameters, SLP contamination, and Applied Sweep
Profile, and runs SPICE simulation.

Simulations happens in two stages:
1) 60 second simulation following user input at faster time step. This
allows system to reach equilibrium potential as a result of charge
accumulation. This data is not saved.
2) High time resolution simulation for two sweep profiles 
(for 10Hz continuous sweep, this will be 0.2 seconds). This data is saved.


Simulation results for stage 2 are parsed into SLP
up and downsweeps, then stored in a .mat file in specified directory. If no
directory is given, local directory "./OutputFiles/ is created.

Requires the following files are located in ./Accessory_Files
1) general_SPICE_netlist.net
2) LTSPice_call_1.bat
3) LTSpice2Matlab.m (https://www.mathworks.com/matlabcentral/fileexchange/23394-fast-import-of-compressed-binary-raw-files-created-with-ltspice-circuit-simulator)


%}
clear, clc, close all

addpath ./Accessory_Files

% ------------------ User Inputs ------------------------------------------
% default settings will run simulations for the given R, C and plasma
% values over a range of frequencies. A folder /OutputFiles/FrequencyAnalysis will be created
% and each frequency will produce a single file


% Contamination Parameters
Contamination_Resistance = 5e4;     % [Ohms]
Contamination_Capacitance = 1e-6;   % [Farads]

% Plasma Parameters
Ni = 5e11;                          % [/m3]
Te = 1000;                          % [K]
% Applied Sweep Profile Parameters
Sweep_Frequency = logspace(-1,3,10); % [Hz]
SweepProfileType = 4;               % 1 = ion dwell, 2 = electron dwell, 3 = alternating,4 = continuous;
dwellratio = 0 ;                    % dwell:sweep time ratio (must be 0 if continuous sweeping)

% Specify local directory to save output files
directory = strcat('./OutputFiles/FrequencyAnalysis/');
% Files will be save in format:
% "GeneralSim_#e##_R#e#_C#e-#_#Hz_SweepProfileType#_dwellratio#

% Simulation Parameters
max_time_step = 1e-6;               % Defines max time step of SPICE simulation
                                    % ^ Larger values will speed up simulation, but increase the voltage step size

% -------------------------------------------------------------------------

for RESISTOR = Contamination_Resistance 
    for CAPACITOR = Contamination_Capacitance
        clearvars -except RESISTOR CAPACITOR Ni Te Sweep_Frequency SweepProfileType dwellratio directory max_time_step 
        for f = Sweep_Frequency 
    
            simulation_time = 60; % will run 60 seconds at fast time step (hardcoded in file), then single sweep at high time res. (max_time_step)
            dwellflag = SweepProfileType;
            if dwellflag == 3
                numofsweeps = 2; 
            else
                numofsweeps = 1;
            end
            
            
            tsweep = 1/f/2;
            legendstring = {};
            for tdwell_ratio = [dwellratio]
                
                clearvars -except dwellflag f simulation_time tsweep legendstring tdwell_ratio filedescrip numofsweeps dwellratio...
                    RESISTOR CAPACITOR max_time_step Ni Te Sweep_Frequency Sweep ProfileType directory max_time_step SweepProfileType
                
               
                Rc = RESISTOR;
                Cc = CAPACITOR;
                
                numofsweeps = 1; % number of sweeps to run at high time res. Set to 1 for all sims except single FB decays!
                number_highres_profiles = 2; % number of sweeps to run at high time res. Set to 2 for all sims except single FB decays!
                
                % File Description
                filedescrip = strcat('/GeneralSim_',sprintf('%0.0e',Ni),'_R',sprintf('%0.3e',Rc),'_C',sprintf('%0.3e',Cc),'_',num2str(f),'Hz_','SweepProfileType',num2str(dwellflag),'_dwellratio',num2str(tdwell_ratio));
           
                
                % find number of profiles
                tdwellval = tsweep *2* tdwell_ratio; 
                if dwellflag == 4 || dwellflag == 3
                    numberofprofiles = round(simulation_time/ (tsweep*2*numofsweeps + tdwellval));
                else
                    numberofprofiles = round(simulation_time/ (tsweep*2*numofsweeps + tdwellval));
                end

                % Initialize Arrays, Start Counters
                tupcut = [];
                tcutdown = [];
                upcount = 1;
                downcount = 1;


                % Create Piecewise Linear Function for SPICE netlist
                % 60 second PWL following user input, allows system to reach
                % equilibrium potential
                
                if mod(numberofprofiles,2) == 0 % 
                    numberofprofiles = numberofprofiles + 1; % force it to end odd
                end
                
                if dwellflag == 1
                    for i = 1:1:numberofprofiles
                        if i ==1
                            pwl = '0, -5,';
                            t = 0;
                            for j = 1:1:2*numofsweeps+1
                                if j == 2*numofsweeps+1 % add dwell
                                    t = t + tdwellval;
                                    pwl = strcat(pwl,num2str(t),', -5,');
                                else
                                    if mod(j,2) == 1 %odd
                                        tupcut(upcount,1) = t;
                                        t = t + tsweep;
                                        pwl = strcat(pwl,num2str(t),', 5,');
                                        tupcut(upcount,2) = t;
                                        upcount = upcount + 1;
                                    elseif mod(j,2) == 0 % even
                                        tdowncut(downcount,1) = t;
                                        t = t + tsweep;
                                        tdowncut(downcount,2) = t;
                                        downcount = downcount + 1;
                                        pwl = strcat(pwl,num2str(t),', -5,');
                                    end
                                    
                                end
                            end
                        else
                            for j = 1:1:2*numofsweeps+1
                                if j == 2*numofsweeps+1 % add dwell
                                    t = t + tdwellval;
            
                                    pwl = strcat(pwl,num2str(t),', -5,');
                                else
                                    if mod(j,2) == 1 %odd
                                        tupcut(upcount,1) = t;
                                        t = t + tsweep;
                                        tupcut(upcount,2) = t;
                                        upcount = upcount + 1;
                                        pwl = strcat(pwl,num2str(t),', 5,');
                                    elseif mod(j,2) == 0 % even
                                        tdowncut(downcount,1) = t;
                                        t = t + tsweep;
                                        tdowncut(downcount,2) = t;
                                        downcount = downcount + 1;
                                        pwl = strcat(pwl,num2str(t),', -5,');
                                    end
                                    
                                end
                            end
            
            
                        end
                    end
                elseif dwellflag == 2
                    for i = 1:1:numberofprofiles
                        if i ==1
                            pwl = '0, 5,';
                            t = 0;
                            for j = 1:1:2*numofsweeps+1
                                if j == 2*numofsweeps+1 % add dwell
                                    t = t + tdwellval;
                                    pwl = strcat(pwl,num2str(t),', 5,');
                                else
                                    if mod(j,2) == 1 %odd
                                        tdowncut(downcount,1) = t;
                                        t = t + tsweep;
                                        tdowncut(downcount,2) = t;
                                        downcount = downcount + 1;
                                        pwl = strcat(pwl,num2str(t),', -5,');
                                    elseif mod(j,2) == 0 % even
                                        tupcut(upcount,1) = t;
                                        t = t + tsweep;
                                        tupcut(upcount,2) = t;
                                        upcount = upcount + 1;
                                        pwl = strcat(pwl,num2str(t),', 5,');
                                    end
                                    
                                end
                            end
                        else
                            for j = 1:1:2*numofsweeps+1
                                if j == 2*numofsweeps+1 % add dwell
                                    t = t + tdwellval;
                                    pwl = strcat(pwl,num2str(t),', 5,');
                                else
                                    if mod(j,2) == 1 %odd
                                        tdowncut(downcount,1) = t;
                                        t = t + tsweep;
                                        tdowncut(downcount,2) = t;
                                        downcount = downcount + 1;
                                        pwl = strcat(pwl,num2str(t),', -5,');
                                    elseif mod(j,2) == 0 % even
                                        tupcut(upcount,1) = t;
                                        t = t + tsweep;
                                        tupcut(upcount,2) = t;
                                        upcount = upcount + 1;
                                        pwl = strcat(pwl,num2str(t),', 5,');
                                    end
                                    
                                end
                            end
            
            
                        end
                    end
                elseif dwellflag == 3 % alternating profile, will always start low
                    for i = 1:1:numberofprofiles
                        if i ==1
                            pwl = '0, -5,';
                            t = 0;
                            for j = 1:1:2*numofsweeps+1
                                if j == 2*numofsweeps+1 % add dwell
                                    t = t + tdwellval;
                                    pwl = strcat(pwl,num2str(t),', -5,');
                                else
                                    if mod(j,2) == 1 %odd
                                        tupcut(upcount,1) = t;
                                        t = t + tsweep;
                                        tupcut(upcount,2) = t;
                                        upcount = upcount + 1;
                                        pwl = strcat(pwl,num2str(t),', 5,');
                                    elseif mod(j,2) == 0 % even
                                        tdowncut(downcount,1) = t;
                                        t = t + tsweep;
                                        tdowncut(downcount,2) = t;
                                        downcount = downcount + 1;
                                        pwl = strcat(pwl,num2str(t),', -5,');
                                    end
                                    
                                end
                            end
                        elseif mod(i,2) == 0 %even , add high dwell
                            for j = 1:1:2*numofsweeps+1
                                if j == 2*numofsweeps +1% add dwell
                                    t = t + tsweep;
                                    pwl = strcat(pwl,num2str(t),', 5,');
                                    t = t + tdwellval;
                                    pwl = strcat(pwl,num2str(t),', 5,');
                                else
                                    if mod(j,2) == 1 %odd
                                        tupcut(upcount,1) = t;
                                        t = t + tsweep;
                                        tupcut(upcount,2) = t;
                                        upcount = upcount + 1;
                                        pwl = strcat(pwl,num2str(t),', 5,');
                                    elseif mod(j,2) == 0 % even
                                        tdowncut(downcount,1) = t;
                                        t = t + tsweep;
                                        tdowncut(downcount,2) = t;
                                        downcount = downcount + 1;
                                        pwl = strcat(pwl,num2str(t),', -5,');
                                    end
                                    
                                end
                            end
                        else % odd, add low dwell
                            for j = 1:1:2*numofsweeps+1
                                if j == 2*numofsweeps+1 % add dwell
                                    t = t + tsweep;
                                    pwl = strcat(pwl,num2str(t),', -5,');
                                    t = t + tdwellval;
                                    pwl = strcat(pwl,num2str(t),', -5,');
                                else
                                    if mod(j,2) == 1 %odd
                                        tdowncut(downcount,1) = t;
                                        t = t + tsweep;
                                        tdowncut(downcount,2) = t;
                                        downcount = downcount + 1;
                                        pwl = strcat(pwl,num2str(t),', -5,');
                                    elseif mod(j,2) == 0 % even
                                        tupcut(upcount,1) = t;
                                        t = t + tsweep;
                                        tupcut(upcount,2) = t;
                                        upcount = upcount + 1;
                                        pwl = strcat(pwl,num2str(t),', 5,');
                                    end
                                    
                                end
                            end
                        end
                    end
                else %dwellflag == 4
                    pwl = '0, -5,';
                    t = 0;
                    for i = 1:1:numberofprofiles
                        for j = 1:1:2
                            if j == 1
                                tupcut(upcount,1) = t;
                                t = t + tsweep;
                                tupcut(upcount,2) = t;
                                upcount = upcount + 1;
                                pwl = strcat(pwl,num2str(t),', 5,');
                            elseif j ==2 % even
                                tdowncut(downcount,1) = t;
                                t = t + tsweep;
                                tdowncut(downcount,2) = t;
                                downcount = downcount + 1;
                                pwl = strcat(pwl,num2str(t),', -5,');
                            end
                        end
                    end
                end

                % Save PWL
                pwl_long = pwl;
                t_long = t;
                clearvars tupcut tdowncut pwl t upcount downcount

                % Create Piecewise Linear Function for SPICE netlist
                % Short duration, high time res  PWL following user input,
                % Time duration = time to sweep twice
                tupcut = [];
                tcutdown = [];
                upcount = 1;
                downcount = 1;
                simulation_time = 1/f*2; % two sweeps worth
                
                numberofprofiles = number_highres_profiles ;
                simulation_time = numberofprofiles *(tsweep*2*numofsweeps + tdwellval);
                % -------------- now, for short sim -----------------
                if dwellflag == 1
                    for i = 1:1:numberofprofiles
                        if i ==1
                            pwl = '0, -5,';
                            t = 0;
                            for j = 1:1:2*numofsweeps+1
                                if j == 2*numofsweeps+1 % add dwell
                                    t = t + tdwellval;
                                    pwl = strcat(pwl,num2str(t),', -5,');
                                else
                                    if mod(j,2) == 1 %odd
                                        tupcut(upcount,1) = t;
                                        t = t + tsweep;
                                        pwl = strcat(pwl,num2str(t),', 5,');
                                        tupcut(upcount,2) = t;
                                        upcount = upcount + 1;
                                    elseif mod(j,2) == 0 % even
                                        tdowncut(downcount,1) = t;
                                        t = t + tsweep;
                                        tdowncut(downcount,2) = t;
                                        downcount = downcount + 1;
                                        pwl = strcat(pwl,num2str(t),', -5,');
                                    end
                                    
                                end
                            end
                        else
                            for j = 1:1:2*numofsweeps+1
                                if j == 2*numofsweeps+1 % add dwell
                                    t = t + tdwellval;
            
                                    pwl = strcat(pwl,num2str(t),', -5,');
                                else
                                    if mod(j,2) == 1 %odd
                                        tupcut(upcount,1) = t;
                                        t = t + tsweep;
                                        tupcut(upcount,2) = t;
                                        upcount = upcount + 1;
                                        pwl = strcat(pwl,num2str(t),', 5,');
                                    elseif mod(j,2) == 0 % even
                                        tdowncut(downcount,1) = t;
                                        t = t + tsweep;
                                        tdowncut(downcount,2) = t;
                                        downcount = downcount + 1;
                                        pwl = strcat(pwl,num2str(t),', -5,');
                                    end
                                    
                                end
                            end
            
            
                        end
                    end
                elseif dwellflag == 2
                    for i = 1:1:numberofprofiles
                        if i ==1
                            pwl = '0, 5,';
                            t = 0;
                            for j = 1:1:2*numofsweeps+1
                                if j == 2*numofsweeps+1 % add dwell
                                    t = t + tdwellval;
                                    pwl = strcat(pwl,num2str(t),', 5,');
                                else
                                    if mod(j,2) == 1 %odd
                                        tdowncut(downcount,1) = t;
                                        t = t + tsweep;
                                        tdowncut(downcount,2) = t;
                                        downcount = downcount + 1;
                                        pwl = strcat(pwl,num2str(t),', -5,');
                                    elseif mod(j,2) == 0 % even
                                        tupcut(upcount,1) = t;
                                        t = t + tsweep;
                                        tupcut(upcount,2) = t;
                                        upcount = upcount + 1;
                                        pwl = strcat(pwl,num2str(t),', 5,');
                                    end
                                    
                                end
                            end
                        else
                            for j = 1:1:2*numofsweeps+1
                                if j == 2*numofsweeps+1 % add dwell
                                    t = t + tdwellval;
                                    pwl = strcat(pwl,num2str(t),', 5,');
                                else
                                    if mod(j,2) == 1 %odd
                                        tdowncut(downcount,1) = t;
                                        t = t + tsweep;
                                        tdowncut(downcount,2) = t;
                                        downcount = downcount + 1;
                                        pwl = strcat(pwl,num2str(t),', -5,');
                                    elseif mod(j,2) == 0 % even
                                        tupcut(upcount,1) = t;
                                        t = t + tsweep;
                                        tupcut(upcount,2) = t;
                                        upcount = upcount + 1;
                                        pwl = strcat(pwl,num2str(t),', 5,');
                                    end
                                    
                                end
                            end
            
            
                        end
                    end
                else % alternating profile, will always start low
                    for i = 1:1:numberofprofiles
                        if i ==1
                            pwl = '0, -5,';
                            t = 0;
                            for j = 1:1:2*numofsweeps+1
                                if j == 2*numofsweeps+1 % add dwell
                                    t = t + tdwellval;
                                    pwl = strcat(pwl,num2str(t),', -5,');
                                else
                                    if mod(j,2) == 1 %odd
                                        tupcut(upcount,1) = t;
                                        t = t + tsweep;
                                        tupcut(upcount,2) = t;
                                        upcount = upcount + 1;
                                        pwl = strcat(pwl,num2str(t),', 5,');
                                    elseif mod(j,2) == 0 % even
                                        tdowncut(downcount,1) = t;
                                        t = t + tsweep;
                                        tdowncut(downcount,2) = t;
                                        downcount = downcount + 1;
                                        pwl = strcat(pwl,num2str(t),', -5,');
                                    end
                                    
                                end
                            end
                        elseif mod(i,2) == 0 %even , add high dwell
                            for j = 1:1:2*numofsweeps+1
                                if j == 2*numofsweeps +1% add dwell
                                    t = t + tsweep;
                                    pwl = strcat(pwl,num2str(t),', 5,');
                                    t = t + tdwellval;
                                    pwl = strcat(pwl,num2str(t),', 5,');
                                else
                                    if mod(j,2) == 1 %odd
                                        tupcut(upcount,1) = t;
                                        t = t + tsweep;
                                        tupcut(upcount,2) = t;
                                        upcount = upcount + 1;
                                        pwl = strcat(pwl,num2str(t),', 5,');
                                    elseif mod(j,2) == 0 % even
                                        tdowncut(downcount,1) = t;
                                        t = t + tsweep;
                                        tdowncut(downcount,2) = t;
                                        downcount = downcount + 1;
                                        pwl = strcat(pwl,num2str(t),', -5,');
                                    end
                                    
                                end
                            end
                        else % odd, add low dwell
                            for j = 1:1:2*numofsweeps+1
                                if j == 2*numofsweeps+1 % add dwell
                                    t = t + tsweep;
                                    pwl = strcat(pwl,num2str(t),', -5,');
                                    t = t + tdwellval;
                                    pwl = strcat(pwl,num2str(t),', -5,');
                                else
                                    if mod(j,2) == 1 %odd
                                        tdowncut(downcount,1) = t;
                                        t = t + tsweep;
                                        tdowncut(downcount,2) = t;
                                        downcount = downcount + 1;
                                        pwl = strcat(pwl,num2str(t),', -5,');
                                    elseif mod(j,2) == 0 % even
                                        tupcut(upcount,1) = t;
                                        t = t + tsweep;
                                        tupcut(upcount,2) = t;
                                        upcount = upcount + 1;
                                        pwl = strcat(pwl,num2str(t),', 5,');
                                    end
                                    
                                end
                            end
                        end
                    end
            
                end
               
                % -------------------------------------------------------%
                %% Call SPICE netlists and run

                % Feed user inputs, piecewise linear function and simulation time to netlist
                % Run spice in background (This can take a while if running at
                % very small timestep to get good resolution around
                % floating potentials)
                
                rawdata = Variable_netlist_Fast(pwl_long,t_long,Ni,Te,Rc,Cc);
                
                % extract four primary voltage nodes
                l = length(rawdata.time_vect);
                ColHeaders = rawdata.variable_name_list;
                Data(:,1) = rawdata.time_vect(:);
                for i = 2:1:length(ColHeaders) 
                    Data(:,i) = rawdata.variable_mat(i-1,:);  
                end
                dat = Data(2:end,:);
                ColHeaders = rawdata.variable_name_list;
                % Find indices
        
                [garb tf1] = max(strcmp(ColHeaders,'V(1)'));
                [garb tf3] = max(strcmp(ColHeaders,'V(3)'));
                [garb tf4] = max(strcmp(ColHeaders,'V(4)'));

                V1end = dat(end,tf1+1);
                V3end = dat(end,tf3+1);
                V4end = dat(end,tf4+1);

                
                ICs = {num2str(V1end), num2str(V3end), num2str(V4end)};
            
                if numberofprofiles == 2
                tend_short = 1/f*numofsweeps+1/f*.1;
                else
                    tend_short = simulation_time;
                end
                clearvars Data ColHeaders l dat 
                rawdata = Variable_netlist(pwl,tend_short,Ni,Te,Rc,Cc,ICs,max_time_step); % one sweep, plus 10% buffer time
        

                %% Take SPICE output, organize and save it to .mat
            
                l = length(rawdata.time_vect);
                ColHeaders = rawdata.variable_name_list;
                Data(:,1) = rawdata.time_vect(:);
                for i = 2:1:length(ColHeaders)             
                    Data(:,i) = rawdata.variable_mat(i-1,:);  
                end
            
                dat = Data(2:end,:);
                ColHeaders = rawdata.variable_name_list;
                % Find indices
                [garb tfc] = max(strcmp(ColHeaders,'I(Cc)'));
                [garb tfcr] = max(strcmp(ColHeaders,'I(Rc)'));
                [garb tf4] = max(strcmp(ColHeaders,'V(4)'));
                [garb tf1] = max(strcmp(ColHeaders,'V(1)'));
                [garb tf33] = max(strcmp(ColHeaders,'V(3)'));
                [garb tfr] = max(strcmp(ColHeaders,'Ix(langmuir:A)'));
            
                V = dat(:,tf1+1) - dat(:,tf4+1) ;
                Vsc = dat(:,tf4+1);

                I = dat(:,tfr+1);
                time = dat(:,1);
                V1 = dat(:,tf1+1);
                V3 = dat(:,tf33+1);
                Vcap = V3-V1;


                % Split into up,down, and dwell

                % Initialize Arrays
                Vup = {};
                Vup_sc = {};
                Vup_fpp = {};
                Iup = {};
                tup = {};
                Vdown = {};
                Vdown_sc = {};
                Vdown_fpp ={};
                tdown = {};
                Idown = {};
                Vdwell = {};
                Idwell = {};
                tdwell = {};
                Vcap_up = {};
                Vcap_down = {};
                for i = 1:1:length(tupcut(:,1))
                    tstart = tupcut(i,1);
                    tend = tupcut(i,2);
            
                    [garb loc_prev] = min(abs(time - tstart));
                    [garb loc] = min(abs(time - tend));
            
                    Vup{i} = V(loc_prev:loc);
                    Vup_sc{i} = Vsc(loc_prev:loc);
                    Iup{i} = I(loc_prev:loc);
                    tup{i} = time(loc_prev:loc);
                    Vcap_up{i} = Vcap(loc_prev:loc);
        
                    % Extract dwell
                    if dwellflag == 2 && i ~= length(tdowncut(:,1)) && tend ~= tdowncut(i+1,1)
                        [garb locdwellend] = min(abs(time - tdowncut(i+1,1)));
                        Vdwell{i} = V(loc:locdwellend);
                        V3dwell{i} = V3(loc:locdwellend);
                        Idwell{i} = I(loc:locdwellend);
                        tdwell{i} = time(loc:locdwellend);
                    end
                end
                
                for i = 1:1:length(tdowncut(:,1))
                    tstart = tdowncut(i,1);
                    tend = tdowncut(i,2);
            
                    [garb loc_prev] = min(abs(time - tstart));
                    [garb loc] = min(abs(time - tend));
            
                    Vdown{i} = V(loc_prev:loc);
                    Vdown_sc{i} = Vsc(loc_prev:loc);
                    Idown{i} = I(loc_prev:loc);
                    tdown{i} = time(loc_prev:loc);
                    Vcap_down{i} = Vcap(loc_prev:loc);
            
                    % Extract dwell
                    if dwellflag ==1 && i ~= length(tupcut(:,1)) && tend ~= tupcut(i+1,1)
                        [garb locdwellend] = min(abs(time - tupcut(i+1,1)));
                        Vdwell{i} = V(loc:locdwellend);
                        V3dwell{i} = V3(loc:locdwellend);
                        Idwell{i} = I(loc:locdwellend);
                        tdwell{i} = time(loc:locdwellend);
                    end
                end
                % Save data as .mat file in one structure
                DAT = struct();
                DAT.Vup = Vup;
                DAT.Vup_sc = Vup_sc;
                DAT.Vcap_up = Vcap_up;
                DAT.Vcap_down = Vcap_down;
                DAT.Vup_fpp = Vup_fpp;
                DAT.Tup = tup;
                DAT.Iup = Iup;
                DAT.Vdown = Vdown;
                DAT.Vdown_sc = Vdown_sc;
                DAT.Vdown_fpp = Vdown_fpp;
                DAT.Tdown = tdown;
                DAT.Idown = Idown;
                if dwellflag ==1 || dwellflag == 2
                    DAT.Vdwell = Vdwell;  
                    DAT.V3dwell = V3dwell;
                    DAT.Idwell = Idwell;
                    DAT.Tdwell = tdwell;
                end
                DAT.ColHeaders = ColHeaders;
                DAT.dat = dat;
                
                % Check if folder exists... if not, create it
                subfolder = strcat(directory,'/FourRegimes_5e+11_R',sprintf('%0.3e',Rc),'_C',sprintf('%0.3e',Cc),'/');
                if ~exist(subfolder,'dir')
                    mkdir (subfolder)
                end
                filename = strcat(subfolder,filedescrip,'.mat');
                save(filename,'DAT')        
            end % dwellratio
        end % f
    end % CAPACITOR
end % RESISTOR




