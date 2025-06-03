%% CMikolaitis @ USA/DISL, 2025
warning('off','MATLAB:table:ModifiedAndSavedVarnames');
warning('off','MATLAB:print:ContentTypeImageSuggested');
clear all
%% Parameters
path = "./Source/";
list = ["USGS*","node*","NOAA*"];
% Mapping Table
idToNode            = table;
idToNode.Name       = ["Orient";"Peconic";"Shelter";"Montauk"];
idToNode.Validation = [01304200; 01304562; 01304650; 8510560];
idToNode.Model      = [48505; 13488; 42773; 80345]; % [48821; 13490; 43039; 80775]
% Offset for differing datums
NGVD29toNAVD88 = -0.95*0.3048; % Based on 2025 Peconic River data
NAVD88toMSL    = -0.101;       % Based on 1983 Epoch Montauk NOAA data
% Validation Window
baseYear = 2021;
baseTime = datetime(baseYear,1,1,'TimeZone','UTC');
t1 = datetime(baseTime,'TimeZone','America/New_York'); 
t2 = datetime(baseYear, 12, 31,'TimeZone','America/New_York');
times = timerange(t1,t2);
% Variable info
varInfo = table;
varInfo.Name = ["Temperature";"Salinity";"Elevation"];
varInfo.Unit = ["°C";"PSU";"m"];
% Figure outputs
saveDir = "./";
if ~exist(saveDir, 'dir')
    mkdir(saveDir);
end
%% File parsing
vTemp = struct(); mTemp = struct();
vSal  = struct(); mSal  = struct();
vElev = struct(); mElev = struct();

boolM = string(idToNode.Model);
boolV = string(idToNode.Validation);
for j = 1:length(list)
    folder = dir(path+list(j));
    for i = 1:length(folder)
        filename = path+folder(i).name;
        parts = split(filename, "_");
        if ~contains(filename,boolM) && ~contains(filename,boolV)
            continue
        end
        % Extract site ID and match to the correct name
        if contains(filename, "USGS") || contains(filename, "NOAA")
            siteID = str2double(parts{2});
            row = find(idToNode.Validation == siteID);
        elseif contains(filename, "node")
            siteID = str2double(parts{2});
            row = find(idToNode.Model == siteID);
        end
        siteName = idToNode.Name{row};

        % Read and process the file
        if contains(filename, 'USGS') % NGVD29
            t = readtable(filename, "FileType", "text", 'Delimiter', '\t');
            time = datetime(t.datetime,'TimeZone','America/New_York');
            out  = table2timetable(t,'RowTimes',time);
            vNames = out.Properties.VariableNames;
            out  = renamevars(out,[string(vNames{7}),string(vNames{9}),string(vNames{5})], ...
                ["temperature","salinity","elevation"]);
            out.elevation = out.elevation*0.3048+NGVD29toNAVD88+NAVD88toMSL;
            vTemp.(char(siteName)) = out(:,7);
            vSal.(char(siteName))  = out(:,9);
            vElev.(char(siteName)) = out(:,5);
    
        elseif contains(filename, 'node') 
            tRaw = readtable(filename, "FileType", "text", 'Delimiter', '\t');
            if contains(filename, 'temperature')
                t = tRaw(tRaw.vgrid_layer == 32, :);  % Filter for top layer
                time = baseTime + seconds(t.time);
                time = datetime(time,'TimeZone','America/New_York');
                out  = table2timetable(t,'RowTimes',time);
                mTemp.(char(siteName)) = out(:,4);
            elseif contains(filename, 'salinity')
                t = tRaw(tRaw.vgrid_layer == 32, :);  % Filter for top layer
                time = baseTime + seconds(t.time);
                time = datetime(time,'TimeZone','America/New_York');
                out  = table2timetable(t,'RowTimes',time);
                mSal.(char(siteName)) = out(:,4);
            elseif contains(filename, 'elevation')
                t = tRaw(tRaw.vgrid_layer == 1, :);
                time = baseTime + seconds(t.time);
                time = datetime(time,'TimeZone','America/New_York');
                out  = table2timetable(t,'RowTimes',time);
                mElev.(char(siteName)) = out(:,4);
            end
        elseif contains(filename, 'NOAA') % MSL
            t    = readmatrix(filename);
            time = datetime(t(:,1), t(:,2), t(:,3), t(:,4), t(:,5), t(:,6),'TimeZone','America/New_York');
            out  = timetable(time,t(:,7));
            if contains(filename, 'temperature')
                out = renamevars(out,'Var1','temperature');
                vTemp.(char(siteName)) = out;
            elseif contains(filename, 'elevation')
                out = renamevars(out,'Var1','elevation');
                vElev.(char(siteName)) = out;
            end
        end
    end
end
%% Generate time filtered paired datasets
pairTemp = caller(vTemp,mTemp);
pairSal  = caller(vSal, mSal);
pairElev = caller(vElev,mElev);
%% Get elevation bias
bias = elevBias(pairElev);
biasCorrection = bias.Montauk;
%% Plot
plotTileComparison(pairTemp,'Temperature','°C', [t1 t2],biasCorrection,saveDir);
plotTileComparison(pairSal, 'Salinity',   'PSU',[t1 t2],biasCorrection,saveDir);
plotTileComparison(pairElev,'Elevation',  'm',  [t1 t2],biasCorrection,saveDir);
%% Kling-Gupta
KGE_Temp = klinggupta(pairTemp);
KGE_Sal  = klinggupta(pairSal);
KGE_Elev = klinggupta(pairElev);
%% Sanity
clearvars -except pair* bias* KGE*
%% Master Func
function pairTable = caller(vStruct,mStruct)
    sites      = fieldnames(vStruct);
    Validation = vStruct.(sites{1});
    Model      = mStruct.(sites{1});
    for i = 2:length(sites)
        Validation = synchronize(Validation,vStruct.(sites{i}),"union",'fillwithmissing');
        Model      = synchronize(Model,mStruct.(sites{i}),"union",'fillwithmissing');
    end
    Validation = renamevars(Validation,Validation.Properties.VariableNames,sites);
    Model      = renamevars(Model,Model.Properties.VariableNames,sites);
    pairTable  = synchronize(Validation,Model,"union",'fillwithmissing');
end
%% Bias Function
function biasTable = elevBias(pairElev)
columns = pairElev.Properties.VariableNames;
    for i = 1:(width(pairElev)/2)
        ia = i+(width(pairElev)/2);
        site = extractBefore(columns{i},"_");
        mData = pairElev.(columns{ia});
        vData = pairElev.(columns{i});
        biasTable.(site) = mean(mData,"omitnan")-mean(vData,"omitnan");
    end
end
%% Plot function
function plotTileComparison(pairTable, variableLabel, yLabel, timeRange, bias, saveDir)
    fig = figure('Name', variableLabel + " Comparison", 'Position', [100, 100, 1200, 600]);
    tiledlayout('flow');
    columns = pairTable.Properties.VariableNames;
    color = [1,0,0,0.25];
    lw = 1;
    sz = 1;
    time = pairTable.Properties.DimensionNames{1};
    for i = 1:(width(pairTable)/2)
        ia   = i+(width(pairTable)/2);
        site = extractBefore(columns{i},"_");
        nexttile;
        scatter(pairTable,time,i, ...
            'filled','MarkerFaceColor','blue','SizeData',sz,'DisplayName', 'Validation');
        hold on
        if contains(variableLabel,'Elevation') && ~contains(site,"Montauk")
            pairTable{:,ia} = pairTable{:,ia}-bias;
            mask = isfinite(pairTable{:,ia});
            plot(pairTable.(time)(mask), pairTable{mask,ia}, ...
                'Color',color,'LineWidth',lw,'DisplayName','Model');
        else
            mask = isfinite(pairTable{:,ia});
            plot(pairTable.(time)(mask), pairTable{mask,ia}, ...
                'Color',color,'LineWidth',lw,'DisplayName','Model');
        end
        title(site + " " + variableLabel);
        xlabel("Time");
        ylabel(yLabel);
        legend('Location', 'best');
        grid on;
        if ~isempty(timeRange)
            xlim(timeRange);
        end
        minY = min(min(pairTable.Variables));
        maxY = max(max(pairTable.Variables));
        ylim([minY maxY]*1.1)
    end
    %filename = fullfile(saveDir, "Comparison_"+variableLabel+".pdf");
    %exportgraphics(fig, filename, 'ContentType', 'vector');
    filename = fullfile(saveDir, "Comparison_"+variableLabel+".fig");
    savefig(fig,filename);
end
%% Kling-Gupta
function output = klinggupta(pairTable)
    output = table('RowNames',{'KGE','r','alpha','beta'});
    columns = pairTable.Properties.VariableNames;
    for i = 1:(width(pairTable)/2)
        ia   = i+(width(pairTable)/2);
        site = extractBefore(columns{i},"_");
        model      = pairTable{:,ia};
        validation = pairTable{:,i};
        mStd  = std(model,"omitnan");
        vStd  = std(validation,"omitnan");
        mMean = mean(model,"omitnan");
        vMean = mean(validation,"omitnan");
        corr  = corrcoef([model,validation],'rows','pairwise'); 
        corr  = corr(1,2);
        pred  = mStd/vStd;   % variability of prediction errors
        bias  = mMean/vMean; % bias
        klinggupta = 1-sqrt(((corr-1)^2) + ((pred-1)^2) + ((bias-1)^2));
        % Save in site specific structs
        output.(site)(1) = klinggupta;
        output.(site)(2) = corr;
        output.(site)(3) = pred;
        output.(site)(4) = bias;
    end
end
