function settings = defaultSettings(user,useSynologyRaw)
    settings.animal = 'FEAA0';
    settings.unit = 'u000';
    settings.experiment = '000';

    settings.username = user;

    if ispc
        settings.filepathSlash = '\';
        settings.rawPath = 'Z:\EphysNew';
        settings.nevFilePrefix = '';
        settings.outShareRoot = 'Z:\EphysNew\outShare';
    else
        settings.filepathSlash = '/';
        settings.rawPath = '/Volumes/NielsenHome/EphysNew';
        settings.nevFilePrefix = [pwd '/'];
        settings.outShareRoot = '/Volumes/NielsenHome/EphysNew/outShare';
    end

    if ~useSynologyRaw
        settings.rawPath = 'raw';
    end
    
    settings.rawAnalyzerPath = [settings.rawPath settings.filepathSlash 'analyzer'];
    settings.rawDataPath = [settings.rawPath settings.filepathSlash 'data'];
    settings.rawlogFilesPath = [settings.rawPath settings.filepathSlash 'log_files'];

    settings.outPath = 'out';
    settings.outSummaryFilePath = [settings.outPath settings.filepathSlash 'summary'];
    settings.outSpikeFilePath = [settings.outPath settings.filepathSlash 'spikeStage'];
    settings.outDataFilePath = [settings.outPath settings.filepathSlash 'dataStage'];
    settings.outPlotsFilePath = [settings.outPath settings.filepathSlash 'plots'];

    settings.outSharePath = [settings.outShareRoot settings.filepathSlash settings.username];
    settings.outShareSummaryFilePath = [settings.outSharePath settings.filepathSlash 'summary'];
    settings.outShareSpikeFilePath = [settings.outSharePath settings.filepathSlash 'spikeStage'];
    settings.outShareDataFilePath = [settings.outSharePath settings.filepathSlash 'dataStage'];
    settings.outSharePlotsFilePath = [settings.outSharePath settings.filepathSlash 'plots'];

    settings.sorterDotSize = 1;

    settings.summaryFileName = 'Ephys Data.xlsx';
end