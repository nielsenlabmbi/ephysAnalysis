settings.animal = 'FEAA0';
settings.unit = 'u000';
settings.experiment = '000';
if ispc
    settings.filepathSlash = '\';
    settings.path = 'Z:\Ephys';
    settings.analyzerPath = 'Z:\Ephys\AnalyzerFiles';
    settings.nevFilePrefix = '';
else
    settings.filepathSlash = '/';
    settings.path = 'raw/data';
    settings.analyzerPath = 'raw/analyzer';
    settings.nevFilePrefix = [pwd '/'];
end

settings.summaryFileFullPath = ['out' settings.filepathSlash 'summary' settings.filepathSlash 'Ephys Data.xlsx'];
settings.spikeFilePath = ['out' settings.filepathSlash 'spikeStage'];
settings.dataFilePath = ['out' settings.filepathSlash 'dataStage'];
settings.sorterDotSize = 1;
