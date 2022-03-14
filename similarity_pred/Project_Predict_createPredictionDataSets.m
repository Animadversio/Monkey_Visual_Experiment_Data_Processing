% function DS = Project_Predict_createPredictionDataSets

myFlags.createDataSet = 0 ;
myFlags.loadDataSet =   1 ;
myFlags.loadEvolvedImages = 0 ; 

net_name = 'alexnet';
clc

if myFlags.createDataSet
    tic
    rootpath = '\\research.files.med.Harvard.edu\neurobio\PonceLab\Stimuli';
    datastore_names_all = {'2022 - Predictions'
        'imagenet-2012\imagenet12\images\train'
        'Caltech256_ObjectCategories'
        'MassiveMemory_Oliva_2008'
        '2021-Natural_Scenes_Dataset\kay-shared1000'
        'Kar_DiCarlo' % make it via Project_General_loadKarDataSet.m
        };
    
    datastores_all = cell(1,length(datastore_names_all));
    for iDS = 1:length(datastore_names_all)
        fprintf('\t loading imagedatastore ''%s''...%.1f min',...
            datastore_names_all{iDS},toc/60);
        datastores_all{iDS} = imageDatastore(fullfile(rootpath,datastore_names_all{iDS}),...
            'IncludeSubfolders',true,'LabelSource','foldernames');
        fprintf(' that took %1.1f min\n',toc/60);
    end
    
    %
    net = alexnet ;
    inputSize = net.Layers(1).InputSize;
    layer = 'fc6';
    
    
    %
    clc; close all;
    allFeatures = [];
    allFeatureNames = [];
    allFeatureDS = cell(length(datastores_all) , 1 ) ;
    max_num_samples = 12e4;
    
    for iDS = 1:length(datastores_all)
        fprintf('%s\n',datastore_names_all{iDS});
        augimdsTrain = augmentedImageDatastore(inputSize(1:2),...
            datastores_all{iDS}, 'ColorPreprocessing','gray2rgb');
        num_files = length(augimdsTrain.Files) ;
        if num_files > max_num_samples
            indices = randsample(num_files,max_num_samples);
            subds = subset(augimdsTrain,indices);
        else
            subds = augimdsTrain;
        end
        
        featuresTrain = activations(net,subds,layer);
        allFeatures= vertcat(allFeatures,squeeze( featuresTrain )' );
        allFeatureNames = vertcat(allFeatureNames,subds.Files);
        allFeatureDS{iDS} = subds;
    end
    
    size(allFeatures)
    figure, imagesc(allFeatures)
    
    
    tic
    fprintf('saving...')
    savefast(fullfile('S:\Data-Ephys-MAT',...
        sprintf('%s_activations_predSet01.mat',net_name)) , 'allFeatures',...
        'allFeatureNames','allFeatureDS');
    fprintf(' that took %1.1f min\n',toc/60);
    winopen('S:\Data-Ephys-MAT')
    
    DS.allFeatures = allFeatures;
    DS.allFeatureNames = allFeatureNames;
    DS.allFeatureDS = allFeatureDS;
    
end % of create datasets

%%

if myFlags.loadDataSet
    DS = load(fullfile('S:\Data-Ephys-MAT',...
        sprintf('%s_activations_predSet01.mat',net_name)));
end

if myFlags.loadEvolvedImages
    %% load up all evolved images
    if ~exist('exampleObject','var')
        exampleObject = matfile('S:\Data-Ephys-MAT\Project_CMA_Alfa_Stats.mat');
        nExp = length(exampleObject.StatsA);
        fieldName = 'StatsA';
        
        exampleObject = matfile('S:\Data-Ephys-MAT\Project_CMA_Beto_Stats.mat');
        nExp = length(exampleObject.StatsB);
        fieldName = 'StatsB';
    end
    
    
    %%
    pics_all = [];
    pics_all_masked = [];
    pics_chan = nan(nExp,1);
    pics_monkey = nan(nExp,1);
    
    figure
    for iExp = 1:nExp
        tmp = exampleObject.(fieldName)(1,iExp) ;
        img =  tmp{1}.picsPerGen_masked(:,:,:,end);
        img2 = tmp{1}.picsPerGen_meanGene(:,:,:,end);
        pics_all= cat(4,pics_all,img2);
        pics_all_masked = cat(4,pics_all_masked,img);
        pics_chan(iExp) = tmp{1}.prefChan;
        pics_monkey(iExp,1) = tmp{1}.monkey(1);
        cla
        montage(pics_all)
    end
    
    
    %%
    path_final = 'N:\Stimuli\2021-Evolutions\Project_CMA_Evolutions_All' ;
    if ~exist(path_final,'dir')
        mkdir(path_final)
    end
    winopen(path_final)
    
    for iExp = 1:nExp
        imwrite(pics_all_masked(:,:,:,iExp),fullfile(path_final,...
            sprintf('%s_%02d_%03d_masked.bmp',pics_monkey(iExp),pics_chan(iExp),iExp )));
        
        imwrite(pics_all(:,:,:,iExp),fullfile(path_final,...
            sprintf('%s_%02d_%03d.bmp',pics_monkey(iExp),pics_chan(iExp),iExp )));
    end
    
    
    
    
%     figure, montage( foo{1}.picsPerGen_masked )
    
    %%
    % testDS = imageDatastore('N:\Stimuli\2019-06-Evolutions\beto-190911b\backup_09_11_2019_14_22_55');
    testDS = imageDatastore('N:\Stimuli\2021-Evolutions\Project_CMA_Evolutions_All');
    
    augimdsTrain = augmentedImageDatastore(inputSize(1:2),...
        testDS, 'ColorPreprocessing','gray2rgb');
    
    isChosenPic = find( contains(augimdsTrain.Files,'B_15_038.bmp') ) ;
    ref_ds = subset(augimdsTrain,isChosenPic);
    pred_images = Project_General_predictImages(ref_ds,allFeatures,allFeatureNames,net,layer);
    
end % of evolved Images