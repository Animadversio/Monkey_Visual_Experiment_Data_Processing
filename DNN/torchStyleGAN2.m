classdef torchStyleGAN2
   % Usage: 
   % G = torchStyleGAN2("model.ckpt-533504.pt");
   % matimg = visualize_codes(G,randn(4,512));
   %   > Elapsed time is 0.820856 seconds.
   % figure;montage(matimg)
   % 
   % First Time Setup Instruction:
   % 
   % # Get a python env with suitable pytorch version
   % Note for Python environment, Pytorch 1.3.1 is recommended. 
   % Pytorch 1.1.0 cannot compile the StyleGAN operators... >1.4.0 will not
   % work with matlab. 
   % assume the name of it is [envname]
   % 
   % # Compile model code (CUDA, C++ code of operators)
   % 1. Open an Anaconda Prompt, activate that env: conda activate [envname]
   % 2. cd to ...\stylegan2-pytorch\op
   % 3. Compile the source code and install the operators as py package `python setup.py install` 
   % 3.1. This is really easy to fail spectacularly... So finger crossed. 
   % 3.2. Make sure you have a Visual Studio 2017 or 2019. VS2015 compiler
   %      will not work
   % 4. Check compilation: cd to ...\stylegan2-pytorch
   % python
   % > from model import Generator
   % 
   % If it succeed you are done on python side! The compilation is successful, 
   % you are half way there. If you have some checkpoints stored there you
   % can generate some pretty images in command line through
   % > 
   % 
   % # Setup Python env in Matlab
   % 
   % For ML2A machine, setup the python env before first time use like this
   % >  setenv('path',['C:\Anaconda3\envs\torch\Library\bin;', getenv('path')]);
   % >  pyenv("Version","C:\Anaconda3\envs\torch\python.exe");
   properties
       Generator
       config
       mean_latent
       random
       Wspace
       batchsize
   end
   methods
   function G = torchStyleGAN2(ckpt, config)
       % Keys have to be a cell array of char strs, values have to be cell
       % array. 
        configMap = containers.Map({'stylegan2-cat-config-f.pt', ...
                                    'ffhq-256-config-e-003810.pt', ...
                                    'ffhq-512-avg-tpurun1.pt', ...
                                    'stylegan2-ffhq-config-f.pt', ...
                                    '2020-01-11-skylion-stylegan2-animeportraits.pt', ...
                                    'stylegan2-car-config-f.pt', ...
                                    'model.ckpt-533504.pt'}, ...
                {struct("size", int32(256), "n_mlp", int32(8), "channel_multiplier", int32(2), "latent", int32(512)),...
                struct("size", int32(256), "n_mlp", int32(8), "channel_multiplier", int32(1), "latent", int32(512)),...
                struct("size", int32(512), "n_mlp", int32(8), "channel_multiplier", int32(2), "latent", int32(512)),...
                struct("size", int32(1024), "n_mlp", int32(8), "channel_multiplier", int32(2), "latent", int32(512)),...
                struct("size", int32(512), "n_mlp", int32(8), "channel_multiplier", int32(2), "latent", int32(512)),...
                struct("size", int32(512), "n_mlp", int32(8), "channel_multiplier", int32(2), "latent", int32(512)),...
                struct("size", int32(512), "n_mlp", int32(8), "channel_multiplier", int32(2), "latent", int32(512))});
       if nargin == 0
           % config = struct("latent",int32(512),"n_mlp",int32(8),"channel_multiplier",int32(2));
           ckpt = "stylegan2-ffhq-config-f.pt";
           config = configMap(ckpt);
       elseif nargin == 1
           config = configMap(ckpt); % struct("latent",int32(512),"n_mlp",int32(8),"channel_multiplier",int32(2));
       end
       % Use the torch 1.3.x and the stylegan2 package like below.
       py.importlib.import_module('torch');
       syspath = py.sys.path(); % add the official stylegan2 repo. 
       switch getenv('COMPUTERNAME')
           case 'DESKTOP-9DDE2RH' % Office 3 Binxu's 
            syspath.append("D:\Github\stylegan2-pytorch"); 
            savedir = "D:\Github\stylegan2-pytorch\checkpoint";  
           case 'DESKTOP-MENSD6S' % Binxu's home work station
            syspath.append("E:\DL_Projects\Vision\stylegan2-pytorch"); 
            savedir = "E:\DL_Projects\Vision\stylegan2-pytorch\checkpoint"; 
           case 'PONCELAB-ML2A' % MLa machine 
            syspath.append("C:\Users\Poncelab-ML2a\Documents\Python\stylegan2-pytorch"); 
            savedir = "C:\Users\Poncelab-ML2a\Documents\Python\stylegan2-pytorch\checkpoint"; 
           case 'PONCELAB-ML2B' % MLb machine 
            % savedir = "C:\Users\Ponce lab\Documents\Python\pytorch-pretrained-BigGAN\weights";
           otherwise
            % savedir = "C:\Users\Poncelab-ML2a\Documents\Python\pytorch-pretrained-BigGAN\weights";
        end
       py.importlib.import_module('model');
       G.Generator = py.model.Generator(config.size, config.latent, config.n_mlp, config.channel_multiplier);
       SD = py.torch.load(fullfile(savedir, ckpt));
       G.Generator.load_state_dict(SD.get('g_ema'));
       G.Generator.to('cuda');G.Generator.eval();
       py.torch.set_grad_enabled(false);
       G.config = config;
       G.mean_latent = G.Generator.mean_latent(int32(4096)); % estimate a mean W latent code by averaging. 
       G.random = false;
       G.Wspace = false; 
       % Tune the batch size depending on the resolution of generated
       % image, secure the memory will not overflow in CUDA. 
       if config.size == 256
           G.batchsize = 40;
       elseif config.size == 512
           G.batchsize = 15;
       elseif config.size == 1024
           G.batchsize = 6;
       end
   end
   function W = style_map(G, Z, truncation)
    if nargin == 2, truncation=1; end
    Wtsr = G.Generator.get_latent(py.torch.tensor(py.numpy.array(Z)).float().cuda());
    if truncation < 1 % truncation shrink the W vector towards center of W distribution
        Wtsr = truncation * (Wtsr - G.mean_latent) + G.mean_latent;
    end
    W = Wtsr.detach.cpu().numpy().single;
   end

   function matimgs = visualize(G, style, truncation)
       if nargin == 2, truncation=0.7;end
       meanlatent = py.None;
       if truncation < 1, meanlatent = G.mean_latent; end% G.Generator.mean_latent(int32(4096));
%        G.batchsize = 10; % decide the batch size in initialization
       samplen = size(style, 1);
       tic
       csr = 1; matimgs = [];
       while csr <= samplen
       cnd = min(samplen,csr+G.batchsize);
       imgs = G.Generator(py.torch.tensor(py.numpy.array(style(csr:cnd,:))).view(py.tuple(int32([1,cnd-csr+1,size(style,2)]))).float().cuda(),...
           pyargs("truncation",truncation,"truncation_latent",meanlatent,"input_is_latent",G.Wspace,"randomize_noise",G.random));
       imgs = imgs{1}; % discard the 2nd term in tuple
       matimg = imgs.detach.cpu().numpy().single; % note range in -1, 1
       matimg = permute(clip((matimg + 1) / 2, 0, 1), [3,4,2,1]);
       matimgs = cat(4, matimgs, matimg);
       csr = cnd + 1;
       end
       toc
       
   end
   function matimg = visualize_singleBatch(G, style, truncation)
       if nargin == 2, truncation=0.7;end
       meanlatent = py.None;
       if truncation < 1, meanlatent = G.mean_latent; end% G.Generator.mean_latent(int32(4096));
       imgs = G.Generator(py.torch.tensor(py.numpy.array(style)).view(py.tuple(int32([1,size(style)]))).float().cuda(),...
           pyargs("truncation",truncation,"truncation_latent",meanlatent,"input_is_latent",G.Wspace,"randomize_noise",G.random));
       imgs = imgs{1}; % discard the 2nd term in tuple
       matimg = imgs.detach.cpu().numpy().single; % note range in -1, 1
       matimg = permute(clip((matimg + 1) / 2, 0, 1),[3,4,2,1]);
   end
   end
end