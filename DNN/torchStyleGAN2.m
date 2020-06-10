classdef torchStyleGAN2
   % Usage 
   % G = torchStyleGAN2("model.ckpt-533504.pt");
   % matimg = visualize_codes(G,randn(4,512));
   %   > Elapsed time is 0.820856 seconds.
   % figure;montage(matimg)
   properties
       Generator
       config
   end
   methods
   function G = torchStyleGAN2(ckpt, config)
       if nargin == 0
           config = struct("latent",int32(512),"n_mlp",int32(8),"channel_multiplier",int32(2));
           ckpt = "stylegan2-ffhq-config-f.pt";
       elseif nargin == 1
           config = struct("latent",int32(512),"n_mlp",int32(8),"channel_multiplier",int32(2));
       end
       switch ckpt
           case "stylegan2-ffhq-config-f.pt"
               config.size = int32(1024);
           case {"model.ckpt-533504.pt", "2020-01-11-skylion-stylegan2-animeportraits.pt"}
               config.size = int32(512);
           otherwise
               config.size = int32(512);
       end
       % Use the torch 1.3.x and the stylegan2 package like below.
       py.importlib.import_module('torch');
       syspath = py.sys.path(); % add the official stylegan2 repo. 
       syspath.append("E:\DL_Projects\Vision\stylegan2-pytorch");
       py.importlib.import_module('model');
       G.Generator = py.model.Generator(config.size, config.latent, config.n_mlp, config.channel_multiplier);
       savedir = "E:\DL_Projects\Vision\stylegan2-pytorch\checkpoint";
       SD = py.torch.load(fullfile(savedir,ckpt));
       G.Generator.load_state_dict(SD.get('g_ema'));
       G.Generator.to('cuda');G.Generator.eval();
       py.torch.set_grad_enabled(false)
       G.config = config;
   end
   
   function matimg = visualize_codes(G, style, truncation)
       if nargin == 2, truncation=0.7;end
       assert(size(style,2)==G.config.latent)
       if truncation < 1
       mean_latent = G.Generator.mean_latent(int32(4096));
       else
       mean_latent = py.None;
       end
       tic
       imgs = G.Generator(py.torch.tensor(py.numpy.array(style)).view(py.tuple(int32([1,size(style)]))).float().cuda(),...
           pyargs("truncation",truncation,"truncation_latent",mean_latent));
       toc
       imgs = imgs{1}; % discard the 2nd term in tuple
       matimg = imgs.detach.cpu().numpy().single; % note range in -1, 1
       matimg = permute(clip((matimg + 1) / 2, 0, 1),[3,4,2,1]);
   end
   
   end
end