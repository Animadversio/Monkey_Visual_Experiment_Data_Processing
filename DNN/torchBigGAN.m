% pe = pyenv('Version','C:\ProgramData\Anaconda3\envs\tf-torch\python.exe'); %
% pe = pyenv('Version','C:\Users\ponce\.conda\envs\caffe36\python.exe'); %
% lab machine
% Set up the python executable before usage
classdef torchBigGAN
   % Usage 
   % BGAN = torchBigGAN("biggan-deep-512")
   % matimgs = BGAN.visualize_class(0.6*randn(5,128),729);figure;montage(matimgs)
   properties
       BGAN
       Embeddings
       Generator
   end
   methods
   function G = torchBigGAN(modelname)
       if nargin == 0
           modelname = "biggan-deep-256";
       end
       % download these weights and definitions and put them in savedir.
        % resolved_config_file = "https://s3.amazonaws.com/models.huggingface.co/biggan/biggan-deep-256-config.json";
        % resolved_model_file = "https://s3.amazonaws.com/models.huggingface.co/biggan/biggan-deep-256-pytorch_model.bin";
       savedir = "C:\Users\binxu\.pytorch_pretrained_biggan";
       savedir = "C:\Users\ponce\.pytorch_pretrained_biggan";
       % install the torch 1.3.x and the biggan package like below.
       py.importlib.import_module('torch');
       py.importlib.import_module('pytorch_pretrained_biggan');
%        import py.pytorch_pretrained_biggan.BigGAN
%        import py.pytorch_pretrained_biggan.BigGANConfig
       cfg = py.pytorch_pretrained_biggan.BigGANConfig();
       cfg = cfg.from_json_file(fullfile(savedir,compose("%s-config.json",modelname)));
       G.BGAN = py.pytorch_pretrained_biggan.BigGAN(cfg);
       G.BGAN.load_state_dict(py.torch.load(fullfile(savedir,compose("%s-pytorch_model.bin",modelname))));
       G.BGAN.to('cuda');G.BGAN.eval()
       py.torch.set_grad_enabled(false)
       tmp = py.list(G.BGAN.named_children);
       G.Embeddings = tmp{1}{2};
       G.Generator = tmp{2}{2};
   end
   
   function matimg = visualize_codes(G, noise, onehot, truncation)
       if nargin == 3, truncation=0.7;end
       tic
       imgs = G.BGAN(py.torch.tensor(py.numpy.array(noise)).view(int32(-1),int32(128)).float().cuda(),...
            py.torch.tensor(py.numpy.array(onehot)).view(int32(-1),int32(1000)).float().cuda(),truncation);
       toc
       matimg = imgs.detach.cpu().numpy().single;
       matimg = permute((matimg + 1) / 2.0,[3,4,2,1]);
   end
   
   function matimg = visualize_class(G, noise, classn, truncation)
       if nargin == 3, truncation=0.7;end
       onehot = zeros(size(noise,1),1000); onehot(:,classn)=1;
       tic
       imgs = G.BGAN(py.torch.tensor(py.numpy.array(noise)).view(int32(-1),int32(128)).float().cuda(),...
            py.torch.tensor(py.numpy.array(onehot)).view(int32(-1),int32(1000)).float().cuda(),truncation);
       toc
       matimg = imgs.detach.cpu().numpy().single;
       matimg = permute((matimg + 1) / 2.0,[3,4,2,1]);
   end
   
   function matimgs = visualize_latent(G, latent, truncation)
       if nargin == 2, truncation=0.7;end
       batchsize = 10;samplen = size(latent,1);csr = 1;
       tic
       matimgs = [];
       while csr <= samplen
       cnd = min(samplen,csr+batchsize);
       imgs = G.Generator(py.torch.tensor(py.numpy.array(latent(csr:cnd,:))).view(int32(-1),int32(256)).float().cuda(), truncation);
       matimg = imgs.detach.cpu().numpy().single;
       matimg = permute((matimg + 1) / 2.0,[3,4,2,1]);
       matimgs = cat(4, matimgs, matimg);
       csr = cnd + 1;
       end
       toc
       
   end
   
   function EmbedVects_mat = get_embedding(G)
       EmbedVects = py.list(G.Embeddings.parameters());
       EmbedVects_mat = EmbedVects{1}.data.cpu().numpy().double;
   end
   end
end
%%

% model = BigGAN.from_pretrained('biggan-deep-256')
% model.to('cuda')
% 
% def BigGAN_render(class_vector, noise_vector, truncation):
%     if class_vector.shape[0] == 1:
%         class_vector = np.tile(class_vector, [noise_vector.shape[0], 1])
%     if noise_vector.shape[0] == 1:
%         noise_vector = np.tile(noise_vector, [class_vector.shape[0], 1])
%     class_vector = torch.from_numpy(class_vector.astype(np.float32)).to('cuda')
%     noise_vector = torch.from_numpy(noise_vector.astype(np.float32)).to('cuda')
%     with torch.no_grad():
%         output = model(noise_vector, class_vector, truncation)
%     imgs = convert_to_images(output.cpu())
%     return imgs