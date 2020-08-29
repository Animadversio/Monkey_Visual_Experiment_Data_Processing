classdef torchImDist
   % Usage: 
   % Visualizing a certian class 
   % 
   % Fix part of the code and Visualizing the other half
   %  G = G.select_space("class", noise_vec);
   %  G.visualize(0.06 * randn(5,128))
   properties
       D
       metric % a variable preset to specify which metric to use
   end
   methods
   function G = torchImDist(metric)
       if nargin == 0
           metric = "squeeze";
       end
       switch getenv('COMPUTERNAME')
           case 'DESKTOP-9DDE2RH' % Office 3 Binxu's 
            repodir = "D:\Github\PerceptualSimilarity";
           case 'DESKTOP-MENSD6S' % Binxu's home work station
            repodir = "E:\Github_Projects\PerceptualSimilarity";
           case 'PONCELAB-ML2A' % MLa machine 
            repodir = "C:\Users\Poncelab-ML2a\Documents\Python\PerceptualSimilarity";
           case 'PONCELAB-ML2B' % MLb machine 
            repodir = "C:\Users\Ponce lab\Documents\Python\PerceptualSimilarity";
           otherwise
            repodir = "C:\Users\Poncelab-ML2a\Documents\Python\PerceptualSimilarity";
        end
       % install the torch 1.3.x and the biggan package like below.
        py.importlib.import_module("sys");
        syspath = py.sys.path(); % add the official stylegan2 repo. 
        syspath.append(repodir);
        py.importlib.import_module("models");
        py.importlib.import_module("torch");
        py.importlib.import_module("numpy");
       
       G.metric = metric;
       switch metric
           case "SSIM"
               G.D = py.models.PerceptualLoss(pyargs("model", "SSIM"));
           case "L2"
               G.D = py.models.PerceptualLoss(pyargs("model", "L2"));
           case "alex"
               G.D = py.models.PerceptualLoss(pyargs("model", "net-lin", "net", "alex"));
           case "vgg"
               G.D = py.models.PerceptualLoss(pyargs("model", "net-lin", "net", "vgg"));
           case "squeeze"
               G.D = py.models.PerceptualLoss(pyargs("model", "net-lin", "net", "squeeze"));
           otherwise
               fprintf("Metric unrecognized, use squeeze-net linear weighted image metric.\n")
               G.D = py.models.PerceptualLoss(pyargs("model", "net-lin", "net", "squeeze"));
               G.metric = "squeeze";
       end
       G.D.requires_grad_(false);
       G.D.eval();%G.BGAN.to('cuda');
       py.torch.set_grad_enabled(false);
   end
   function G = select_metric(G, metric)
       G.metric = metric;
       switch metric
           case "SSIM"
               G.D = py.models.PerceptualLoss(pyargs("model", "SSIM"));
           case "L2"
               G.D = py.models.PerceptualLoss(pyargs("model", "L2"));
           case "alex"
               G.D = py.models.PerceptualLoss(pyargs("model", "net-lin", "net", "alex"));
           case "vgg"
               G.D = py.models.PerceptualLoss(pyargs("model", "net-lin", "net", "vgg"));
           case "squeeze"
               G.D = py.models.PerceptualLoss(pyargs("model", "net-lin", "net", "squeeze"));
           otherwise
               fprintf("Metric unrecognized, use squeeze-net linear weighted image metric.\n")
               G.D = py.models.PerceptualLoss(pyargs("model", "net-lin", "net", "squeeze"));
               G.metric = "squeeze";
       end
   end
   
   function dists = distance(G, im1, im2)
       if max(im1,[],'all')>1.2, im1 = single(im1) / 255.0; end
       if max(im2,[],'all')>1.2, im2 = single(im2) / 255.0; end
       % interface with generate integrated code, cmp to FC6GAN
       im1_tsr = py.torch.tensor(py.numpy.array(permute(im1,[4,3,1,2])));
       im2_tsr = py.torch.tensor(py.numpy.array(permute(im2,[4,3,1,2])));
       dists = G.D.forward(im1_tsr, im2_tsr).squeeze().cpu().detach().numpy().double;
   end
   
   function distMat = distmat(G, imgs, B)
       if nargin==2, B = 100;end
       if max(imgs,[],'all')>1.2, imgs = single(imgs) / 255.0; end
       % interface with generate integrated code, cmp to FC6GAN
       imgn = size(imgs,4);
       distMat = zeros(imgn,imgn);
       for i=1:imgn
           csr=1;
           dist_row = [];
           while csr <= imgn
           csr_end = min(imgn, csr+B-1);
           dists = G.distance(imgs(:,:,:,i), imgs(:,:,:,csr:csr_end));
           dist_row = cat(2, dist_row, dists);
           csr = csr_end+1;
           end
           distMat(i, :) = dist_row;
       end
   end
   
   end
end