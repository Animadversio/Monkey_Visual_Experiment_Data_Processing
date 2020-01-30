classdef FC6Generator
   properties
      BGRMean = reshape([104.0, 117.0, 123.0],[1,1,3,1]);
      DeconvNet
      LinNet
   end
   methods
      function G = FC6Generator(loadpath)
         %D:\Github\Monkey_Visual_Experiment_Data_Processing
         data = load(loadpath,'DeconvNet','LinNet','BGRMean');
         G.DeconvNet = data.DeconvNet;
         G.LinNet = data.LinNet;
      end
      function imgs = visualize(G, codes)
        if size(codes,2)==4096 
            codes = codes';
        end
        assert(size(codes,1)==4096, "Code shape error")% first dimension should have size 4096, second the batch
        Z = dlarray(reshape(codes,1,1,4096,size(codes,2)), 'SSCB');
        hiddenout = G.LinNet.predict(Z);
        hiddenout_r = dlarray(hiddenout.reshape(4,4,256,[]).permute([2,1,3,4]),"SSCB"); % Note to flip your hiddenoutput
        out = G.DeconvNet.predict(hiddenout_r);
        imgs = extractdata(out(:,:,:,:));
        imgs = uint8(clip(imgs + G.BGRMean, 0, 255));
        imgs = imgs(:,:,[3,2,1],:);
      end
   end
end