# Analysis Pipeline Documentation

Neural code is intrinsically multifaceted. To parse out useful information and make sure to extract informative aspects as exhaustive as possible, we categorize the analysis applied on the data. 

## Data Management

* All the manifold experiments have a evolution  (`generate_parallel`) component and a exploration (`selectivity_basic` experiment) component. 
* We store all evolution experiments in `Manifold_Evolv_Exps.mat`, and all the manifold experiment in `Manifold_Exps.mat`, arranged in time order, numbered by `Expi`
* `Project_Manifold_Beto_loadRaw` codes the name and path of source file (`pl2`,`bhv2`,images path) for each experiment. 
* Use the file `Set_Exp_Specs.m` to load the specific information of each experiment. 

## Analysis Pipeline for Manifold Study

* Pre-Processing 
  * `Project_Manifold_Beto_loadRaw` extract the spikes and lfps and segment them into trial structure. 
  * Montage the images into a map! 
* General Characterization of Manifold Data
  * Compute the F statistics of ANOVA across images. 
  * Compute the t statistics between baseline and response period
  * Fit tuning function (*Kent* or other unimodal function) to the tuning landscape and extract statistics 
  * Plot the tuning map (hot spot) across channels: see if they tile the image space.
  * Un-parametrized characterization of the landscape. 
    * Flood fill / level set method (size of hot spot)
  * Correlation between channel
    * in baseline period
    * in response period, score period
    * noise correlation (variability across trial)
    * limited in some image subspace (PC23, PC4950, RND12) 
    * Network analysis
  * Montage the image place maps together. 
* Focus on **Fluctuation** 
  * **Baseline fluctuation** and correlation  
* Focus on Temporal Dynamics
  * Is there any tuning dynamics?
* Focus on **Sequence**?
  * Is there sequential effect affecting the response to a single image. 
* Focus on **LFP and frequency** aspects
  * Is LFP modulating the firing? 



## Analysis Pipeline for Evolution Trajectory

* General Characterization 
  * Successfulness of an evolution 
  * **Norm Trajectory**, Slope and intercept of *Norm^2 ~ geni* curve. 
* 
* Focus on temporal dynamics
  * : (What's the time period in PSTH that contribute most to the increase of firing rate?)  