# Analysis Pipeline Documentation

Neural code is intrinsically multifaceted. To parse out useful information and make sure to extract informative aspects as exhaustive as possible, we categorize the analysis applied on the data. 

## Analysis Pipeline for Manifold Study

* Pre-Processing 
  * `Project_Manifold_Beto_loadRaw`
  * Montage the images into a map! 
* General Characterization of Manifold Data
  * Compute the F statistics of ANOVA across images. 
  * Compute the t statistics between baseline and response period
  * Fit tuning function (*Kent* or other unimodal function) to the tuning landscape 
  * Un-parametrized characterization of the landscape. ??? Flood fill level set method
  * Correlation between channel
    * in baseline period
    * in response period, score period
    * noise correlation (variability across trial)
    * limited in some image subspace (PC23, PC4950, RND12) 
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
  * **Norm Trajectory** 
* 
* Focus on temporal dynamics
  * : (What's the time period in PSTH that contribute most to the increase of firing rate?)  