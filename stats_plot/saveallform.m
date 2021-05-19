function saveallform(figdir,fignm,h,sfxlist)
% Save a (current) figure with all suffices in a figdir. 
if nargin <=3, h=gcf; end
if nargin <=4, sfxlist = ["fig","pdf","png"]; end
for sfx = sfxlist
if strcmp(sfx, "fig")
   savefig(h,fullfile(figdir,fignm+"."+sfx))
else
   saveas(h,fullfile(figdir,fignm+"."+sfx))
end
end
end