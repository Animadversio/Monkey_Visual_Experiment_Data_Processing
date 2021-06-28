function saveallform(figdir,fignm,h,sfxlist)
% Save a (current) figure with all suffices in a figdir.
% signature:
%   saveallform(figdir,fignm,h,sfxlist)
% 
if nargin <=3, h=gcf; end
if nargin <=4, sfxlist = ["fig","png","pdf"]; end
for sfx = sfxlist
if strcmp(sfx, "fig")
   savefig(h,fullfile(figdir,fignm+"."+sfx))
elseif strcmp(sfx, "pdf")
%    set(h,'Units','Inches');
%    pos = get(h,'Position');
%    set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
%    print(h,fullfile(figdir,fignm+"."+sfx),'-dpdf','-bestfit')
   exportgraphics(h,fullfile(figdir,fignm+"."+sfx),'ContentType','vector')
else
   saveas(h,fullfile(figdir,fignm+"."+sfx))
end
end
end