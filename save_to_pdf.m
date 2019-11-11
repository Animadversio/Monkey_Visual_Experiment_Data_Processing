% function to take an image handle and save to the correct size pdf
function save_to_pdf(fig_h, path)
set(fig_h,'Units','Inches');
pos = get(fig_h,'Position');
set(fig_h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(fig_h, path,'-dpdf','-r0')
end