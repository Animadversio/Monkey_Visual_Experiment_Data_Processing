function setMatlabTitle(name)
if nargin == 0
    com.mathworks.mlservices.MatlabDesktopServices.getDesktop.getMainFrame.getTitle()
else
    com.mathworks.mlservices.MatlabDesktopServices.getDesktop.getMainFrame.setTitle(name)
end
end