function area = area_map(chan,array_layout)
% defaultly use the old layout pattern 
if nargin == 1, array_layout = "Alfa"; end
area = strings(size(chan));
if strcmp(array_layout,"Beto") || strcmp(array_layout,"Alfa")
    area(chan<=48 & chan>=33) = "V1";
	area(chan>48) = "V4";
	area(chan<33) = "IT";
	% if chan<=48 & chan>=33
	% area = "V1";
	% elseif chan>48
	% area = "V4";
	% elseif chan<33
	% area = "IT";
	% end
elseif strcmp(array_layout,"Beto_new")
	area(chan<=32 & chan>=17) = "V1";
	area(chan<17) = "V4";
	area(chan>=33) = "IT";
	% if chan<=32 & chan>=17
	% area = "V1";
	% elseif chan<17
	% area = "V4";
	% elseif chan>=33
	% area = "IT";
	% end
else 
	error("array layout not recognized")
end
end