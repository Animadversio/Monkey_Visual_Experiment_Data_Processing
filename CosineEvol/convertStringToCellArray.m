function updatedStructArray = convertStringToCellArray(inputStructArray)
    % Recursively goes through a struct array and converts any string variable to a cell array of chars.
    % Note this function is highly useful for saving a file in mat and load it using mat73 in python. 

    % Create a copy of the input struct array
    updatedStructArray = inputStructArray;

    % Check if the input is a struct array
    if length(inputStructArray) > 1
        % Iterate over each struct in the struct array
        for i = 1:length(inputStructArray)
            updatedStructArray(i) = convertStringToCellArray(inputStructArray(i));
        end
    else
        % Loop through each field of the struct
        fields = fieldnames(updatedStructArray);
        for i = 1:length(fields)
            field = fields{i};
            % Check if the field value is a struct or struct array
            if isstruct(updatedStructArray.(field))
                % Recursively call the function for sub-structures or struct arrays
                updatedStructArray.(field) = convertStringToCellArray(updatedStructArray.(field));
            elseif isstring(updatedStructArray.(field)) || ischar(updatedStructArray.(field))
                % Convert string or char array to cell array of chars
                updatedStructArray.(field) = cellstr(updatedStructArray.(field));
            elseif isdatetime(updatedStructArray.(field))
                updatedStructArray.(field) = cellstr(string(updatedStructArray.(field),'yyyyMMdd'));
                %string(datetime(updatedStructArray.(field),'Format','yyyyMMdd'));
            end
        end
    end
end