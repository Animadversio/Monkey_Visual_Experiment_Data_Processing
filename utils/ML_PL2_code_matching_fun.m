load("ML_PL2_code_matching.mat")
%%
strfind(words',allCodes(1:3610)');
%%
strfind(words',allCodes(3611:4252)');
%%
strfind(words',allCodes(14695:end)');
%%
csr = 1;
csr_end = 1;
targ_csr = 1;
targ_csr_end = 1;
subseq = {};
targsubseq = {};
baseLen = 50; % length of the default probe sequence
while true
    if targ_csr_end > numel(words) || csr_end > numel(allCodes)
        subseq{end+1} = [csr,csr_end-1];
        targsubseq{end+1} = [targ_csr,targ_csr_end-1];
        break
    end
    if words(targ_csr_end) == allCodes(csr_end)
        csr_end = csr_end + 1;
        targ_csr_end = targ_csr_end + 1;
    else
        subseq{end+1} = [csr,csr_end-1];
        targsubseq{end+1} = [targ_csr,targ_csr_end-1];
        disp(subseq)
        csr = csr_end;
        targ_csr = targ_csr_end;
        curlen = baseLen;
        while true
            relidx = strfind(words(targ_csr:end)',allCodes(csr:csr+curlen)');
            if ~isempty(relidx)
                targ_csr = targ_csr + relidx(1) - 1;
                targ_csr_end = targ_csr + curlen;
                csr_end = csr+curlen;
                break
            else
                curlen = max(10,round(curlen / 2)); % shorten the probe sequence length
            end
        end
    end
end
%% 
assert(sum(cellfun(@(A)A(2)-A(1)+1,subseq)) == numel(allCodes))
assert(sum(cellfun(@(A)A(2)-A(1)+1,targsubseq)) == numel(allCodes))
%%
missing_words = words(cellfun(@(A)A(1)-1, targsubseq(2:end)));

