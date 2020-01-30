%% NM-Co-F
% Cf. code from https://github.com/kimjingu/nonnegfac-matlab
par.init = [];
par.verbose = 0;
par.max_iter = 100;
par.reg_w = [0 0]; 
par.reg_h = [0 0]; 

    [m,n] = size(A); [m2,n2] = size(B); assert(n==n2, "A B should share the 2nd dimension. "); 
    par = params.Results;
    par.m = m;
    par.n = n;
    par.k = k;
    if isempty(par.init)
        W = rand(m,k); H = rand(k,n); V = rand(m2, k); 
    else
        W = par.init.W; H = par.init.H; 
    end

    initializer= hals_initializer;% str2func([par.method,'_initializer']);
    iterSolver = hals_iterSolver;% str2func([par.method,'_iterSolver']);
    iterLogger = hals_iterLogger;% str2func([par.method,'_iterLogger']);
    [W,H,par,val,ver] = feval(initializer,A,W,H,par);

%     if par.verbose && ~isempty(ver)
%         tTemp = cputime;
%         REC.HIS = saveHIS(1,ver,REC.HIS);
%         tPrev = tPrev+(cputime-tTemp);
%     end
% 
%     REC(1).par = par;
%     REC.start_time = datestr(now);
%     display(par);

    tStart = cputime;
    tTotal = 0;
    if par.track_grad
        initSC = getInitCriterion(par.stop_criterion,A,W,H,par);
    end
    SCconv = 0; SC_COUNT = 3;

    for iter=1:par.max_iter

        % Actual work of this iteration is executed here.
        [W,H,gradW,gradH,val] = feval(iterSolver,A,W,H,iter,par,val);

        if par.verbose          % Collect information for analysis/debugging
            elapsed = cputime-tPrev;
            tTotal = tTotal + elapsed;

%             clear('ver');
%             ver = prepareHIS(A,W,H,prev_W,prev_H,init,par,iter,elapsed,gradW,gradH);
% 
%             ver = feval(iterLogger,ver,par,val,W,H,prev_W,prev_H);
%             REC.HIS = saveHIS(iter+1,ver,REC.HIS);
% 
%             if par.track_prev, prev_W = W; prev_H = H; end
%             if par.verbose == 2, display(ver);, end
            tPrev = cputime;
        end

        if (iter > par.min_iter)
            if (par.verbose && (tTotal > par.max_time)) || (~par.verbose && ((cputime-tStart)>par.max_time))
                break;
            elseif par.track_grad
                SC = getStopCriterion(par.stop_criterion,A,W,H,par,gradW,gradH);
                if (SC/initSC <= par.tol)
                    SCconv = SCconv + 1;
                    if (SCconv >= SC_COUNT), break;, end
                else
                    SCconv = 0;
                end
            end
        end
    end
    [m,n]=size(A);
    [W,H]=normalize_by_W(W,H); % normalization scheme
    



% 'hals' : Hierarchical Alternating Least Squares Method
% Reference (See Algorithm 2):
%    Cichocki, A. and Phan, A.H.
%    Fast local algorithms for large scale nonnegative matrix and tensor factorizations.
%    IEICE Trans. Fundam. Electron. Commun. Comput. Sci. E92-A(3), 708–721 (2009)

function [W,H,par,val,ver] = hals_initializer(A,W,H,par)
    [W,H]=normalize_by_W(W,H);

    val = struct([]);
    ver = struct([]);
end

function [W,H,gradW,gradH,val] = hals_iterSolver(A,W,H,iter,par,val)
    epsilon = 1e-16;

    WtA = W'*A;
    WtW = W'*W;
    WtW_reg = applyReg(WtW,par,par.reg_h);
    for i = 1:par.k
        H(i,:) = max(H(i,:) + WtA(i,:) - WtW_reg(i,:) * H,epsilon);
    end

    AHt = A*H';
    HHt_reg = applyReg(H*H',par,par.reg_w);
    for i = 1:par.k
        W(:,i) = max(W(:,i) * HHt_reg(i,i) + AHt(:,i) - W * HHt_reg(:,i),epsilon);
        if sum(W(:,i))>0
            W(:,i) = W(:,i)/norm(W(:,i));
        end
    end

    if par.track_grad
        gradW = W*HHt_reg - AHt;
        gradH = getGradientOne(W'*W,W'*A,H,par.reg_h,par);
    else
        gradH = 0; gradW = 0;
    end
end

function [ver] = hals_iterLogger(ver,par,val,W,H,prev_W,prev_H)
end
%-------------------------------------------------------------------------------
function AtA = applyReg(AtA,par,reg)
    % Frobenius norm regularization
    if reg(1) > 0
        AtA = AtA + 2 * reg(1) * eye(par.k);
    end
    % L1-norm regularization
    if reg(2) > 0
        AtA = AtA + 2 * reg(2) * ones(par.k,par.k);
    end
end

function [grad] = modifyGradient(grad,X,reg,par)
    if reg(1) > 0
        grad = grad + 2 * reg(1) * X;
    end
    if reg(2) > 0
        grad = grad + 2 * reg(2) * ones(par.k,par.k) * X;
    end
end

function [grad] = getGradientOne(AtA,AtB,X,reg,par)
    grad = AtA*X - AtB;
    grad = modifyGradient(grad,X,reg,par);
end

function [gradW,gradH] = getGradient(A,W,H,par)
    HHt = H*H';
    HHt_reg = applyReg(HHt,par,par.reg_w);

    WtW = W'*W;
    WtW_reg = applyReg(WtW,par,par.reg_h);

    gradW = W*HHt_reg - A*H';
    gradH = WtW_reg*H - W'*A;
end

%-------------------------------------------------------------------------------
function pGradF = projGradient(F,gradF)
    pGradF = gradF(gradF<0|F>0);
end

%-------------------------------------------------------------------------------
function [W,H,weights] = normalize_by_W(W,H)
    norm2=sqrt(sum(W.^2,1));
    toNormalize = norm2>0;

    if any(toNormalize)
        W(:,toNormalize) = W(:,toNormalize)./repmat(norm2(toNormalize),size(W,1),1);
        H(toNormalize,:) = H(toNormalize,:).*repmat(norm2(toNormalize)',1,size(H,2));
    end

    weights = ones(size(norm2));
    weights(toNormalize) = norm2(toNormalize);
end
