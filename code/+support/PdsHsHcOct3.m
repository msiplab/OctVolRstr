classdef PdsHsHcOct3 < matlab.System
    % PDSHSHCOCT3: Primal-dual splitting for hierachical sparisty and 
    % range constraint
    %
    % Output:
    %
    %    Distribution of reflactive index
    %
    %
    % Requirements: MATLAB R2020b
    % Copyright (c) 2018-2021, Shogo MURAMATSU
    %
    % All rights reserved.
    %
    % Contact address: Shogo MURAMATSU,
    %                Faculty of Engineering, Niigata University,
    %                8050 2-no-cho Ikarashi, Nishi-ku,
    %                Niigata, 950-2181, JAPAN
    %
    % http://msiplab.eng.niigata-u.ac.jp/
    %
    
    % Public, tunable properties
    properties (Nontunable)
        Observation
        Lambda  = 0.01     % Regularization parameter
        Eta     = 0.01     % Regularization parameter
        Gamma1  = 0.01     % Step size
        Gamma2  = []
        Beta    = 0.0      % Lipschitz constant of gradient of fidelity term
        VRange  = [ 1.00 1.50 ]  % Range constraint
        PhiMode = 'Linear'       % Mode of map ƒ³
        MeasureProcess
        Dictionary
        GaussianDenoiser
        %
        SplitFactor = []
        PadSize     = [ 0 0 0 ]
        %
    end
    
    properties (GetAccess = public, SetAccess = private)
        Result
        LambdaCompensated
        EtaCompensated
    end
    
    properties(Nontunable, Access = private)
        dltFcn
        grdFcn
        parProc
    end
    
    properties(Nontunable, Logical)
        IsIntegrityTest = true
        IsSizeCompensation = false
        UseParallel = false
        UseGpu = false
    end
    
    properties(Nontunable,Logical, Hidden)
        Debug = false
    end
    
    properties(Access = private)
        y1
        y2
        xpre
        scls
    end
    
    properties (Hidden)
        PhiModeSet = ...
            matlab.system.StringSet(...
            {'Reflection','Linear','Signed-Quadratic','Identity'});
    end
    
    properties(DiscreteState)
        Iteration
    end
    
    methods
        function obj = PdsHsHcOct3(varargin)
            setProperties(obj,nargin,varargin{:})
            %
            import support.Sobel3d
            import support.RefractIdx2Reflect
            import support.CostEvaluator
            %
            obj.dltFcn = Sobel3d(...
                'KernelMode','Normal',...
                'UseGpu',obj.UseGpu);
            phi_ = RefractIdx2Reflect(...
                'PhiMode',obj.PhiMode,...
                'VRange',obj.VRange,...
                'UseGpu',obj.UseGpu);
            obj.grdFcn = CostEvaluator(...
                'Observation',obj.Observation,...
                'MeasureProcess',obj.MeasureProcess,...
                'RefIdx2Ref',phi_,...
                'OutputMode','Gradient',...
                'UseGpu',obj.UseGpu);
            %
            if isempty(obj.Gamma2)
                tauSqd     = obj.dltFcn.LambdaMax + 1;
                obj.Gamma2 = 1/(1.05*tauSqd)*(1/obj.Gamma1-obj.Beta/2);
            end
        end
    end
    
    methods(Access = protected)
        
        function s = saveObjectImpl(obj)
            s = saveObjectImpl@matlab.System(obj);
            s.Result = obj.Result;
            s.y1 = obj.y1;
            s.y2 = obj.y2;
            s.xpre = obj.xpre;
            s.scls = obj.scls;
            s.dltFcn = matlab.System.saveObject(obj.dltFcn);
            s.grdFcn = matlab.System.saveObject(obj.grdFcn);
            s.parProc = matlab.System.saveObject(obj.parProc);
            s.Dictionary{1} = matlab.System.saveObject(obj.Dictionary{1});
            s.Dictionary{2} = matlab.System.saveObject(obj.Dictionary{2});
            if isLocked(obj)
                s.Iteration = obj.Iteration;
            end
        end
        
        function loadObjectImpl(obj,s,wasLocked)
            if wasLocked
                obj.Iteration = s.Iteration;
            end
            obj.Dictionary{1} = matlab.System.loadObject(s.Dictionary{1});
            obj.Dictionary{2} = matlab.System.loadObject(s.Dictionary{2});
            obj.dltFcn = matlab.System.loadObject(s.dltFcn);
            obj.grdFcn = matlab.System.loadObject(s.grdFcn);
            obj.parProc = matlab.System.loadObject(s.parProc);
            obj.Result = s.Result;
            obj.xpre = s.xpre;
            obj.scls = s.scls;
            obj.y1 = s.y1;
            obj.y2 = s.y2;
            loadObjectImpl@matlab.System(obj,s,wasLocked);
        end
        
        
        function setupImpl(obj)
            % Perform one-time calculations, such as computing constants
            if obj.UseGpu
                obj.Observation = gpuArray(obj.Observation);
            end
            vObs = obj.Observation;
            %
            msrProc = obj.MeasureProcess;            
            fwdDic  = obj.Dictionary{1};
            adjDic  = obj.Dictionary{2};
            %
            if obj.IsSizeCompensation
                sizeM = numel(vObs); % Observation data size
                src   = msrProc.step(vObs,'Adjoint');
                sizeN = numel(src); % Size of refractive index distribution
                coefs = adjDic.step(src); % Size of coefficients
                sizeL = numel(coefs);
                obj.LambdaCompensated = obj.Lambda*(sizeM^2/sizeL);
                obj.EtaCompensated    = obj.Eta*(sizeM^2/sizeN);
            else
                obj.LambdaCompensated = obj.Lambda;
                obj.EtaCompensated    = obj.Eta;
            end
            %
            if obj.UseGpu
                obj.Gamma1 = gpuArray(obj.Gamma1);
                obj.Gamma2 = gpuArray(obj.Gamma2);
                obj.EtaCompensated = gpuArray(obj.EtaCompensated);
                obj.LambdaCompensated = gpuArray(obj.LambdaCompensated);
            end

            lambda_ = obj.LambdaCompensated;
            eta_    = obj.EtaCompensated;
            gamma1_ = obj.Gamma1;
            %gamma2_ = obj.Gamma2;
            %
            obj.GaussianDenoiser{1}.release();
            obj.GaussianDenoiser{1}.Sigma = sqrt(gamma1_*lambda_);
            obj.GaussianDenoiser{2}.release();
            %obj.GaussianDenoiser{2}.Sigma = sqrt(eta_/gamma2_);
            obj.GaussianDenoiser{2}.Sigma = sqrt(eta_); % Simplified by S. Muramatsu on Mar. 29 2022

            % Initialization
            %vmax = obj.VRange(2);
            vmin = obj.VRange(1);
            obj.y1 = zeros(1,'like',vObs);
            obj.y2 = zeros(1,'like',vObs);
            res0 = vmin*ones(size(vObs),'like',vObs); %zeros(size(vObs),'like',vObs); % Modified by S. Muramatsu  on 29 Mar. 2022
            obj.Result = res0; % zeros(1,'like',vObs); % Modified by S. Muramatsu on 29 Mar. 2022
            if isempty(obj.SplitFactor) % Normal process
                obj.parProc = [];
                %
                fwdDic.release();
                obj.Dictionary{1} = fwdDic.clone();
                adjDic.release();
                obj.Dictionary{2} = adjDic.clone();
                %
                [obj.xpre,obj.scls] = adjDic(res0); % Initial values of Coefs.
            else
                import saivdr.restoration.*
                gdn = obj.GaussianDenoiser{1};
                cm = CoefsManipulator(...
                    'Manipulation',...
                    @(t,cpre)  gdn.step(cpre-gamma1_*t));
                obj.parProc = OlsOlaProcess3d(...
                    'Synthesizer',fwdDic,...
                    'Analyzer',adjDic,...
                    'CoefsManipulator',cm,...
                    'SplitFactor',obj.SplitFactor,...
                    'PadSize',obj.PadSize,...
                    'UseParallel',obj.UseParallel,...
                    'UseGpu',obj.UseGpu,...
                    'IsIntegrityTest',obj.IsIntegrityTest,...
                    'Debug',obj.Debug);
                obj.xpre = obj.parProc.analyze(res0); % Initial values of Coefs.
                obj.parProc.InitialState = obj.xpre;
            end
            
        end
        
        function varargout = stepImpl(obj)
            % Implement algorithm. Calculate y as a function of input u and
            % discrete states.
            dltFcn_ = obj.dltFcn;
            grdFcn_ = obj.grdFcn;
            %
            vmin = obj.VRange(1);
            vmax = obj.VRange(2);
            gamma2_ = obj.Gamma2;
            gdnFcnH = obj.GaussianDenoiser{2};
            %
            y1_  = obj.y1;
            y2_  = obj.y2;
            up = obj.Result;
            % Primal step
            prx_ = grdFcn_.step(up) + (-dltFcn_.step(y1_)) + y2_;
            if isempty(obj.SplitFactor) % Normal process
                fwdDic  = obj.Dictionary{1};
                adjDic  = obj.Dictionary{2};
                gdnFcnG = obj.GaussianDenoiser{1};
                %
                scls_ = obj.scls;
                xpre_ = obj.xpre;
                gamma1_ = obj.Gamma1;
                %
                t_ = adjDic.step(prx_); % Analysis process
                xc = gdnFcnG.step(xpre_-gamma1_*t_); % Manipulation for Coefs.
                uc = fwdDic.step(xc,scls_); % Synthesis process
                %
                obj.xpre = xc;
            else % OLS/OLA Analysis-Synthesis process
                uc = obj.parProc.step(prx_);
            end
            % line 6
            q_ = 2*uc - up;
            % lines 7-8
            y1_ = y1_ + gamma2_*dltFcn_.step(q_);
            y2_ = y2_ + gamma2_*q_;
            % line 9
            %y1_ = y1_ - gamma2_*gdnFcnH.step( y1_/gamma2_ );
            y1_ = y1_ - gdnFcnH.step( y1_ );  % Simplified by S. Muramatsu on Mar. 29 2022
            % line 10
            pcy2 = y2_/gamma2_;
            pcy2((y2_/gamma2_)<vmin) = vmin;
            pcy2((y2_/gamma2_)>vmax) = vmax;
            y2_ = y2_ - gamma2_*pcy2;
            % 
            %r = (q_+up)/2; % FIXED by S. Muramatsu on 29 Mar. 2022

            % Output
            if nargout > 0
                varargout{1} = uc;
            end
            if nargout > 1
                rmse = norm(uc(:)-up(:),2)/norm(uc(:),2);
                varargout{2} = rmse;
            end
            
            % Update states
            obj.y1 = y1_;
            obj.y2 = y2_;
            obj.Result = uc;

            % line 11
            obj.Iteration = obj.Iteration + 1;
        end
        
        function resetImpl(obj)
            % Initialize / reset discrete-state properties
            obj.Iteration = 0;
        end
    end
  
end