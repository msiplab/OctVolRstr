classdef Pds1dtv < matlab.System
    % PDS1DTV: Primal-dual splitting for 1-D TVD
    %
    % Output:
    %
    %    Distribution of reflactive index
    %
    % Requirements: MATLAB R2020b
    %
    % Reference:
    %
    % - M. Shamouilian and I. Selesnick, "Total Variation Denoising for Optical
    %   Coherence Tomography," 2019 IEEE Signal Processing in Medicine and Biology
    %   Symposium (SPMB), 2019, pp. 1-5, doi: 10.1109/SPMB47826.2019.9037832.
    %
    % Copyright (c) 2022, Shogo MURAMATSU
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
        Eta     = 0.01     % Regularization parameter
        Gamma1  = 0.01     % Step size
        Gamma2  = []
        Beta    = 0.0      % Lipschitz constant of gradient of fidelity term
        %
        MeasureProcess
    end

    properties (GetAccess = public, SetAccess = private)
        Result
        EtaCompensated
    end

    properties(Nontunable, Access = private)
    end

    properties(Nontunable, Logical)
        IsSizeCompensation = false
        UseGpu = false
    end

    properties(Nontunable,Logical, Hidden)
        Debug = false
    end

    properties(Access = private)
        y
    end

    properties (Hidden)
    end

    properties(DiscreteState)
        Iteration
    end

    methods
        function obj = Pds1dtv(varargin)
            setProperties(obj,nargin,varargin{:})
            %
            if isempty(obj.Gamma2)
                tauSqd = 4; % max(abs(fftn([1 -1]).^2,[],'all'); % (ƒÐmax(G))^2
                obj.Gamma2 = 1/(1.05*tauSqd)*(1/obj.Gamma1-obj.Beta/2);
            end
        end
    end

    methods(Access = protected)

        function s = saveObjectImpl(obj)
            s = saveObjectImpl@matlab.System(obj);
            s.Result = obj.Result;
            s.yp = obj.y;
            s.xpre = obj.xpre;
            s.scls = obj.scls;
            if isLocked(obj)
                s.Iteration = obj.Iteration;
            end
        end

        function loadObjectImpl(obj,s,wasLocked)
            if wasLocked
                obj.Iteration = s.Iteration;
            end
            obj.Result = s.Result;
            obj.xpre = s.xpre;
            obj.scls = s.scls;
            obj.y = s.yp;
            loadObjectImpl@matlab.System(obj,s,wasLocked);
        end

        function setupImpl(obj)
            % Perform one-time calculations, such as computing constants
            if obj.UseGpu
                obj.Observation = gpuArray(obj.Observation);
            end
            
            vObs = obj.Observation;
            msrProc = obj.MeasureProcess;
            %
            if obj.IsSizeCompensation
                sizeM = numel(vObs); % Observation data size
                src   = msrProc.step(vObs,'Adjoint');
                sizeN = numel(src); % Size of refractive index distribution
                obj.EtaCompensated    = obj.Eta*(sizeM^2/sizeN);
            else
                obj.EtaCompensated    = obj.Eta;
            end

            %Initialize
            fdf = shiftdim([1 -1],-1);
            if obj.UseGpu            
                fdf = gpuArray(fdf);
                obj.Gamma1 = gpuArray(obj.Gamma1);
                obj.Gamma2 = gpuArray(obj.Gamma2);
                obj.EtaCompensated = gpuArray(obj.EtaCompensated);
            end            
            diffproc_ = @(x) imfilter(x,fdf,'conv','circ');    
            obj.y = diffproc_(vObs);
            obj.Result = vObs;

        end

        function varargout = stepImpl(obj)
            % Implement algorithm. Calculate y as a function of input u and
            % discrete states.
            %
            rmin = -1;
            rmax = 1;
            eta_ = obj.EtaCompensated;
            gamma1_ = obj.Gamma1;
            gamma2_ = obj.Gamma2;
            %
            msrproc_ = obj.MeasureProcess;
            fdf = shiftdim([1 -1],-1);
            if obj.UseGpu            
                fdf = gpuArray(fdf);
            end
            diffproc_ = @(x) imfilter(x,fdf,'conv','circ');
            diffadjp_ = @(x) circshift(imfilter(x,fdf,'corr','circ'),[0 0 1]);             
            softthresh_ = @(x,t) sign(x).*max(abs(x)-t,0);
            distproj_ = @(x,xmin,xmax) min(max(x,xmin),xmax);
            % 
            yp_ = obj.y;
            rp_ = obj.Result;
            v_  = obj.Observation;
            % Primal step
            rg_ = msrproc_.step(msrproc_.step(rp_,'Forward')-v_,'Adjoint');
            rc_ = distproj_(rp_-gamma1_*(rg_ + diffadjp_(yp_)),rmin,rmax);
            % Dual step
            yt_ = yp_ + gamma2_*diffproc_(2*rc_ - rp_);
            yc_ = yt_ - softthresh_(yt_,eta_);

            % Output
            if nargout > 0
                varargout{1} = rc_;
            end
            if nargout > 1
                rmse = norm(rc_(:)-rp_(:),2)/norm(rc_(:),2);
                varargout{2} = rmse;
            end

            % Update states
            obj.y = yc_;
            obj.Result = rc_;
            obj.Iteration = obj.Iteration + 1;
        end

        function resetImpl(obj)
            % Initialize / reset discrete-state properties
            obj.Iteration = 0;
        end
    end
  
end