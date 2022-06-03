classdef Coherence3d < matlab.System
    % COHERENCD3D: 3-D coherence function 
    %
    % Requirements: MATLAB R2020b
    %
    % Copyright (c) 2021, Ruiki KOBAYASHI and Shogo MURAMATSU
    %
    % All rights reserved.
    %
    % Contact address: Shogo MURAMATSU,
    %    Faculty of Engineering, Niigata University,
    %    8050 2-no-cho Ikarashi, Nishi-ku,
    %    Niigata, 950-2181, JAPAN
    %
    % http://msiplab.eng.niigata-u.ac.jp/
    %

    % Public, tunable properties
    properties (Nontunable)
        Scale       = 1.0;
        Sigma       = 8.0;
        Sigmaxy     = 8.0;
        Frequency   = 1/4;
        Sinc        = 1.5e-2;
        Kernel
        KernelZ
        KernelY
        KernelX
        %{
        %Kernel3d
        %Sinckernel
        SinckernelZ
        %Sigmakernel
        SigmazkernelZ
        %Sigmaxykernel
        SigmaxykernelXY
        SigmaxykernelZ
        %Freqkernel
        FreqkernelZ
        %Scalekernel
        ScalekernelZ
        %}
    end

    properties(DiscreteState)

    end

    % Pre-computed constants
    properties(Access = private)
    end

    methods
        function obj = Coherence3d(varargin)
            setProperties(obj,nargin,varargin{:});
            %
            len = 2*round(obj.Sigma/(obj.Sinc)^0.3/2)+1; % Set an odd number
            nz = -floor(len/2):floor(len/2);
            len = round(obj.Sigmaxy*12)+1; %2*round(obj.Sigmaxy*3)+1;
            nx = -floor(len/2):floor(len/2);
            gcx = exp(-nx.^2./(2*obj.Sigmaxy.^2));
            gcxy = gcx(:)*gcx(1,:);
            gcz = obj.Scale...
                *exp(-nz.^2./(2*obj.Sigma.^2)).*cos(2*pi*obj.Frequency*nz).*sinc(nz*obj.Sinc);
            gck = permute(gcz(:),[2 3 1]);
            obj.Kernel = pagemtimes(gck,gcxy);
            obj.KernelZ = gck;
            %obj.KernelXY = gcxy;
            obj.KernelX = gcx(:).';
            obj.KernelY = gcx(:);
            %{
            % Scale
            gcz = exp(-nz.^2./(2*obj.Sigma.^2)).*cos(2*pi*obj.Frequency*nz).*sinc(nz*obj.Sinc);
            gck = permute(gcz(:),[2 3 1]);
            %obj.Scalekernel = pagemtimes(gck,gcxy);
            obj.ScalekernelZ = gck;
            
            % Sinc
            gcz = obj.Scale*exp(-nz.^2./(2*obj.Sigma.^2)).*cos(2*pi*obj.Frequency*nz)...
                .*(cos(pi()*nz*obj.Sinc)/obj.Sinc-sinc(nz*obj.Sinc)/obj.Sinc);
            gck = permute(gcz(:),[2 3 1]);
            %obj.Sinckernel = pagemtimes(gck,gcxy);
            obj.SinckernelZ = gck;
            
            % Sigma
            gcz = obj.Scale*cos(2*pi*obj.Frequency*nz).*sinc(nz*obj.Sinc)...
                .*(nz.^2./obj.Sigma.^3).*exp(-nz.^2./(2*obj.Sigma.^2));
            gck = permute(gcz(:),[2 3 1]);
            %obj.Sigmakernel = pagemtimes(gck,gcxy);
            obj.SigmazkernelZ = gck;

            % Sigmaxy
            gcz = obj.Scale...
                *exp(-nz.^2./(2*obj.Sigma.^2)).*cos(2*pi*obj.Frequency*nz).*sinc(nz*obj.Sinc);
            gck = permute(gcz(:),[2 3 1]);
            nxx = nx(:).';
            nyy = nx(:);
            gcDiffxy = ((nxx.^2+nyy.^2)./(Sigmaxy.^3)).*exp(-(nxx.^2+nyy.^2)./(2*Sigmaxy.^2));
            %obj.Sigmaxykernel = pagemtimes(gck,gcDiffxy);
            obj.SigmaxykernelXY = gcDiffxy;
            obj.SigmaxykernelZ = gck;

            gcz = obj.Scale*exp(-nz.^2./(2*obj.Sigma.^2)).*sinc(nz*obj.Sinc)...
                .*(-2*pi*nz.*sin(2*pi*obj.Frequency*nz));
            gck = permute(gcz(:),[2 3 1]);
            %obj.Freqkernel = pagemtimes(gck,gcxy);
            obj.FreqkernelZ = gck;
            %}
        end
    end
    
    methods(Access = protected)
        function setupImpl(obj)
            % Perform one-time calculations, such as computing constants
        end

        function y = stepImpl(obj,u,dir)
            % Implement algorithm. Calculate y as a function of input u and
            % discrete states.
            if strcmp(dir,'Forward')
                v = imfilter(u,obj.KernelY,'conv','circ');
                w = imfilter(v,obj.KernelX,'conv','circ');
                y = imfilter(w,obj.KernelZ,'conv','circ');     
                %y_ = imfilter(u,obj.Kernel,'conv','circ');
                %norm(y(:)-y_(:))
            elseif strcmp(dir,'Adjoint')
                v = imfilter(u,obj.KernelY,'corr','circ');
                w = imfilter(v,obj.KernelX,'corr','circ');                
                y = imfilter(w,obj.KernelZ,'corr','circ');                                
                %y_ = imfilter(u,obj.Kernel,'corr','circ');
                %norm(y(:)-y_(:))
                %{
            elseif strcmp(dir,'Sinc')
                v = imfilter(u,obj.KernelY,'corr','circ');
                w = imfilter(v,obj.KernelX,'corr','circ');
                y = imfilter(w,obj.SinckernelZ,'corr','circ');
                %y_ = imfilter(u,obj.Sinckernel,'corr','circ');
                %norm(y(:)-y_(:))
            elseif strcmp(dir,'Freq')
                v = imfilter(u,obj.KernelY,'corr','circ');
                w = imfilter(v,obj.KernelX,'corr','circ');                
                y = imfilter(w,obj.FreqkernelZ,'corr','circ');                                
                %y_ = imfilter(u,obj.Freqkernel,'corr','circ');
                %norm(y(:)-y_(:))
            elseif strcmp(dir,'Sigma')
                v = imfilter(u,obj.KernelY,'corr','circ');
                w = imfilter(v,obj.KernelX,'corr','circ');
                y = imfilter(w,obj.SigmazkernelZ,'corr','circ');
                %y_ = imfilter(u,obj.Sigmakernel,'corr','circ');
                %norm(y(:)-y_(:))
            elseif strcmp(dir,'Sigmaxy')
                v = imfilter(u,obj.SigmaxykernelXY,'corr','circ');              
                y = imfilter(v,obj.SigmaxykernelZ,'corr','circ');
                %y_ = imfilter(u,obj.Sigmaxykernel,'corr','circ');
                %norm(y(:)-y_(:))
            elseif strcmp(dir,'Scale')
                v = imfilter(u,obj.KernelY,'corr','circ');
                w = imfilter(v,obj.KernelX,'corr','circ');
                y = imfilter(w,obj.ScalekernelZ,'corr','circ');
                %y_ = imfilter(u,obj.Scalekernel,'corr','circ');
                %norm(y(:)-y_(:))
                %}
            else
                error('Invalid mode')
            end
        end

        function resetImpl(obj)
            % Initialize / reset discrete-state properties
        end
    end
end
