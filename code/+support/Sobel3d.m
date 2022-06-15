classdef Sobel3d < matlab.System % codegen
    % SOBEL3D: 3-D Sobel filter class
    %
    % Requirements: MATLAB R2020b
    %
    % Copyright (c) 2021-2022, Ruiki KOBAYASHI and Shogo MURAMATSU
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
        Kernel
        KernelX
        KernelY
        KernelZ
        KernelMode = 'Normal'
        LambdaMax
    end
    
    properties (Nontunable, Logical)
        UseGpu = false
    end
    
    properties (Hidden, Transient)
        KernelModeSet = ...
            matlab.system.StringSet({'Normal','Absolute'});
        
    end
    
    properties(DiscreteState)

    end

    % Pre-computed constants
    properties(Access = private)

    end

    
    methods
        function obj = Sobel3d(varargin)
            setProperties(obj,nargin,varargin{:});
            %
            obj.KernelX = [ 1 2 1 ].'/4;
            obj.KernelY = [ 1 2 1 ]/4;
            obj.KernelZ  = permute([ 1 0 -1 ].',[ 2 3 1 ])/2;
            KernelXY_ = kron([ 1 2 1 ].', [ 1 2 1 ])/16;
            obj.Kernel = convn(KernelXY_,obj.KernelZ);
            if strcmp(obj.KernelMode,'Absolute')
                obj.KernelZ = abs(obj.KernelZ);
                obj.Kernel = abs(obj.Kernel);
            end            
            %
            obj.LambdaMax = sum(abs(obj.Kernel(:)))^2; 
        end
            
    end
    
    methods(Access = protected)
        function setupImpl(obj)
            % Perform one-time calculations, such as computing constants 
        end

        function y = stepImpl(obj,u)
            % Implement algorithm. Calculate y as a function of input u and
            % discrete states.
            if obj.UseGpu
                u = gpuArray(u);
            end
            v = imfilter(u,obj.KernelX,'conv','circ');
            w = imfilter(v,obj.KernelY,'conv','circ');
            y = imfilter(w,obj.KernelZ,'conv','circ');
            %y_ = imfilter(u,obj.Kernel,'conv','circ');
            %norm(y(:)-y_(:))
            if obj.UseGpu
                y = gather(y);
            end
        end

        function resetImpl(obj)
            % Initialize / reset discrete-state properties
        end
    end
end

