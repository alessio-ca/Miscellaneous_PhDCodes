classdef Check
    % Check : Input validation
    %
    % Check methods:
	%   samesize        -   (Static) Validate same size
    %   isnumeric       -   (Static) Validate number
    %   isreal          -   (Static) Validate real number
    %   isinteger       -   (Static) Validate integer number
    %   isa             -   (Static) Validate object class
    %   samedim         -   (Static) Validate same number

	%   Author: Giovanni Volpe
    %   Modification: Alessio Caciagli
    %   Revision: 1.0.0  
    %   Date: 2015/01/01

    methods (Static)
        function samesize(msg,x,y,varargin)
            % SAMESIZE Validate same size
            %
            % SAMESIZE(MSG,X,Y) returns the error message MSG 
            %   if the sizes of X and Y are not equal.
            %
            % More than two variables can be checked at once, e.g., 
            %   SAMESIZE(MSG,X,Y,Z,W,T)
            %
            % See also Check.
            
            if size(x)~=size(y)
                error(msg);
            end
            for n = 1:1:length(varargin)
                if size(x)~=size(varargin{n})
                    error(msg);
                end
            end
            
        end
        function isnumeric(msg,x,varargin)
            % ISNUMERIC Validate number
            %
            % ISNUMERIC(MSG,X) returns the error message MSG if X is not a number.
            %
            % ISNUMERIC(MSG,X,'~=',Y) returns the error message MSG if X~=Y is not verified.
            %
            % ISNUMERIC(MSG,X,'<',Y) returns the error message MSG if X<Y is not verified.
            %
            % ISNUMERIC(MSG,X,'<=',Y) returns the error message MSG if X<=Y is not verified.
            %
            % ISNUMERIC(MSG,X,'>',Y) returns the error message MSG if X>Y is not verified.
            %
            % ISNUMERIC(MSG,X,'>=',Y) returns the error message MSG if X>=Y is not verified.
            %
            % Multiple couples of relation and value are also allowed, e.g., 
            %   ISNUMERIC(MSG,X,'>',0,'<',10)
            %
            % See also Check.

            if ~isnumeric(x)
                error(msg);
            end
            
            for n = 1:2:length(varargin)
                if strcmpi(varargin{n},'~=')
                    if ~(x~=varargin{n+1})
                        error(msg);
                    end
                end
                if strcmpi(varargin{n},'>')
                    if ~(x>varargin{n+1})
                        error(msg);
                    end
                end
                if strcmpi(varargin{n},'>=')
                    if ~(x>=varargin{n+1})
                        error(msg);
                    end
                end
                if strcmpi(varargin{n},'<')
                    if ~(x<varargin{n+1})
                        error(msg);
                    end
                end
                if strcmpi(varargin{n},'<=')
                    if ~(x<=varargin{n+1})
                        error(msg);
                    end
                end
            end
        end
        function isreal(msg,x,varargin)
            % ISREAL Validate real number
            %
            % ISREAL(MSG,X) returns the error message MSG if X is not a real number.
            %
            % ISREAL(MSG,X,'~=',Y) returns the error message MSG if X~=Y is not verified.
            %
            % ISREAL(MSG,X,'<',Y) returns the error message MSG if X<Y is not verified.
            %
            % ISREAL(MSG,X,'<=',Y) returns the error message MSG if X<=Y is not verified.
            %
            % ISREAL(MSG,X,'>',Y) returns the error message MSG if X>Y is not verified.
            %
            % ISREAL(MSG,X,'>=',Y) returns the error message MSG if X>=Y is not verified.
            %
            % Multiple couples of relation and value are also allowed, e.g., 
            %   ISREAL(MSG,X,'>',0,'<',10)
            %
            % See also Check.
            
            Check.isnumeric(msg,x,varargin{:});
            
            if ~isreal(x)
                error(msg);
            end
        end
        function isinteger(msg,x,varargin)
            % ISINTEGER Validate integer number
            %
            % ISINTEGER(MSG,X) returns the error message MSG if X is not an integer number.
            %
            % ISINTEGER(MSG,X,'~=',Y) returns the error message MSG if X~=Y is not verified.
            %
            % ISINTEGER(MSG,X,'<',Y) returns the error message MSG if X<Y is not verified.
            %
            % ISINTEGER(MSG,X,'<=',Y) returns the error message MSG if X<=Y is not verified.
            %
            % ISINTEGER(MSG,X,'>',Y) returns the error message MSG if X>Y is not verified.
            %
            % ISINTEGER(MSG,X,'>=',Y) returns the error message MSG if X>=Y is not verified.
            %
            % Multiple couples of relation and value are also allowed, e.g., 
            %   ISINTEGER(MSG,X,'>',0,'<',10)
            %
            % See also Check.
            
            Check.isnumeric(msg,x,varargin{:});
            
            if x~=round(x)
                error(msg);
            end
        end
        function isa(msg,obj,class,varargin)
            % ISA Validate same size
            %
            % ISA(MSG,OBJ,CLASS) returns the error message MSG 
            %   if the class of the object OBJ is not CLASS.
            %
            % More than one class can be checked at once, e.g., 
            %   ISA(MSG,OBJ,CLASS1,CLASS2,CLASS3)
            %
            % See also Check.

            check = false;
            if isa(obj,class)
                check = true;
            end
            for n = 1:1:length(varargin)
                if isa(obj,varargin{n})
                    check = true;
                end
            end
            if ~check
                error(msg);
            end
        end
        function samedim(msg,x,y)
            % SAMEDIM   Validate same dimension
            %
            % SAMEDIM(MSG,X,Y) returns the error message MSG
            %   if the dimension of Y is not X.
            %
            % See also Check.
            
            switch x
                case 3
                    if ~(length(size(y))==3)
                        error(msg)
                    end
                case 2
                    if ~(length(size(y))==2)
                        error(msg)
                    end
                    if ~all(size(y) - [1 1])
                        error(msg')
                    end
                case 1
                    if ~(length(size(y))==2)
                        error(msg)
                    end
                    if all(size(y) - [1 1])
                        %2D matrix.
                        error(msg)
                    end
                    if ~any(size(y) - [1 1])
                        %0D matrix.
                        error(msg)
                    end
                otherwise
                    error(msg)
            end
        end
        function samemesh(msg,X,Y)
            % SAMEMESH  Validate that quantity X (ComplexVector) has same mesh as quantity Y (ComplexVector)
            %
            % SAMEMESH (MSG,X,Y) returns the error message MSG
            %   if the quantity X does not have the same mesh as Y.
            %
            %
            % See also Check.
            
            if ~all(all([X.X,X.Y,X.Z]==[Y.X,Y.Y,Y.Z]))
                error(msg)
            end
        end
        function outofbound(msg,X,Y,varargin)
            % OUTOFBOUND  Validate that quantity Y (ComplexVector) is within
            % bound of X (ComplexVector)
            %
            % OUTOFBOUND(MSG,X,Y,'dr',DR) sets the tolerance to DR [default = 1e-10 m].
            %
            %
            %
            % See also Check.
            
            % increment [m]
            dr = 1e-10;
            for n = 1:2:length(varargin)
                if strcmpi(varargin{n},'dr')
                    dr = varargin{n+1};
                    Check.isnumeric('dr must be a positive real number',dr,'>',0)
                end
            end
            
            if(any([max(max((X.X)))-dr max(max((X.Y)))-dr max(max((X.Z)))-dr]>[max(max((Y.X))) max(max((Y.Y))) max(max((Y.Z)))]))
                error(msg)
            end
            if(any([min(min((X.X)))+dr min(min((X.Y)))+dr min(min((X.Z)))+dr]<[min(min((Y.X))) min(min((Y.Y))) min(min((Y.Z)))]))
                error(msg)
            end
        end
    end
end