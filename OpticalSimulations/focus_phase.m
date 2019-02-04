function phase = focus_phase(obj,f,X,Y,Z,varargin)
            % FOCUS Focal fields
            % 
            % E = FOCUS(B,F,X,Y,Z) calculates the focal fields of B 
            %   at positions X, Y, Z for a lens with focal lenght F
            % 
            % E = FOCUS(B,F,X,Y,Z,'PropertyName',PropertyValue) sets the property
            %   PropertyName to PropertyValue. The properties listed below
            %   can be used:
            %       er  -   relative electric permittivity [default: B.er]
            %       mr  -   relative magnetic permeability [default: B.mr]
            %
            % See also Beam.

            % Relative dielectric permittivity
            er1 = obj.er; % before objective
            er2 = er1; % after objective
            for n = 1:2:length(varargin)
                if strcmpi(varargin{n},'er')
                    er2 = varargin{n+1};
                    Check.isnumeric('er must be a number',er2)
                end
            end

            % Relative magnetic permeability
            mr1 = obj.mr; % before objective
            mr2 = mr1; % after objective
            for n = 1:2:length(varargin)
                if strcmpi(varargin{n},'mr')
                    mr2 = varargin{n+1};
                    Check.isnumeric('mr must be a number',mr2)
                end
            end
            
            % Numerically more stable version
            phi = obj.phi;
            
            theta = asin(obj.r/f);
            
            k = 2*pi*sqrt(er2*mr2)/obj.lambda0;
            X = [X,X-1e-12,X+1e-12,X,X,X,X];
            Y = [Y,Y,Y,Y-1e-12,Y+1e-12,Y,Y];
            Z = [Z,Z,Z,Z,Z,Z-1e-12,Z+1e-12];
            rho = sqrt(X.^2+Y.^2);
            varphi = atan2(Y,X);
            phase=zeros(size(theta,1),size(theta,2),7);
            for n=1:7
                phase(:,:,n) = exp(1i*(k*cos(theta)*Z(n)+k*sin(theta)*rho(n).*cos(varphi(n)-phi)));
            end
end
