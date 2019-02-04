function E = focus_fast_E(obj,f,phase,varargin)
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
            dphi = phi(2,1)-phi(1,1);
            
            theta = asin(obj.r/f);
            dr = obj.r(1,2)-obj.r(1,1);
            R = dr*size(obj.r,2);
            dtheta = ones(size(obj.r,1),1)*(asin((dr:dr:R)/f)-asin((0:dr:R-dr)/f));
            
            utheta(:,:,1) = cos(phi).*cos(theta);
            utheta(:,:,2) = sin(phi).*cos(theta);
            utheta(:,:,3) = -sin(theta);
            uphi(:,:,1) = -sin(phi);
            uphi(:,:,2) = cos(phi);
            uphi(:,:,3) = zeros(size(phi));
            
            k = 2*pi*sqrt(er2*mr2)/obj.lambda0;
            
            n1 = sqrt(er1);
            n2 = sqrt(er2);
            Et = (repmat(obj.Ephi,[1,1,3]).*uphi + repmat(obj.Er,[1,1,3]).*utheta).*(sqrt(n1/n2)*sqrt(cos(repmat(theta,[1,1,3]))));
            Et = reshape(Et,size(Et,1)*size(Et,2),3);
            thetaC = sin(theta).*dphi.*dtheta;
            
            
            E = zeros(7,3);
            
            for n = 1:7
                phaseTemp = phase(:,:,n);
                E(n,:) = sum(phaseTemp(:).*thetaC(:).*[Et(:,1),Et(:,2),Et(:,3)]);
            end
            E = 1i*k*f*exp(-1i*k*f)/(2*pi)*E;

        end