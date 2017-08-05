function B = focus_fast_B(ef,E,varargin)
            % B Magnetic field [T]
            %
            % B = B(EF,R) calculates the magnetic field at positions R (Point).
            %   B is a ComplexVector.
            %
            % B = B(EF,R,'dr',DR) sets the increment in the calcualtion of
            %   the derivatives to DR [default = 1e-12 m]. 
            %
            % See also EField, Point, ComplexVector.
            

            % increment [m]
            dr = 1e-12;
            for n = 1:2:length(varargin)
                if strcmpi(varargin{n},'dr')
                    dr = varargin{n+1};
                    Check.isnumeric('dr must be a positive real number',dr,'>',0)
                end
            end
            
            dxE = (E(3,:) - E(2,:))./(2*dr);
            dyE = (E(5,:) - E(4,:))./(2*dr);
            dzE = (E(7,:) - E(6,:))./(2*dr);
            
            B = -1i/ef.omega() * [ ...
                dyE(3)-dzE(2), ...
                -dxE(3)+dzE(1), ...
                dxE(2)-dyE(1) ...
                ];