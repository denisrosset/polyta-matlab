classdef (InferiorClasses = {?sym, ?double}) Polyta
% Describes a polytope with rational coefficents.
%
% Polyta is a class representing a polytope with rational
% coefficents, in the two standard representations
% (intersection of half-spaces and convex hull of set of points).
%
% All computations are exact.
    properties
        % INEQUALITIES
        A = [];
        b = [];
        % EQUALITIES
        Aeq = [];
        beq = [];
        % EXTREMAL POINTS
        E = []; % list of extremal points, one per row, rational
                % coefficients
        R = []; % list of rays, one per row, rational coefficents
        % FLAGS (not used now)
        Hminrep = false;
        Vminrep = false;
        % CACHED COMPUTATIONS
        cache = struct('x', []);
    end
    methods % Specialized methods with PolytaLibrary intelligent
            % selection algorithm
        function P = vertexEnumeration(P)
            assertHrepr(P);
            LM = PolytaLibraryManager.instance();
            for lib = LM.activated_libs
                if lib{:}.canDoVertexEnumeration ...
                        && (lib{:}.canHandleDegenerateH || ...
                            isequal(P.isfulldim(true), true)) ...
                        && (lib{:}.canHandleUnboundedH || ...
                            isequal(P.isbounded(true), true)) ...
                        && (~lib{:}.mustHaveValidPoint || ...
                            isequal(P.hasvalidpoint, true)) ...
                        && (lib{:}.canHandleOriginOutsideH || ...
                            isequal(P.hasorigin(true), true))
                    E = [];
                    R = [];
                    try
                        if LM.options.verbose
                            disp(['Vertex enumeration with ' ...
                                  class(lib{:})]);
                        end
                        [E, R] = ...
                            lib{:}.vertexEnumeration(P.A, P.b, ...
                                                     P.Aeq, P.beq, ...
                                                     P.cache.x, ...
                                                     LM.options);
                    catch me
                        warning(['Error using ' class(lib{:})]);
                        E = [];
                        R = [];
                    end
                    if ~isequal(E, []) || ~isequal(R, [])
                        P.E = E;
                        P.R = R;
                        P.Vminrep = true;
                        return
                    end
                end
            end
            error(['No library available to compute vertices of this ' ...
                   'polytope.']);
        end
        function P = facetEnumeration(P)
            assertVrepr(P);
            LM = PolytaLibraryManager.instance();
            for lib = LM.activated_libs
                if lib{:}.canDoVertexEnumeration ...
                        && (lib{:}.canHandleDegenerateV || ...
                            isequal(P.isfulldim(true), true)) ...
                        && (lib{:}.canHandleUnboundedV || ...
                            isequal(P.isbounded(true), true)) ...
                        && (lib{:}.canHandleOriginOutsideV || ...
                            isequal(P.hasorigin(true), true))
                    A = [];
                    b = [];
                    Aeq = [];
                    beq = [];
                    x = [];
                    try
                        if LM.options.verbose
                            disp(['Facet enumeration with ' ...
                                  class(lib{:})]);
                        end
                        [A, b, Aeq, beq, x] = ...
                            lib{:}.facetEnumeration(P.E, P.R, ...
                                                     LM.options);
                    catch me
                        warning(['Error using ' class(lib{:})]);
                        A = [];
                        b = [];
                        Aeq = [];
                        beq = [];
                        x = [];
                    end
                    if ~isequal(A, []) || ~isequal(Aeq, [])
                        P.A = A;
                        P.b = b;
                        P.Aeq = Aeq;
                        P.beq = beq;
                        P.Hminrep = true;
                        if ~isequal(x, [])
                            P.cache.x = x;
                        end
                        return
                    end
                end
            end
            error(['No library available to compute facets of this ' ...
                   'polytope.']);
        end
        function P = findInteriorPoint(P)
            LM = PolytaLibraryManager.instance();
            error('Not implemented');
        end
        function P = linearProgramming(P, f)
            LM = PolytaLibraryManager.instance();
            error('Not implemented');
        end
        function P = projection(P, dims)
            LM = PolytaLibraryManager.instance();
            error('Not implemented');
        end
        function P = removeRedundantH()
            error('Not implemented');
        end
        function P = removeRedundantV()
            error('Not implemented');
        end
    end
    methods
        function P = Polyta(Hminrep, A, b, Aeq, beq, Vminrep, E, R)
            P.Hminrep = Hminrep;
            P.A = A;
            P.b = b;
            P.Aeq = Aeq;
            P.beq = beq;
            P.Vminrep = Vminrep;
            P.E = E;
            P.R = R;
        end
    end
    methods
        function assert(P)
        % Verifies that the polytope P is well formed.
            if ~isequal(P.cache.x, []) && P.hasHrepr
                assert(isinside(P, P.x));
            end
            % TODO: add more safety checks
        end
        function assertHrepr(P)
            if ~hasHrepr(P)
                error(['The polytope must have an H representation for this ' ...
                       'algorithm to work.']);
            end
        end
        function assertVrepr(P)
            if ~hasVrepr(P)
                error(['The polytope must have an V representation for this ' ...
                       'algorithm to work.']);
            end
        end
        function assertbounded(P)
            if ~isbounded(P)
                error(['This algorithm works only on bounded ' ...
                       'polytopes. Maybe you can remove the ' ...
                       'unconstrained dimensions.']);
            end
        end
    end
    methods
        function b = hasHrepr(P)
            b = ~isequal(P.A, []) || ~isequal(P.Aeq, []);
        end
        function b = hasVrepr(P)
            b = ~isequal(P.E, []) || ~isequal(P.R, []);
        end
        function b = hasvalidpoint(P)
            b = ~isequal(P.cache.x, []);
        end
        function b = isunbounded(P, fast)
            if nargin < 2
                fast = false;
            end
            b = ~P.isbounded(fast);
        end
        function b = isbounded(P, fast)
            if nargin < 2
                fast = false;
            end
            if P.hasVrepr
                b = isequal(P.R, []);
                return
            end
            if fast
                b = [];
                return
            end
            error('TODO: implement this using Minkowski theorem.');
        end
        function b = isdegenerate(P, fast)
            if nargin < 2
                fast = false;
            end
            b = ~P.isfulldim(fast);
        end
        function b = isfulldim(P, fast)
            if nargin < 2
                fast = false;
            end
            if P.hasHrepr && P.Hminrep
                b = isequal(P.Aeq, []);
                return
            end
            if fast
                b = [];
                return
            end
            error('TODO: implement using rank.');
        end
        function b = hasorigin(P, fast)
        % Checks if the origin (0,0...0) is inside the polytope P.
            if nargin < 2
                fast = false;
            end
            b = P.isinside(P.origin, fast);
        end
        function b = isinside(P, x0, fast)
        % Checks if the given point x0 lies in the interior of the
        % polytope P.
            if nargin < 3
                fast = false;
            end
            if P.hasHrepr
                b = true;
                if ~isequal(P.A, [])
                    b = b && all(P.A * x0(:) - P.b <= 0);
                end
                if ~isequal(P.Aeq, [])
                    b = b && all(P.Aeq * x0(:) - P.beq == 0);
                end
                return
            end
            if fast
                b = [];
                return
            end
            error('TODO: implement using linear programming.');
        end
        function b = hasinteriorpoint(P)
            b = ~isequal(P.cache.x, []);
        end
    end
    methods % from mpt/polytope
        function [xc, rc2, lambda] = chebyball(P)
        % [xc, rc2, lambda] = chebyball(P)
        % Computes center and squared radius of the largest ball inscribed
        % in a polytope
            error('Not implemented');
        % TODO: using 10.1.1.47.9124.pdf
        end
        
        function nx = dimension(P)
        % nx = dimension(P)
        % Returns dimension of the given polytope P.
            nx = max([size(P.A, 2) size(P.Aeq, 2) size(P.E, 2) ...
                        size(P.R, 2)]);
        end
        function x = origin(P)
        % Origin point in R^n, where n is the dimension of P.
            x = rational(zeros(1, P.dimension), 1);
        end
    end
    methods
        function x = random_point_boundary(P, d)
        % x = random_point_boundary(P, denom)
        % Generate a random point on the boundary of P
        %
        % INPUTS
        %     P   - Polyta, bounded, with both representations
        %     d   - Denominator for rational random numbers,
        %           default d=65536
        % OUTPUTS
        %     x   - Random point on P surface
            assertbounded(P);
            assertHrepr(P);
            assertVrepr(P);
            if nargin < 2
                d = 65536;
            end
            % select facet at random
            facet = randi(size(P.A,1));
            % find extreme points who are on selected facet
            onfacet = (P.A(facet, :) * P.E' - P.b(facet)) == 0;
            coeff = rational(zeros(1, size(P.E, 1)));
            % point is a convex combination with random coefficents
            % of extremal points on facet
            coeff(onfacet) = rational(randi(d, 1, ...
                                            size(find(onfacet), 2)))/d;
            x = coeff * P.E/sum(coeff);
        end
        function P = Hrepr(P)
            if ~hasHrepr(P)
                P = P.facetEnumeration;
            end
        end
        function P = Vrepr(P)
            if ~hasVrepr(P)
                P = P.vertexEnumeration;
            end
        end
        function [P, Aeq, beq] = affinespace(P)
        % [P, Aeq, beq] = affinespace(P)
        % Find the affine space spanned by P
        %
        % INPUTS
        %     P   - Polyta
        % OUTPUTS
        %     P   - Polyta with cached H representation
        %     Aeq,
        %     beq - Linear system describing affine space spanned
        %           by P, that is all the points x such that 
        %                   Aeq * x == b
            if P.hasHrepr && P.Hminrep
                Aeq = P.Aeq;
                beq = P.beq;
            else
                if ~P.hasVrepr
                    P = Vrepr(P);
                end
                error('Not implemented with V repr.');
                % TODO: implement using linear algebra
            end
        end
        function [P, A, b, Aeq, beq] = facets(P)
        % [P, A, b, Aeq, beq] = facets(P)
        % Find all facets of polytope P.
            P = Hrepr(P);
            assert(~~P.Hminrep);
            % TODO: implement even if H representation is not minimal
            A = P.A;
            b = P.b;
            Aeq = P.Aeq;
            beq = P.beq;
        end
        function [P, x] = feasible(P, x)
        % [P, x] = feasible(P)
        % Find a feasible point for polytope P
        %
        % INPUTS
        %     P   - Polyta
        % OUTPUTS
        %     P   - Polyta with cached feasible point
        %     x   - Feasible point satisfying the constraints
        % TODO: redo this function
            if hasVrepr(P)
                if ~isequal(P.E, [])
                    x = P.E(1,:);
                else
                    x = P.origin;
                end
            else
                x = P.findInteriorPoint();
            end
        end
        function Pproj = project(P, dims)
        % Pproj = project(P, dims)
        % PROJECT projects the polytope on specified dimensions.
        %
        % Preserves both representations if they are available.
        %error('Does not work without polytalib global.'); % TODO
            global polytalib;
            Pproj = P;
            Pproj.cache.x = [];
            if hasVrepr(Pproj)
                Pproj.E(:, dims) = 0;
                if ~isequal(Pproj.R, [])
                    Pproj.R(:, dims) = 0;
                end
            end
            if hasHrepr(Pproj)
                [Pproj.A, Pproj.b, Pproj.Aeq, Pproj.beq, Pproj.cache.x] = ...
                    polytalib.projection(P.A, P.b, ...
                                         P.Aeq, P.beq, ...
                                         P.cache.x, dims);
            end
        end
        function Pcut = slice(P, cut_dim, x)
        % SLICE cuts a polytope at specified dimensions
        %
        % Pcut = slice(P, cut_dim)
        %
        % INPUTS
        %     P       - Polyta
        %     cut_dim - vector of dimensions to cut at 0
        % OUTPUTS
        %     Pcut    - Cut Polyta
            cut_dim = cut_dim(:);
            Pcut = facets(P);
            Aeqs = zeros(size(cut_dim,1), size(Pcut.A, 2));
            Aeqs(1:end, cut_dim) = eye(size(cut_dim, 1));
            Aeq = Pcut.Aeq;
            beq = Pcut.beq;
            Pcut.Aeq = [Aeq; Aeqs];
            Pcut.beq = [beq; zeros(size(cut_dim))];
            Pcut.cache.x = x;
            Pcut.E = [];
            Pcut.Hminrep = false;
        end
        function [P, E, R] = extreme(P)
        % [P, E, R] = extreme(P)
            P = Vrepr(P);
            E = P.E;
            R = P.R;
        end
    end
    methods(Static)
        function P = from_minimal_facets(A, b, Aeq, beq)
            if nargin < 3
                Aeq = [];
                beq = [];
            end
            if ~isequal(A, [])
                A = rational(A);
            end
            if ~isequal(b, [])
                b = rational(b);
            end
            if ~isequal(Aeq, [])
                Aeq = rational(Aeq);
            end
            if ~isequal(beq, [])
                beq = rational(beq);
            end
            P = Polyta(true, A, b, Aeq, beq, [], [], []);
        end
        function P = from_redundant_facets(A, b, Aeq, beq)
            if nargin < 3
                Aeq = [];
                beq = [];
            end
            if ~isequal(A, [])
                A = rational(A);
            end
            if ~isequal(b, [])
                b = rational(b);
            end
            if ~isequal(Aeq, [])
                Aeq = rational(Aeq);
            end
            if ~isequal(beq, [])
                beq = rational(beq);
            end
            P = Polyta(false, A, b, Aeq, beq, [], [], []);
        end
        function P = from_minimal_vertices(E, R)
            if nargin < 2
                R = [];
            end
            if ~isequal(E, [])
                E = rational(E);
            end
            if ~isequal(R, [])
                R = rational(R);
            end
            P = Polyta([], [], [], [] ,[], true, E, R);
        end
        function P = from_redundant_vertices(E, R)
            if nargin < 2
                R = [];
            end
            if ~isequal(E, [])
                E = rational(E);
            end
            if ~isequal(R, [])
                R = rational(R);
            end
            P = Polyta([], [], [], [] ,[], false, E, R);
        end
        function P = from_ext_file(filename, minimal)
            if nargin < 2
                minimal = false;
            end
            fid = fopen(filename, 'r');
            [E, R] = PolytaCddLibrary.parse_ext(fid);
            fclose(fid);
            P = Polyta([], [], [], [], [], minimal, E, R);
        end
        function P = from_ine_file(filename, minimal)
            if nargin < 2
                minimal = false;
            end
            fid = fopen(filename, 'r');
            [A, b, Aeq, beq] = PolytaCddLibrary.parse_ine(fid);
            fclose(fid);
            P = Polyta(minimal, A, b, Aeq, beq, [], [], []);
        end
        function P = from_poi_file(filename, minimal)
            if nargin < 2
                minimal = false;
            end
            fid = fopen(filename, 'r');
            [E, R] = PolytaPortaLibrary.parse_poi(fid);
            fclose(fid);
            P = Polyta([], [], [], [], [], minimal, E, R);
        end
        function P = from_ieq_file(filename, minimal)
            if nargin < 2
                minimal = false;
            end
            fid = fopen(filename, 'r');
            [A, b, Aeq, beq, x] = parse_ieq(fid);
            fclose(fid);
            P = Polyta(minimal, A, b, Aeq, beq, [], [], []);
            P.cache.x = x;
        end
    end
end
