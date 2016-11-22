classdef PolytaLibrary
    properties(Abstract)
        canDoVertexEnumeration;
        canDoFacetEnumeration;
        canDoLinearProgramming;
        canFindInteriorPoint;
        canDoSliceAtOrigin;
        canDoSliceAnywhere;
        canDoProjection;
        canRemoveRedundancyInH;
        canRemoveRedundancyInV;
        
        mustHaveValidPoint;
        canHandleUnboundedV;
        canHandleUnboundedH;
        canHandleDegenerateH;
        canHandleDegenerateV;
        canHandleOriginOutsideH;
        canHandleOriginOutsideV;
    end
    methods
        function ok = isAvailable(L)
            A = rational([-1  0
                           1  0
                           0 -1
                           0  1]);
            b = rational([1;1;1;1]);
            x = rational([0 0]);
            E = [];
            R = [];
            options = struct();
            options.verbose = false;
            options.delete = true;
            try
                [E, R] = L.vertexEnumeration(A, b, [], [], x, options);
            catch me
                % do nothing
            end
            ok = isequal(sortrows(double(E)), ...
                         [-1 -1
                          -1  1
                           1 -1
                           1  1]);
        end
        function [E, R] = vertexEnumeration(L, A, b, Aeq, beq, x, options)
            error('Not implemented.');
        end
        function [x] = findInteriorPoint(L, A, b, Aeq, beq, options)
            error('Not implemented.');
        end
        function [A, b, Aeq, beq, x] = facetEnumeration(L, E, R, options)
            error('Not implemented.');
        end
        function [x, fval, lambda] = ...
                linearProgramming(L, f, A, b, Aeq, beq, x0, options)
            % f is a row vector
            error('Not implemented.');
        end
        function [E, R] = sliceAtOrigin(L, E, R, dims, options)
            error('Not implemented.');
        end
        function [E, R] = slice(L, E, R, dims, val, options)
            error('Not implemented.');
        end
        function [A, b, Aeq, beq, x] = projection(L, A, b, ...
                                                  Aeq, beq, ...
                                                  x, dims, options)
            error('Not implemented.');
        end
        function [E, R] = removeRedundantV(L, E, R, options)
            error('Not implemented.');
        end
        function [A, b, Aeq, beq] = removeRedundantH(L, A, b, ...
                                                     Aeq, beq, options)
            error('Not implemented.');
        end
    end
    methods(Static)
        function [n, d] = parse_ratnum(c)
        % Parse a string representing a rational number.
        % 
        % INPUTS
        %     c   - char string representing a number of the form
        %           -2/3, +342/54, +2, -1, 2 or 0
        % OUTPUTS
        %     n   - rational number numerator, integer stored as double
        %     d   - rational number denominator, integer stored as double
            assert(isa(c, 'char'));
            c = strtrim(c);
            if isequal(c, '+')
                n = 1;
                d = 1;
                return
            end
            if isequal(c, '-')
                n = -1;
                d = 1;
                return
            end
            [s,e,t,m,toks] = regexp(c, '([+-]?\s*\d+)(\s*/\s*\d+)?', 'once');
            n = str2num(toks{1});
            if isequal(toks{2}, '')
                d = 1;
                return
            end
            d = str2num(strrep(toks{2}, '/', ''));
            assert(d~=0);
        end
    end
end
