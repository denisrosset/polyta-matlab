classdef PolytaPplLibrary < PolytaLibrary
    properties
        canDoVertexEnumeration = true;
        canDoFacetEnumeration = true;
        canDoLinearProgramming = false;
        canFindInteriorPoint = false;
        canDoSliceAtOrigin = false;
        canDoSliceAnywhere = false;
        canDoProjection = false;
        canRemoveRedundancyInH = false;
        canRemoveRedundancyInV = false;

        mustHaveValidPoint = false;
        canHandleUnboundedV = true;
        canHandleUnboundedH = true;
        canHandleDegenerateH = true;
        canHandleDegenerateV = true;
        canHandleOriginOutsideH = true;
        canHandleOriginOutsideV = true;
    end
    methods
        function [E, R] = vertexEnumeration(L, A, b, Aeq, beq, x, options)
            tmp = tempname;
            if nargin < 7
                p = struct();
            else
                p = options;
            end
            p.infilename = [tmp '.ine'];
            p.outfilename = [tmp '.ext'];
            p.cmd = ['ppl_lcdd -o ' p.outfilename ' ' p.infilename];
            [E, R] = PolytaCddLibrary.vertexEnumerationDriver(A, b, ...
                                                              Aeq, ...
                                                              beq, p);
        end
        function [A, b, Aeq, beq, x] = facetEnumeration(L, E, R, options)
                        tmp = tempname;
            if nargin < 4
                p = struct();
            else
                p = options;
            end
            p.infilename = [tmp '.ext'];
            p.outfilename = [tmp '.ine'];
            p.cmd = ['ppl_lcdd -o ' p.outfilename ' ' p.infilename];
            [A, b, Aeq, beq] = PolytaCddLibrary.facetEnumerationDriver(E, R, p);
            x = [];
        end
    end
end
