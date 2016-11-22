classdef PolytaLrsLibrary < PolytaLibrary
    properties
        canDoVertexEnumeration = true;
        canDoFacetEnumeration = true;
        canDoLinearProgramming = true;
        canFindInteriorPoint = true;
        canDoSliceAtOrigin = false;
        canDoSliceAnywhere = false;
        canDoProjection = false;
        canRemoveRedundancyInH = true;
        canRemoveRedundancyInV = true;
        
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
            p.cmd = ['lrs ' p.infilename ' > ' p.outfilename];
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
            p.cmd = ['lrs ' p.infilename ' > ' p.outfilename];
            [A, b, Aeq, beq] = ...
                PolytaCddLibrary.facetEnumerationDriver(E, R, p);
            x = [];
        end
        function [x, fval, lambda] = ...
                linearProgramming(L, f, A, b, Aeq, beq, x0, options)
            tmp = tempname;
            if nargin < 8
                p = struct();
            else
                p = options;
            end
            p.infilename = [tmp '.ine'];
            p.outfilename = [tmp '.lps'];
            p.cmd = ['lrs ' p.infilename ' > ' p.outfilename];
            p.inverse_lambda = true;
            [x, fval, lambda] = ...
                PolytaCddLibrary.linearProgrammingDriver(f, A, b, ...
                                                         Aeq, beq, p);
            exitflag = 1; % TODO : check for errors
            output = [];
        end
        function x = findInteriorPoint(L, A, b, Aeq, beq, options)
        % BUGS : does not work with unbounded polytopes, because
        % LRS can output a ray before an extreme point
            tmp = tempname;
            if nargin < 6
                p = struct();
            else
                p = options;
            end
            p.infilename = [tmp '.ine'];
            p.outfilename = [tmp '.ext'];
            p.inoptions = {'maxoutput 1'};
            p.cmd = ['lrs ' p.infilename ' > ' p.outfilename];
            [E, R] = PolytaCddLibrary.vertexEnumerationDriver(A, b, ...
                                                              Aeq, beq, p);
            x = E(1,:);
        end
        function [E, R] = removeRedundantV(L, E, R, options)
            tmp = tempname;
            infilename = [tmp '.ext'];
            outfilename = [tmp '.ext.ext'];
            fid = fopen(infilename, 'w');
            PolytaCddLibrary.write_ext(fid, E, R);
            fclose(fid);
            cmd = ['redund ' p.infilename ' > ' p.outfilename];
            if nargin < 4 || ~isfield(options, 'verbose') || options.verbose
                status = unix(cmd);
            else
                [status, out] = unix(cmd);
            end
            assert(status == 0);
            fid = fopen(p.outfilename, 'r');
            [E, R] = PolytaCddLibrary.parse_ext(fid);
            fclose(fid);
            if nargin < 4 || ~isfield(options, 'delete') || options.delete
                delete(infilename);
                delete(outfilename);
            end
        end
        function [A, b, Aeq, beq] = removeRedundantH(L, A, b, ...
                                                     Aeq, beq, options)
            tmp = tempname;
            infilename = [tmp '.ine'];
            outfilename = [tmp '.ine.ine'];
            fid = fopen(infilename, 'w');
            PolytaCddLibrary.write_ine(fid, A, b, Aeq, beq);
            fclose(fid);
            cmd = ['redund ' p.infilename ' > ' p.outfilename];
            if nargin < 6 || ~isfield(options, 'verbose') || options.verbose
                status = unix(cmd);
            else
                [status, out] = unix(cmd);
            end
            assert(status == 0);
            fid = fopen(p.outfilename, 'r');
            [A, b, Aeq, beq] = PolytaCddLibrary.parse_ine(fid);
            fclose(fid);
            if nargin < 6 || ~isfield(options, 'delete') || options.delete
                delete(infilename);
                delete(outfilename);
            end
        end
    end
end
