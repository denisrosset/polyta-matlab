classdef PolytaCddLibrary < PolytaLibrary
    properties
        canDoVertexEnumeration = true;
        canDoFacetEnumeration = true;
        canDoLinearProgramming = true;
        canFindInteriorPoint = true;
        canDoSliceAtOrigin = false;
        canDoSliceAnywhere = false;
        canDoProjection = false;
        canRemoveRedundancyInH = false;
        canRemoveRedundancyInV = false;
        
        mustHaveValidPoint = false;
        canHandleUnboundedV = true;
        canHandleUnboundedH = true;
        canHandleDegenerateH = false;
        canHandleDegenerateV = false;
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
            p.inoptions = {'dynout_off'};
            p.outfilename = [tmp '.ext'];
            p.logfilename = [tmp '.ddl'];
            p.cmd = ['cddr+ ' p.infilename];
            [E, R] = PolytaCddLibrary.vertexEnumerationDriver(A, b, ...
                                                              Aeq, beq, ...
                                                              p);
        end
        function [A, b, Aeq, beq, x] = facetEnumeration(L, E, R, ...
                                                        options)
            tmp = tempname;
            if nargin < 4
                p = struct();
            else
                p = options;
            end
            p.infilename = [tmp '.ext'];
            p.inoptions = {'dynout_off'};
            p.outfilename = [tmp '.ine'];
            p.logfilename = [tmp '.ddl'];
            p.cmd = ['cddr+ ' p.infilename];
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
            p.inoptions = {'dynout_off'};
            p.outfilename = [tmp '.lps'];
            p.logfilename = [tmp '.ddl'];
            p.cmd = ['cddr+ ' p.infilename];
            p.inverse_lambda = true;
            [x, fval, lambda] = ...
                PolytaCddLibrary.linearProgrammingDriver(f, A, b, ...
                                                         Aeq, beq, p);
            exitflag = 1; % todo : check for errors
            output = [];
        end
        function x = findInteriorPoint(L, A, b, Aeq, beq, options)
            tmp = tempname;
            if nargin < 6
                p = struct();
            else
                p = options;
            end
            p.infilename = [tmp '.ine'];
            p.inoptions = {'dynout_off'};
            p.outfilename = [tmp '.lps'];
            p.logfilename = [tmp '.ddl'];
            p.cmd = ['cddr+ ' p.infilename];
            x = PolytaCddLibrary.findInteriorPointDriver(A, b, Aeq, ...
                                                         beq, p);
            % if x(end) == 0, the polytope has no interior point
            % if x(end) < 0, the polytope is empty
            if x(end) < 0
                x = [];
            else
                x = x(1:end-1);
            end
        end
    end
    methods(Static)
        function [E, R] = vertexEnumerationDriver(A, b, Aeq, beq, p)
            fid = fopen(p.infilename, 'w');
            PolytaCddLibrary.write_ine(fid, A, b, Aeq, beq);
            if isfield(p, 'inoptions')
                for l = p.inoptions
                    fprintf(fid, [char(l) '\n']);
                end
            end
            fclose(fid);
            if nargin < 5 || ~isfield(p, 'verbose') || p.verbose
                status = unix(p.cmd);
            else
                [status, out] = unix(p.cmd);
            end
            assert(status == 0);
            fid = fopen(p.outfilename, 'r');
            [E, R] = PolytaCddLibrary.parse_ext(fid);
            fclose(fid);
            if ~isfield(p, 'delete')
                delete(p.infilename);
                delete(p.outfilename);
                if isfield(p, 'logfilename')
                    delete(p.logfilename);
                end
            end
        end
        function [A, b, Aeq, beq] = facetEnumerationDriver(E, R, p)
            fid = fopen(p.infilename, 'w');
            PolytaCddLibrary.write_ext(fid, E, R);
            if isfield(p, 'inoptions')
                for l = p.inoptions
                    fprintf(fid, [char(l) '\n']);
                end
            end
            fclose(fid);
            if ~isfield(p, 'verbose') || p.verbose
                status = unix(p.cmd);
            else
                [status, out] = unix(p.cmd);
            end
            assert(status == 0);
            fid = fopen(p.outfilename, 'r');
            [A, b, Aeq, beq] = PolytaCddLibrary.parse_ine(fid);
            fclose(fid);
            if ~isfield(p, 'delete')
                delete(p.infilename);
                delete(p.outfilename);
                if isfield(p, 'logfilename')
                    delete(p.logfilename);
                end
            end
        end

        function [x, fval, lambda] = ...
                linearProgrammingDriver(f, A, b, Aeq, beq, p)
            fid = fopen(p.infilename, 'w');
            PolytaCddLibrary.write_ine(fid, A, b, Aeq, beq);
            fprintf(fid, 'minimize\n%d', 0);
            for i = 1:size(f, 2)
                fprintf(fid, ' %s', char(f(i)));
            end
            if isfield(p, 'inoptions')
                for l = p.inoptions
                    fprintf(fid, [char(l) '\n']);
                end
            end
            fclose(fid);
            if ~isfield(p, 'verbose') || p.verbose
                status = unix(p.cmd);
            else
                [status, out] = unix(p.cmd);
            end
            assert(status == 0);
            fid = fopen(p.outfilename, 'r');
            [x, fval, lmb] = PolytaCddLibrary.parse_lps(fid);
            if isfield(p, 'inverse_lambda') && p.inverse_lambda
                lmb = -lmb;
            end
            lambda = struct();
            lambda.eqlin = lambda(size(A, 1) + 1:...
                                  size(A, 1) + size(Aeq, 1));
            lambda.ineqlin = lambda(1:size(A, 1));
            fclose(fid);
            if ~isfield(p, 'delete')
                delete(p.infilename);
                delete(p.outfilename);
                if isfield(p, 'logfilename')
                    delete(p.logfilename);
                end
            end
        end
        function x = findInteriorPointDriver(A, b, Aeq, beq, p)
            fid = fopen(p.infilename, 'w');
            PolytaCddLibrary.write_ine(fid, A, b, Aeq, beq);
            fprintf(fid, 'find_interior\n');
            if isfield(p, 'inoptions')
                for l = p.inoptions
                    fprintf(fid, [char(l) '\n']);
                end
            end
            fclose(fid);
            if ~isfield(p, 'verbose') || p.verbose
                status = unix(p.cmd);
            else
                [status, out] = unix(p.cmd);
            end
            assert(status == 0);
            fid = fopen(p.outfilename, 'r');
            x = PolytaCddLibrary.parse_lps(fid);
            fclose(fid);
            if ~isfield(p, 'delete')
                delete(p.infilename);
                delete(p.outfilename);
                if isfield(p, 'logfilename')
                    delete(p.logfilename);
                end
            end
        end
        function [x, fval, lambda] = parse_lps(fid)
            delim = [char(9) ' ' ':' '|' char(10) char(13)];
            state = 0; % state = 0, beginning of file
                       % state = 1, encountered begin
                       % state = 2, reading primal_solution
                       % state = 3, reading dual_solution
                       % state = 4, end read
            x = rational([]);
            lambda = rational([]);
            fval = [];
            while 1
                tl = fgetl(fid);
                if isa(tl, 'double') && tl == -1
                    assert(state == 4);
                    break
                end
                [word, remain] = strtok(tl, delim);
                switch state
                  case 0
                    if isequal(word, 'begin')
                        state = 1;
                    end
                  case {1, 2, 3}
                    switch word
                      case 'end'
                        state = 4;
                      case 'primal_solution'
                        state = 2;
                      case 'dual_solution'
                        state = 3;
                      case 'optimal_value'
                        [word, remain] = strtok(remain, delim);
                        [n, d] = PolytaLibrary.parse_ratnum(word);
                        fval = rational(n) / rational(d);
                      otherwise
                        switch state
                          case 1
                            error('Invalid identifier.');
                          case 2
                            index = str2num(word);
                            assert(isa(index, 'double'));
                            [word, remain] = strtok(remain, delim);
                            [n, d] = PolytaLibrary.parse_ratnum(word);
                            x(1, index) = rational(n) / rational(d);
                          case 3
                            index =  str2num(word);
                            assert(isa(index, 'double'));
                            [word, remain] = strtok(remain, delim);
                            [n, d] = PolytaLibrary.parse_ratnum(word);
                            lambda(1, index) = rational(n) / rational(d);
                        end
                    end
                  case 4
                    % do nothing
                end
            end
            if size(x, 2) == 0
                x = [];
            end
            if size(lambda, 2) == 0
                lambda = [];
            end
        end
        function [E, R] = parse_ext(fid)
            delim = [char(9) ' ' ':' '|' char(10) char(13)];
            state = 0; % state = 0, beginning of file
                       % state = 1, encountered begin, but not dimensions
                       % state = 2, reading points
                       % state = 3, waiting for end
                       % state = 4, end read
            N = 0;
            d = 0;
            MN = [];
            ND = [];
            vertex = [];
            while 1
                tl = fgetl(fid);
                if isa(tl, 'double') && tl == -1
                    assert(state == 4);
                    break
                end
                [word, remain] = strtok(tl, delim);
                switch state
                  case 0
                    switch word
                      case 'H-representation'
                        error(['Parsing a H-representation file with ' ...
                               'V-representation routine.']);
                      case 'begin'
                        state = 1;
                    end
                  case 1
                    if isequal(word, '*****')
                        N = [];
                    else
                        N = str2num(word);
                    end
                    [word, remain] = strtok(remain, delim);
                    d = str2num(word) - 1;
                    [word, remain] = strtok(remain, delim);
                    assert(isequal(word, 'rational') || ...
                           isequal(word, 'integer'));
                    if ~isequal(N, [])
                        MN = zeros(N, d);
                        MD = zeros(N, d);
                        vertex = zeros(N, 1) == 1;
                    end
                    i = 1;
                    state = 2;
                  case 2
                    if isequal(word, 'end')
                        assert(isequal(N, []));
                        N = i - 1;
                        state = 4;
                    else
                        vertex(i) = str2num(word) == 1;
                        for j = 1:d
                            [word, remain] = strtok(remain, delim);
                            [MN(i,j) MD(i,j)] = ...
                                PolytaLibrary.parse_ratnum(word);
                        end
                        if ~isequal(N, []) && i == N
                            state = 3;
                        else
                            i = i + 1;
                        end
                    end
                  case 3
                    assert(isequal(word, 'end'));
                    state = 4;
                  otherwise
                    % do nothing
                end
            end
            % post process output
            E = rational(MN(~~vertex, :)) ./ rational(MD(~~vertex, :));
            R = rational(MN(~vertex, :)) ./ rational(MD(~vertex, :));
            if size(E, 1) == 0
                E = [];
            end
            if size(R, 1) == 0
                R = [];
            end
        end
        function [A, b, Aeq, beq] = parse_ine(fid)
            delim = [char(9) ' ' ':' '|' char(10) char(13)];
            state = 0; % state = 0, beginning of file
                       % state = 1, encountered begin, but not dimensions
                       % state = 2, reading points
                       % state = 3, waiting for end
                       % state = 4, end read
            eq = [];
            N = 0;
            d = 0;
            MN = [];
            ND = [];
            while 1
                tl = fgetl(fid);
                if isa(tl, 'double') && tl == -1
                    assert(state == 4);
                    break
                end
                [word, remain] = strtok(tl, delim);
                switch state
                  case 0
                    switch word
                      case 'linearity'
                        assert(isequal(eq, []));
                        [n, remain] = strtok(remain, delim);
                        n = str2num(n);
                        eq = [];
                        for i = 1:n
                            [val, remain] = strtok(remain, delim);
                            eq = [eq str2num(val)];
                        end
                      case 'V-representation'
                        error(['Parsing a V-representation file with ' ...
                               'H-representation routine.']);
                      case 'begin'
                        state = 1;
                      otherwise
                        % do nothing
                    end
                  case 1
                    if isequal(word, '*****')
                        N = [];
                    else
                        N = str2num(word);
                    end
                    [word, remain] = strtok(remain, delim);
                    d = str2num(word) - 1;
                    [word, remain] = strtok(remain, delim);
                    assert(isequal(word, 'rational') || ...
                           isequal(word, 'integer'));
                    if ~isequal(N, [])
                        MN = zeros(N, d + 1);
                        MD = zeros(N, d + 1);
                    end
                    i = 1;
                    state = 2;
                  case 2
                    if isequal(word, 'end')
                        assert(isequal(N, []));
                        N = i - 1;
                        state = 4;
                    else
                        remain = tl;
                        for j = 1:d+1
                            [word, remain] = strtok(remain, delim);
                            [MN(i,j) MD(i,j)] = ...
                                PolytaLibrary.parse_ratnum(word);
                        end
                        if ~isequal(N, []) && i == N
                            state = 3;
                        else
                            i = i + 1;
                        end
                    end
                  case 3
                    assert(N == i);
                    assert(isequal(word, 'end'));
                    state = 4;
                  otherwise
                    % do nothing
                end
            end
            ineq = setdiff(1:N, eq);
            if size(eq, 2) > 0
                Aeq = -rational(MN(eq, 2:end)) ./ ...
                      rational(MD(eq, 2:end));
                beq = rational(MN(eq, 1)) ./ rational(MD(eq, 1));
            else
                Aeq = [];
                beq = [];
            end
            A = -rational(MN(ineq, 2:end)) ./ ...
                rational(MD(ineq, 2:end));
            b = rational(MN(ineq, 1)) ./ rational(MD(ineq, 1));
            if size(A, 1) == 0
                A = [];
            end
            if size(b, 1) == 0
                b = [];
            end
        end            
        function write_ine(fid, A, b, Aeq, beq)
        % PolytaCddLibrary.write_ine(fid, A, b, Aeq, beq)
        % Write .ine cdd file (1997 version)
        %
        % INPUTS
        %     fid     - opened file descriptor (write mode)
        %     A,b     - linear inequalities A * x <= b
        %     Aeq,beq - linear equalities Aeq * x == beq
            fprintf(fid, 'H-representation\n');
            d = size(A, 2);
            if size(Aeq, 1) > 0
                M1 = [b -A];
                M2 = [beq -Aeq];
                M = [M1;M2];
                inds = (size(M1, 1)+1) : (size(M1, 1)+size(M2, 1));
                fprintf(fid, 'linearity %d', size(inds, 2));
                for i = inds
                    fprintf(fid, ' %d', i);
                end
                fprintf(fid, '\n');
            else
                M = [b -A];
            end
            m = size(M, 1);
            fprintf(fid, 'begin\n');
            fprintf(fid, '%d %d rational\n', m, d + 1);
            for i = 1:size(M, 1)
                str = '';
                delim = '';
                for j = 1:size(M, 2)
                    str = [str delim char(M(i,j))];
                    delim = ' ';
                end
                str = [str '\n'];
                fprintf(fid, str);
            end
            fprintf(fid, 'end\n');
        end
        function write_ext(fid, E, R)
        % PolytaCddLibrary.write_ext(fid, E, R)
        % Write a .ext cdd file (1997 version)
        %
        % INPUTS
        %      fid   - opened file descriptor (write mode)
        %      E     - matrix of extremal points (one per row)
        %      R     - matrix of rays (one per row)
            fprintf(fid, 'V-representation\n');
            fprintf(fid, 'begin\n');
            n = size(E, 1);
            s = size(R, 1);
            d = size(E, 2);
            E1 = [ones(size(E,1), 1) E];
            R1 = [zeros(size(R, 1), 1) R];
            M = [E1;R1];
            fprintf(fid, '%d %d rational\n', n + s, d + 1);
            for i = 1:size(M, 1)
                fmt = '%s';
                for j = 1:size(M, 2)
                    fprintf(fid, fmt, char(M(i,j)));
                    fmt = ' %s';
                end
                fprintf(fid, '\n');
            end
            fprintf(fid, 'end\n');
        end
    end
end
