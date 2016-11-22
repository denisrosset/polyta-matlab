classdef PolytaPortaLibrary < PolytaLibrary
% BUGS
%
% PORTA must be compiled with support for long filenames (array
% sizes in function main of porta.c must be > 100 instead of 20).
    properties
        canDoVertexEnumeration = true;
        canDoFacetEnumeration = true;
        canDoLinearProgramming = false;
        canFindInteriorPoint = false;
        canDoSliceAtOrigin = true;
        canDoSliceAnywhere = false;
        canDoProjection = true;
        canRemoveRedundancyInH = false;
        canRemoveRedundancyInV = false;
        
        mustHaveValidPoint = true;
        canHandleUnboundedV = true;
        canHandleUnboundedH = true;
        canHandleDegenerateH = true;
        canHandleDegenerateV = true;
        canHandleOriginOutsideH = true;
        canHandleOriginOutsideV = true;
    end
    methods
        function [E, R] = vertexEnumeration(L, A, b, Aeq, beq, x, options)
            fn = [tempname '.ieq'];
            f = fopen(fn, 'w');
            PolytaPortaLibrary.write_ieq(f, A, b, Aeq, beq, x);
            fclose(f);
            if nargin < 7 || ~isfield(options, 'verbose') || options.verbose
                result = unix(['traf -o ' fn]);
            else
                [result, out] = unix(['traf -o ' fn]);
            end
            assert(result == 0);
            f = fopen([fn '.poi'], 'r');
            [E, R] = PolytaPortaLibrary.parse_poi(f);
            fclose(f);
            if nargin < 7 || ~isfield(options, 'delete') || options.delete
                delete(fn);
                delete([fn '.poi']);
            end
        end
        function [A, b, Aeq, beq, x] = facetEnumeration(L, E, R, options)
            fn = [tempname '.poi'];
            f = fopen(fn, 'w');
            PolytaPortaLibrary.write_poi(f, E, R);
            fclose(f);
            if nargin < 4 || ~isfield(options, 'verbose') || options.verbose
                result = unix(['traf -o ' fn]);
            else
                [result, out] = unix(['traf -o ' fn]);
            end
            assert(result == 0);
            f = fopen([fn '.ieq'], 'r');
            [A, b, Aeq, beq, x] = PolytaPortaLibrary.parse_ieq(f);
            fclose(f);
            if nargin < 4 || ~isfield(options, 'delete') || options.delete
                delete(fn);
                delete([fn '.ieq']);
            end
        end
        function [A, b, Aeq, beq, x] = projection(L, A, b, ...
                                                  Aeq, beq, ...
                                                  x, dims, options)
            fn = [tempname '.ieq'];
            f = fopen(fn, 'w');
            N = size(x, 2);
            EO = zeros(1, N);
            EO(dims) = 1:size(dims, 2);
            PolytaPortaLibrary.write_ieq(f, A, b, Aeq, beq, x, EO);
            fclose(f);
            if nargin < 8 || ~isfield(options, 'verbose') || options.verbose
                result = unix(['fmel ' fn]);
            else
                [result, out] = unix(['fmel ' fn]);
            end
            f = fopen([fn '.ieq'], 'r');
            [A, b, Aeq, beq, x] = PolytaPortaLibrary.parse_ieq(f);
            Aeqsup = zeros(size(dims, 2), size(A, 2));
            beqsup = zeros(size(dims, 2), 1);
            Aeqsup(:, dims) = eye(size(dims, 2));
            Aeq = rational([Aeq; Aeqsup]);
            beq = rational([beq; beqsup]);
            fclose(f);
            if nargin < 8 || ~isfield(options, 'delete') || options.delete
                delete(fn);
                delete([fn '.ieq']);
            end
        end
    end
    methods(Static)
        function [E, R] = parse_poi(fid)
        % E = Polyta.parse_poi(fid)
        % Parse a .poi PORTA file
        %
        % INPUTS
        %     fid - opened file descriptor (read mode)
        % OUTPUTS
        %     E   - Extremal points as a matrix of row vectors
        %     R   - Rays as a matrix of row vectors
            DIM = [];
            CONV = {};
            CONE = {};
            section = '';
            while 1
                tl = fgetl(fid);
                if ~ischar(tl)
                    break
                end
                tl = strtrim(tl);
                if ~isequal(tl, '')
                    switch strtok(tl)
                      case 'DIM'
                        DIM = tl;
                      case {'CONV_SECTION', 'CONE_SECTION'}
                        section = strtok(tl);
                      case 'END'
                        break
                      otherwise
                        switch section
                          case 'CONV_SECTION'
                            CONV = vertcat(CONV, {tl});
                          case 'CONE_SECTION'
                            CONE = vertcat(CONE, {tl});
                        end
                    end
                end
            end
            % now parse inside sections
            [s, e, t, m, toks] = regexp(DIM, 'DIM\s+=\s+(\d+)');
            assert(isequal(size(toks), [1 1]));
            % dimension
            D = str2num(toks{1}{1});
            E = parse_section(CONV);
            R = parse_section(CONE);
            if size(E, 1) == 0
                E = [];
            end
            if size(R, 1) == 0
                R = [];
            end
            function M = parse_section(SECT)
            % parse section SECT
                MN = zeros(size(SECT, 1), D);
                MD = ones(size(SECT, 1), D);
                for i = 1:size(SECT, 1)
                    pts = '[+-]?\d+\s*(/\s*\d+)?';
                    [s, e, t, m, toks] = regexp(SECT{i}, ['(\(\s*\d+\))'...
                                        '((\s*' pts ')+)'], 'once');
                    assert(isequal(size(toks), [1 2]));
                    num = toks{1};
                    points = toks{2};
                    [s, e, t, m, vars] = regexp(points, ['(' pts ')']);
                    for j = 1:size(vars, 2)
                        [MN(i,j), MD(i,j)] = PolytaLibrary.parse_ratnum(vars{j}{1});
                    end
                end
                M = rational(MN)./rational(MD);
            end
        end
        function [A, b, Aeq, beq, x] = parse_ieq(fid)
        % [A, b, Aeq, beq, x] = Polyta.parse_ieq(fid)
        % Parse a .ieq PORTA file
        %
        % INPUTS
        %     fid - opened file descriptor (read mode)
        % OUTPUTS
        %     A   - Inequality matrix
        %     b   - Inequality col. vector
        % such that A x <= b is a constraint.
        %     Aeq - Equality matrix
        %     beq - Equality col. vector
        % such that Aeq x = b is a constraint.
        %     x   - Feasible point
        % BUGS
        %     Does not support (for now) lower and upper bounds

            % Parse text file into sections
            VALID = {};
            DIM = [];
            INEQ = {};
            section = '';
            while 1
                tl = fgetl(fid);
                if ~ischar(tl)
                    break
                end
                tl = strtrim(tl);
                if ~isequal(tl, '')
                    switch strtok(tl)
                      case 'DIM'
                        DIM = tl;
                      case {'VALID', 'INEQUALITIES_SECTION'}
                        section = strtok(tl);
                      case 'END'
                        break
                      case {'LOWER_BOUNDS', 'UPPER_BOUNDS', ...
                            'ELIMINATION_ORDER'}
                        section = 'IGNORE';
                      otherwise
                        switch section
                          case 'VALID'
                            VALID = vertcat(VALID, {tl});
                          case 'INEQUALITIES_SECTION'
                            INEQ = vertcat(INEQ, {tl});
                        end
                    end
                end
            end
            
            % Parse DIM section
            [s, e, t, m, toks] = regexp(DIM, 'DIM\s+=\s+(\d+)');
            assert(isequal(size(toks), [1 1]));
            D = str2num(toks{1}{1}); % dimension
            
            % Parse VALID point section
            if isequal(size(VALID), [1 1])
                VALID = VALID{1};
                [s, e, t, m, toks] = regexp(VALID, '(-?\d+\s*(/\s*\d+)?)');
                xN = double(zeros(size(toks)));
                xD = double(ones(size(toks)));
                for i = 1:size(toks, 2)
                    [xN(i), xD(i)] = PolytaLibrary.parse_ratnum(toks{i}{1});
                end
                assert(size(xN, 2) == D); % dim. is consistent
                x = rational(xN)./xD;
            else
                x = [];
            end
            
            % Parse INEQUALITIES_SECTION section
            % we store A and b as integers, and do the division
            % as a last resort (sym variable creation is very slow)
            AieN = double(zeros(size(INEQ, 1), D));
            AieD = double(ones(size(INEQ, 1), D));
            bieN = double(zeros(size(INEQ, 1), 1));
            bieD = double(ones(size(INEQ, 1), 1));
            % typ(i) = 0 if line i is an equality,
            % typ(i) = 1 if line i is a a <= b inequality,
            % typ(i) = -1 if line i is a a >= b inequality.
            typ = double(zeros(size(INEQ, 1), 1));
            for i = 1:size(INEQ, 1)
                % (num) coeffs {<=,==,=>} cte
                ineq = '([+-]\s*(\d+(\s*/\s*\d+)?)?)\s*x(\d+)';
                [s, e, t, m, toks] = regexp(INEQ{i}, ['(\(\s*\d+\))'...
                                    '((\s*' ineq ')+)'...
                                    '\s*(<=|>=|==)\s*(-?\d+\s*(/\s*\d+)?)'], ...
                                            'once');
                assert(isequal(size(toks), [1 4]));
                num = toks{1};
                coeffs = toks{2};
                switch toks{3}
                  case '<='
                    typ(i) = 1;
                  case '>='
                    typ(i) = -1;
                  case '=='
                    typ(i) = 0;
                  otherwise
                    assert(0);
                end
                cte = toks{4};
                [bieN(i), bieD(i)] = PolytaLibrary.parse_ratnum(cte);
                [s, e, t, m, vars] = regexp(coeffs, ['(' ineq ')']);
                for v = 1:size(vars, 2)
                    [s,e,t,m,toks] = regexp(vars{v}, ineq, 'once');
                    coeff = toks{1}{1};
                    j = str2num(toks{1}{2}); % x_j
                    assert (j <= D);
                    [AieN(i,j), AieD(i,j)] = PolytaLibrary.parse_ratnum(coeff);
                end
            end
            AieN(typ == -1, :) = -AieN(typ == -1, :);
            bieN(typ == -1, :) = -bieN(typ == -1, :);
            A = rational(AieN(typ ~= 0, :)) ./ AieD(typ ~= 0, :);
            b = rational(bieN(typ ~= 0)) ./ bieD(typ ~= 0);
            if size(find(typ == 0), 1) > 0
                Aeq = rational(AieN(typ == 0, :)) ./ AieD(typ == 0, :);
                beq = rational(bieN(typ == 0)) ./ bieD(typ == 0);
            else
                Aeq = [];
                beq = [];
            end
            if size(A, 1) == 0
                A = [];
            end
            if size(b, 1) == 0
                b = [];
            end
        end
        function write_ieq(fid, A, b, Aeq, beq, x, EO)
        % PolytaPortaLibrary.write_ieq(fid, A, b, Aeq, beq, x, EO)
        % Write a .ieq PORTA file
        %
        % INPUTS
        %     fid - opened file descriptor (write mode)
        %     A   - Inequality matrix
        %     b   - Inequality col. vector
        % such that A x <= b is a constraint.
        %     Aeq - Equality matrix
        %     beq - Equality col. vector
        % such that Aeq x = b is a constraint.
        %     x   - Feasible point
        %     EO  - (optional) Elimination order in PORTA format
        %           EO must be a row vector.
        %           EO(i) = 0 if variable i must be eliminated
        %           {EO(i) | EO ~= 0} = {1, 2,.. N} describes the
        %           elimination order
        % BUGS
        %     Does not support (for now) lower and upper bounds
            fprintf(fid, 'DIM = %d\n\n', size(x, 2));
            fprintf(fid, 'VALID\n');
            fmt = '%s';
            for j = 1:size(x, 2)
                fprintf(fid, fmt, char(x(j)));
                fmt = ' %s';
            end
            fprintf(fid, '\n\n');
            if nargin == 7
                fprintf(fid, 'ELIMINATION_ORDER\n');
                fmt = '%d';
                for j = 1:size(EO, 2)
                    fprintf(fid, fmt, EO(j));
                    fmt = ' %d';
                end
                fprintf(fid, '\n');
            end
            fprintf(fid, 'INEQUALITIES_SECTION\n');
            Aie = [A;Aeq];
            bie = [b;beq];
            typ = zeros(size(bie));
            typ(1:size(b, 1)) = 1;
            % find widths
            Ac = cell(size(Aie));
            wA = zeros(1, size(Aie, 2));
            matt = num(Aie);
            for i = 1:size(Aie, 1)
                nonzero = find(matt(i,:) ~= 0);
                for j = setdiff(1:size(Aie,2), nonzero)
                    Ac{i,j} = '';
                end
                for j = nonzero
                    sgn = '+';
                    if abs(Aie(i,j)) == Aie(i,j)
                        sgn = '+';
                    else
                        sgn = '';
                    end
                    var = ['x' num2str(j)];
                    Ac{i,j} = [sgn char(Aie(i,j)) var];
                    wA(j) = max(wA(j), size(Ac{i,j},2));
                end
            end
            bc = cell(size(bie));
            wb = 0;
            for i = 1:size(bie, 1)
                bc{i} = char(bie(i));
                wb = max(wb, size(bc{i},2));
            end
            % write
            pad = @(n) char(32*ones(1, n));
            padnum = @(f, w) [pad(w - length(num2str(f))) num2str(f)];
            for i = 1:size(Aie, 1)
                str = sprintf('(% 4d)', i);
                for j = 1:size(Aie, 2)
                    str = [str padnum(Ac{i,j}, wA(j))];
                end
                switch typ(i)
                  case 1
                    str = [str ' <= ' padnum(bc{i}, wb) '\n'];
                  case 0
                    str = [str ' == ' padnum(bc{i}, wb) '\n'];
                end
                fprintf(fid, str);
            end
            fprintf(fid, '\nEND\n');
        end
        function write_poi(fid, E, R)
        % PolytaPortaLibrary.write_poi(fid, E, R)
        % Write a .poi PORTA file
        %
        % INPUTS
        %     fid - opened file descriptor (write mode)
        %     E   - matrix of extremal points (one per row)
        %     R   - matrix of rays (one per row)
        % BUGS
        %     Does not support unbounded polytopes.
            fprintf(fid, 'DIM = %d\n\n', size(E, 2));
            fprintf(fid, 'CONV_SECTION\n');
            write_section(E);
            if ~isequal(R, [])
                fprintf(fid, 'CONE_SECTION\n');
                write_section(R);
            end
            fprintf(fid, '\nEND\n');
            function write_section(M)
                Mc = cell(size(M));
                wM = zeros(1, size(M, 2));
                for j = 1:size(M, 2)
                    for i = 1:size(M, 1)
                        Mc{i,j} = char(M(i, j));
                        wM(j) = max(wM(j), size(Mc{i,j}, 2));
                    end
                end
                for i = 1:size(M, 1)
                    fprintf(fid, '(% 4d)', i);
                    for j = 1:size(M, 2)
                        fprintf(fid, [' % ' num2str(wM(j)) 's'], Mc{i,j});
                    end
                    fprintf(fid, '\n');
                end
            end
        end
    end
end
