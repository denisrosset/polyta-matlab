classdef PolytaLibraryManager
    properties
        supported_libs = {...
            PolytaPortaLibrary() ...
            PolytaCddLibrary() ...
            PolytaLrsLibrary() ...
            PolytaPdLibrary() ...
            PolytaPplLibrary() ...
                         };
        installed_libs = {};
        activated_libs = {};
        options = struct('verbose', true, 'delete', true);
    end
    methods
        function PL = PolytaLibraryManager()
            for lib = PL.supported_libs
                if lib{:}.isAvailable()
                    PL.installed_libs = horzcat(PL.installed_libs, lib);
                end
            end
            PL.activated_libs = PL.installed_libs;
        end
    end
    methods(Static)
        function obj = instance()
            persistent uniqueInstance;
            if isempty(uniqueInstance)
                uniqueInstance = PolytaLibraryManager();
            end
            obj = uniqueInstance;
        end
    end
end
