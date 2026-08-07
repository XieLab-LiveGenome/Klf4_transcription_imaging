function check_python_env()
    pe = pyenv;
    if pe.Status == "NotLoaded" && pe.Executable == ""
        error(['No Python interpreter configured. Run: ' ...
               'pyenv(''Version'', ''/path/to/env/bin/python'')']);
    end
    try
        v = char(py.importlib.metadata.version('cellpose'));
    catch
        error('cellpose not importable from %s. See README.', pe.Executable);
    end
    if sscanf(v, '%d', 1) >= 4
        warning('Cellpose %s found; pipelines were validated on 3.x.', v);
    end
    fprintf('Python: %s | cellpose %s\n', pe.Executable, v);
end