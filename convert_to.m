function A = convert_to(A, cl)
    if isequal(class(A), cl)
        return
    end
    switch cl
      case 'rational'
        A = rational(A);
      case 'sym'
        A = sym(A);
      case 'double'
        A = double(A);
      case 'single'
        A = single(A);
      otherwise
        error(['Unknown numeric class ' cl]);
    end
end
