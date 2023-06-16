function  [funcs] = createEdgeFuncs(edge,l,w)

param = [l,w,-l,-w];
param_y = [w,l,-w,-l];
for i = 1:length(edge.shape)
    switch edge.shape{i}
        case 'sin'
            B = 2*pi*edge.period/abs(param(i));
            funcs{i} = @(x) edge.coef{i}*sin(B*x) + param_y(i)/2;
        case 'parabol'
            coef = 1;
            funcs{i} = @(x) edge.coef{i}*(((2*x/abs(param(i))).^2)-1) + param_y(i)/2;
        case 'cos'
            B = pi/abs(param(i));
            coef = 1;
            funcs{i} = @(x) edge.coef{i}*cos(B*x) + param_y(i)/2;
        case 'line'
            funcs{i} = @(x) 0*x + param_y(i)/2;
    end
end