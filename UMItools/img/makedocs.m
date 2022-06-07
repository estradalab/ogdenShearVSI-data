f = dir('*.m')
diary ('help.txt')
for n=1:length(f)
    str = sprintf('help %s',f(n).name);
    eval(str)
end
diary off