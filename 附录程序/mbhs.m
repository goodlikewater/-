% Ä¿±êº¯Êý
function y=mbhs(g,f)
 
y=0;
for n=1:length(g)
    y=y+(g(n)-f(n))^2;
end
end
