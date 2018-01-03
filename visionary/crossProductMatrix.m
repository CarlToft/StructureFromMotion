% A = crossProductMatrix(a)
% 
% Takes as input a vector a in R^3 and outputs the corresponding cross
% product matrix, i.e., the matrix A such that A*u = a x u for all vectors
% u. 
% 
function A = crossProductMatrix(a)
    A = [0, -a(3), a(2); a(3), 0, -a(1); -a(2), a(1), 0];    
end
