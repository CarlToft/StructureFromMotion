
function R=getrotation(alfa,beta,gamma)



R = [cos(alfa) -sin(alfa) 0;sin(alfa) cos(alfa) 0;0 0 1];
R = R*[cos(beta) 0 -sin(beta);0 1 0;sin(beta) 0 cos(beta)];
R = R*[1 0 0;0 cos(gamma) sin(gamma);0 -sin(gamma) cos(gamma)];







