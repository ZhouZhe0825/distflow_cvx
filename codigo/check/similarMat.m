function [disc] = similarMat(A, B, tol)
% si disc = [], las matrices son similares, abs(A-B) < tol para todos los
% elementos
% si disc = [i j v], los elementos (i,j) difieren en v
% si disc = [inf inf inf] las matrices son de distintas dimensiones


disc = [];
Af = full(A);
Bf = full(B);

if isequal(size(Af), size(Bf))
    C = abs(Af-Bf);
    [i, j] = find(C > tol);
    if size(i,1) ~= length(i)
        i = i';
        j = j';
    end
    if ~isempty(i)
        v = C(sub2ind(size(C),i,j));
        if size(v,1) ~= length(v)
            v = v';
        end
        disc = [i j v];
    end
else
    disc = [inf inf inf];
end