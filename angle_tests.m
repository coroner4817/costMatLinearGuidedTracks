% Simple tests for a few basic angles

% test precision
prec = 10^-5;

a = [0 1];
b = [1 1];
assert(abs(vvAngle(a, b) - 45) < prec);

a = [0 1];
b = [1 0];
assert(abs(vvAngle(a, b) - 90) < prec);

a = [0 1];
b = [1 -1];
assert(abs(vvAngle(a, b) - 135) < prec);

a = [1 1];
b = [-1 -1];
assert(abs(vvAngle(a, b) - 180) < prec);

a = [1 1];
b = [1 1];
assert(abs(vvAngle(a, b) - 0) < prec);

a = [1 1];
b = [0 0];
assert(isnan(vvAngle(a, b)));

a = [0 0];
b = [0 0];
assert(isnan(vvAngle(a, b)));

fprintf(1, 'all good.\n');