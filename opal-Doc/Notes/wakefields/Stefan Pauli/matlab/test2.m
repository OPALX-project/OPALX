function test2()
% just used for tests

load ('wakefieldMathematica.mat')

compareData('wakefield from matlab','wakefield from mathematica','s[um]', [V/pC/m], wake,wakefieldMathematica, xStart, xEnd)


end