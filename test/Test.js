require("../src/optimization");
require("../src/util");
var math = require('mathjs');

describe('Optimization', function () {

    // Test case obtained from Pag 443, Chap 8.
    it('Nelder-Mead', function () {
        var vertices = [[0, 0], [1.2, 0.0], [0.0, 0.8]];
        var f = math.eval("f(x) = (subset(x, index(1)) * subset(x, index(1))) - 4 * subset(x, index(1)) + " +
                          "(subset(x, index(2)) * subset(x, index(2))) - subset(x, index(2)) - subset(x, index(1)) * subset(x, index(2))");
        var result = Optimization.nelderMead(f, vertices, 33, 35, 10e-4);

        var optimumVector = [2.99996456, 1.99983839];

        for(var i = 0; i < vertices[0].length; ++i) {
            (math.subset(result.vector, math.index(i))).should.be.approximately(optimumVector[i], 10e-4);
        }

        (result.eval).should.be.approximately(-6.99999998, 10e-5);
    });
});