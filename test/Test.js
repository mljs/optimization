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

    it('Gradiend descent', function () {
        var p0 = [3.5, 2.5];
        var f = math.eval("f(x) = (subset(x, index(1)) * subset(x, index(1))) - 4 * subset(x, index(1)) + " +
        "(subset(x, index(2)) * subset(x, index(2))) - subset(x, index(2)) - subset(x, index(1)) * subset(x, index(2))");
        var g = function(x) {
            var z = math.zeros(1, 2);
            var g = math.eval("f(x) = matrix([subtract(multiply(2, subset(x, index(1))), subtract(4, subset(x, index(2))))"  +
                                            ", subtract(multiply(2, subset(x, index(2))), subtract(1, subset(x, index(1))))])");
            var valG = g(x);
            var result = math.multiply(math.divide(1, math.norm(valG)), valG);
            return math.multiply(-1, result);
        };

        var result = Optimization.gradientDescent(f, g, p0, 10000, 10e-5, 10e-5);

        var optimumVector = [3.1126511, 1.99983839];

        for(var i = 0; i < optimumVector.length; ++i) {
            (math.subset(result.vector, math.index(i))).should.be.approximately(optimumVector[i], 10e-2);
        }

        (result.eval).should.be.approximately(-6.99999998, 10e-2);

    });
});