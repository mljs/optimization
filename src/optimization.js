var math = require("mathjs");
require("util");

Optimization = {

    /*
    * Function that transform a "Matrix" vector into a Vector
    *
    * @param {Matrix} vect - Vector that mathjs suppose as a matrix.
    * @param {Vector} normalVector - Vector transform.
    *
    * */
    _toNormalVector : function(vect) {
        var normalVector = [];
        var size = math.subset(math.size(vect), math.index(1));
        for (var i = 0; i < size; ++i) {
            normalVector.push(math.subset(vect, math.index(0, i)));
        }
        return normalVector;
    },

    /*
    * Function that obtains the local minimum vector of a function with a vertices of a Simplex.
    *
    * @param {Function} f - Function to get the local minimum, this function require one argument and should be a vector.
    * @param {Vector} vertices - Vertices of the Simplex shape.
    * @param {Number} minimumIterations - Minimum number of iterations we want.
    * @param {Number} maximumIterations - Maximum number of iterations we want.
    * @param {Number|BigNumber} epsilon - Local minimum tolerance.
    * @return {Object} result - return an object that contains a minimum vector of the function and the evaluation of the function in this vector.
    *
    * */
    nelderMead : function (f, vertices, minimumIterations, maximumIterations, epsilon) {
        var size = math.subset(math.size(vertices), math.index(1));
        var y = [];
        for(var i = 0; i < size + 1; ++i) {
            var vector = vertices[i];
            y.push(f(vector));
        }

        var low = Util.getIndexRow(y, math.min(y));
        var high = Util.getIndexRow(y, math.max(y));
        var li = high;
        var ho = low;

        for(var i = 0; i <= size; ++i) {
            if(!math.equal(i, low) && !math.equal(i, high) && math.smallerEq(y[i], y[li])) {
                li = i;
            }
            if(!math.equal(i, low) && !math.equal(i, high) && math.largerEq(y[i], y[ho])) {
                ho = i;
            }
        }

        var counter = 0;

        while( (math.larger(y[high], math.add(y[low], epsilon)) && counter < maximumIterations) || counter < minimumIterations ) {

            var S = math.zeros(1, size);
            for(var i = 0; i < size + 1; ++i) {
                S = math.add(S, math.subset(vertices, math.index(i, [0, size])));
            }

            var hiVector = math.subset(vertices, math.index(high, [0, size]));
            var M = math.divide(math.subtract(S, hiVector), size);
            var R = math.subtract(math.multiply(2, M), hiVector);
            var yR = f(this._toNormalVector(R));

            if(math.smaller(yR, y[ho])) {
                if(math.smaller(y[li], yR)) {
                    vertices = math.subset(vertices, math.index(high, [0, size]), R);
                    y[high] = yR;
                } else {
                    var E = math.subtract(math.multiply(2, R), M);
                    yE = f(this._toNormalVector(E));

                    if(math.smaller(yE, y[li])) {
                        vertices = math.subset(vertices, math.index(high, [0, size]), E);
                        y[high] = yE;
                    } else {
                        vertices = math.subset(vertices, math.index(high, [0, size]), R);
                        y[high] = yR;
                    }
                }
            } else {
                if(math.smaller(yR, y[high])) {
                    vertices = math.subset(vertices, math.index(high, [0, size]), R);
                    y[high] = yR;
                }

                var actual = math.subset(vertices, math.index(high, [0, size]));
                var C = math.divide(math.add(actual, M), 2);
                var yC = f(this._toNormalVector(C));

                var C2 = math.divide(math.add(M, R), 2);
                var yC2 = f(this._toNormalVector(C2));

                if(math.smaller(yC2, yC)) {
                    C = C2;
                    yC = yC2;
                }
                if(math.smaller(yC, y[high])) {
                    vertices = math.subset(vertices, math.index(high, [0, size]), C);
                    y[high] = yC;
                } else {
                    for(var i = 0; i <= size; ++i) {
                        if(i != low) {
                            var vectorj = math.subset(vertices, math.index(j, [0, size]));
                            var vectorlo = math.subset(vertices, math.index(low, [0, size]));

                            vertices = math.subset(vertices, math.index(j, [0, size]), math.divide(math.add(vectorj, vectorlo),2));

                            var Z = math.subset(vertices, math.index(j, [0, size]));
                            y[i] = f(this._toNormalVector(Z));
                        }
                    }
                }
            }

            low = Util.getIndexRow(y, math.min(y));
            high = Util.getIndexRow(y, math.max(y));
            li = high;
            ho = low;

            for(var i = 0; i <= size; ++i) {
                if(!math.equal(i, low) && !math.equal(i, high) && math.smallerEq(y[i], y[li])) {
                    li = i;
                }
                if(!math.equal(i, low) && !math.equal(i, high) && math.largerEq(y[i], y[ho])) {
                    ho = i;
                }
            }
            counter++
        }

        var V0 = math.subset(vertices, math.index(low, [0, size]));
        var y0 = f(this._toNormalVector(V0));
        return {vector: this._toNormalVector(V0), eval: y0};
    }
};
