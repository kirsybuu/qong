/*
Copyright (c) 2017 Michael Joseph Coulombe
*/

let assert = function(b) {
    if (b == false) {
        console.log("Assertion Failure!");
        //var stop = false;
        //while(!stop) {}
    }
}

let percent = function(p) {
    return Math.floor(p * 10000) / 100;
};

let bitflip = function(strbit) {
    if (strbit == "0") return "1";
    return "0";
}

let intBits = function(v, numBits) {
    var a = [];
    for (var i = 0; i < numBits; i++) {
        if (v % 2 == 1) a.push("1");
        else            a.push("0");
        v = Math.floor(v/2);
    }
    return a.join("");
};

let intFromBits = function(bits) {
    var v = 0;
    for (var i = bits.length - 1; i >= 0; i--) {
        v *= 2;
        if (bits[i] == "1") {
            v += 1;
        }
    }
    return v;
};

let feq = function(a,b) {
    return Math.floor(a * 1024) == Math.floor(b * 1024);
}

let charsAt = function(str, indices) {
    var restricted = "";
    for (var i of indices) {
        restricted = restricted.concat(str.charAt(i));
    }
    return restricted
};

// complex numbers
let complex = function(real,im) {
    return {
        r: real,
        i: im,
        add: function(z) {
            return complex(this.r + z.r, this.i + z.i);
        },
        sub: function(z) {
            return complex(this.r - z.r, this.i - z.i);
        },
        mult: function(z) {
            return complex(this.r * z.r - this.i * z.i, this.r * z.i + this.i * z.r)
        },
        div: function(z) {
            let denom = z.r * z.r + z.i * z.i;
            return complex(
                (this.r * z.r + this.i * z.i) / denom,
                (this.i * z.r - this.r * z.i) / denom
            );
        },
        mag: function() {
            return this.r * this.r + this.i * this.i;
        },
        equal: function(z) {
            return this.r == z.r && this.i == z.i;
        },
        copy: function() {
            return complex(this.r, this.i);
        }
    };
};
let cZero = function() { return complex(0,0); }
let cOne = function() { return complex(1,0); }
let cIm = function() { return complex(0,1); }

// qubits and extanged groups
let qubit = function(amp0, amp1) {
    return {
        z0: amp0,
        z1: amp1
    };
};

let isEmpty = function(map) {
    for (var k in map) return false;
    return true;
};

let totalAmplitude = function(terms) {
    let amp = 0;
    for (var tensor in terms) {
        amp += terms[tensor].mag();
    }
    return amp;
};

let normalize = function(terms) {
    if (isEmpty(terms)) return;
    while (true) {
        let amp = totalAmplitude(terms);
        if (feq(amp,1)) return;
        if (amp == 0) return;
        //console.log("normalize: %f", amp);
        let cAmp = complex(1 / Math.sqrt(amp),0);
        for (var tensor in terms) {
            terms[tensor] = terms[tensor].mult(cAmp);
        }
    }
};

let qstate = function() {
    var init = {
        n: 0,
        terms: {},
        cond: [],
        // set qstate to an initial non-entangled state
        resetTo: function(tensor, phase) {
            this.n = tensor.length;
            this.terms = {};
            this.terms[tensor] = phase;
            this.cond.length = this.n;
            this.cond.fill(false);
        },
        /*
        // insert a new qubit at position i
        insertQubit: function(q, i) {
            var newTerms = {};
            for (var tensor in terms) {
                let before = tensor.slice(0, i);
                let after = tensor.slice(i);
                if (q.z0 > 0) {
                    newTerms[before.concat("0",after)] = q.z0.mult(terms[tensor]);
                }
                if (q.z1 > 0) {
                    newTerms[before.concat("1",after)] = q.z1.mult(terms[tensor]);
                }
            }
            this.n = this.n + 1;
            this.terms = newTerms;
            this.cond.splice(i, 0, false);
        },
        */
        // mark qubits as being tested in a conditional gate
        requireOne: function(indices) {
            for (var i of indices) {
                this.cond[i] = true;
            }
        },
        // unmark qubits as being tested in a conditional gate
        unrequireOne: function(indices) {
            for (var i of indices) {
                this.cond[i] = false;
            }
        },
        passesCond: function(tensor) {
            for (var i = 0; i < this.n; i++) {
                if (this.cond[i] && tensor.charAt(i) == "0") return false;
            }
            return true;
        },
        // apply f to each non-entangled state
        map: function(f) {
            var newTerms = {};
            for (var tensor in this.terms) {
                let amp = this.terms[tensor];
                if (this.passesCond(tensor)) {
                    // apply f and merge its scaled results
                    let result = f(tensor);
                    for (let newTensor in result) {
                        if (! (newTensor in newTerms)) {
                            newTerms[newTensor] = cZero();
                        }
                        newTerms[newTensor] = newTerms[newTensor].add( result[newTensor].mult(amp) );
                    }
                }
                else {
                    // don't apply f
                    newTerms[tensor] = amp;
                }
            }
            this.terms = newTerms;
        },
        flip: function(indices) {
            for (var i of indices) {
                assert(this.cond[i] == false);
            }
            this.map(function(tensor) {
                var arr = [];
                var last = 0;
                for (var i of indices) {
                    arr.push(tensor.slice(last, i));
                    arr.push(bitflip(tensor.charAt(i)));
                    last = i+1;
                }
                arr.push(tensor.slice(last));
                let flipped = arr.join("");
                /*
                let before = tensor.slice(0, i);
                let after = tensor.slice(i+1);
                let flipped = before.concat(bitflip(tensor.charAt(i)), after);
                */
                var s = {};
                s[flipped] = cOne();
                return s;
            });
        },
        hadamard: function(i) {
            assert(this.cond[i] == false);
            let posAmp = complex(1 / Math.sqrt(2), 0);
            let negAmp = posAmp.mult(complex(-1,0));
            
            this.map(function(tensor) {
                let before = tensor.slice(0, i);
                let val = tensor.charAt(i);
                let after = tensor.slice(i+1);
                let flipped = before.concat(bitflip(val), after);
                
                var s = {};
                if (val == "0") {
                    s[tensor] = posAmp;  // |0>
                    s[flipped] = posAmp; // |1>
                }
                else {
                    s[tensor] = negAmp;  // |1>
                    s[flipped] = posAmp; // |0>
                }
                return s;
            });
        },
        ygate: function(i) {
            assert(this.cond[i] == false);
            let zeroAmp = complex(1 / Math.sqrt(2), 0); 
            let posAmp = complex(0, zeroAmp.r);
            let negAmp = posAmp.mult(complex(-1,0));
            
            this.map(function(tensor) {
                let before = tensor.slice(0, i);
                let val = tensor.charAt(i);
                let after = tensor.slice(i+1);
                let flipped = before.concat(bitflip(val), after);
                
                var s = {};
                if (val == "0") {
                    s[tensor] = zeroAmp;  // |0>
                    s[flipped] = posAmp; // |1>
                }
                else {
                    s[tensor] = negAmp;  // |1>
                    s[flipped] = zeroAmp; // |0>
                }
                return s;
            });
        },
        /*
        swap: function(i,j) {
            assert(i < j);
            assert(this.cond[i] == false && this.cond[j] == false);
            this.map(function(tensor) {
                let vi = tensor.charAt(i);
                let vj = tensor.charAt(j);
                if (vi == vj) {
                    var s = {};
                    s[tensor] = cOne();
                    return s;
                }
                else {
                    let before  = tensor.slice(0,   i);
                    let between = tensor.slice(i+1, j);
                    let after   = tensor.slice(j+1);
                    let swapped = before.concat(vj, between, vi, after);
                    var s = {};
                    s[swapped] = cOne();
                    return s;
                }
            });
        },
        */
        forceSet: function(i, val) {
            this.map(function(tensor) {
                var s = {};
                let vi = tensor.charAt(i);
                if (vi == val) {
                    s[tensor] = cOne();
                }
                else {
                    let before = tensor.slice(0, i);
                    let after = tensor.slice(i+1);
                    let flipped = before.concat(val, after);
                    s[flipped] = cOne();
                }
                return s;
            });
        },
        probOf: function(f) {
            var prob = 0;
            for (var tensor in this.terms) {
                if (f(tensor)) {
                    prob += this.terms[tensor].mag();
                }
            }
            return prob;
        },
        /*
        distOf: function(indices) {
            var dist = {};
            for (var tensor in this.terms) {
                let amp = this.terms[tensor];
                // apply f and merge its scaled results
                var restricted = charsAt(tensor, indices);
                if (! (restricted in dist)) {
                    dist[restricted] = amp.mag();
                }
                else {
                    dist[restricted] = dist[restricted] + amp.mag();
                }
            }
            return dist;
        },
        */
        forEach: function(f) {
            for (var tensor in this.terms) {
                f(tensor, this.terms[tensor].copy());
            }
        },
        postSelect: function(f) {
            var missingAmp = 0;
            var newTerms = {};
            for (var tensor in this.terms) {
                let amp = this.terms[tensor];
                if (f(tensor)) {
                    // keep tensor
                    newTerms[tensor] = amp;
                }
                else {
                    missingAmp += amp.mag();
                    // remove tensor
                }
            }
            if (newTerms.length == 0) {
                missingAmp = 1;
            }
            if (missingAmp > 0) {
                //console.log("post-select elim %f", missingAmp);
            }
            /*let scale = complex(1.0 / Math.sqrt(1 - missingAmp), 0);
            for (var tensor in newTerms) {
                newTerms[tensor] = newTerms[tensor].mult(scale);
            }
            */
            normalize(newTerms);
            this.terms = newTerms;
            return 1 - missingAmp;
        },
        normalize: function() {
            normalize(this.terms);
        }
    };
    return init;
};

let circuitLogger = function(n, logFunc) {
    let blank = " ";
    let condChar = "c";
    let flipChar = "f";
    let hChar = "h";
    let yChar = "y";
    
    let initLog = new Array(n);
    
    for (var i = 0; i < n; i++) {
        initLog[i] = " ";
    }
    
    return {
        log: initLog,
        anyLoggedGates: function() {
            for (var i = 0; i < n; i++) {
                if (this.log[i] != blank && this.log[i] != condChar) {
                    return true;
                }
            }
            return false;
        },
        clearLog: function() {
            var removedAny = false;
            for (var i = 0; i < n; i++) {
                if (this.log[i] != condChar) {
                    this.log[i] = blank;
                    removedAny = true;
                }
            }
            return removedAny;
        },
        outputLog: function() {
            // output
            logFunc(this.log);
            // clear any gates that were only applied previously
            this.clearLog();
        },
        // mark qubits as being tested in a conditional gate
        requireOne: function(indices) {
            if (this.anyLoggedGates()) {
                this.outputLog();
            }
            for (var i of indices) {
                this.log[i] = condChar;
            }
        },
        // unmark qubits as being tested in a conditional gate
        unrequireOne: function(indices) {
            if (this.anyLoggedGates()) {
                this.outputLog();
            }
            for (var i of indices) {
                assert(this.log[i] == condChar);
                this.log[i] = blank;
            }
        },
        pushGates: function(indices, ch) {
            // ensure that each index is blank and ready for a gate
            for (var i of indices) {
                assert(this.log[i] != condChar);
                if (this.log[i] != blank) {
                    this.outputLog();
                    break;
                }
            }
            // add the blanks
            for (var i of indices) {
                assert(this.log[i] == blank);
                this.log[i] = ch;
            }
        },
        flip: function(indices) {
            this.pushGates(indices, flipChar);
        },
        hadamard: function(i) {
            this.pushGates([i], hChar);
        },
        ygate: function(i) {
            this.pushGates([i], yChar);
        },
    };
};

/////////////////////////////// quantum interface //////////////////////////////

// qubit gates - transfer values from input to output wires

let qflip = function(indices) {
    assert(indices instanceof Array);
    return function(qs) {
        qs.flip(indices);
    };
};

let qseq = function(gates) {
    assert(gates instanceof Array);
    return function(qs) {
        for (var g of gates) {
            g(qs);
        }
    };
};

let qctrl = function(indices, gate) {
    assert(indices instanceof Array);
    return function(qs) {
        //for (var i of indices) {
        //    qs.requireOne(i);
        //}
        qs.requireOne(indices);
        gate(qs);
        //for (var i of indices) {
        //    qs.unrequireOne(i);
        //}
        qs.unrequireOne(indices);
    };
};

let qHadamard = function(indices) {
    assert(indices instanceof Array);
    return function(qs) {
        for (var i of indices) {
            qs.hadamard(i);
        }
    };
};

let qYGate = function(indices) {
    assert(indices instanceof Array);
    return function(qs) {
        for (var i of indices) {
            qs.ygate(i);
        }
    };
};

////////////////////////////////// qong circuit ////////////////////////////////

var Qong = {};

Qong.xorBits = function(x, val, bit) {
    let n = x.length;
    let valBits = intBits(val, n);
    var flipGates = [];
    for (var i = 0; i < n; i++) {
        if (valBits[i] == bit) {
            flipGates.push(qflip(x.slice(i,i+1)));
        }
    }
    return qseq(flipGates);
};

// increment
Qong.inc = function(indices) {
    let n = indices.length;
    if (n == 1) {
        return qflip(indices);
    }
    let lowIndices = indices.slice(0, n-1);
    let highIndex = indices.slice(n-1);
    return qseq([
        qctrl(lowIndices, qflip(highIndex)),
        Qong.inc(lowIndices)
    ]);
};
// decrement = inverse of increment
Qong.dec = function(indices) {
    let n = indices.length;
    if (n == 1) {
        return qflip(indices);
    }
    let lowIndices = indices.slice(0, n-1);
    let highIndex = indices.slice(n-1);
    return qseq([
        Qong.dec(lowIndices),
        qctrl(lowIndices, qflip(highIndex)),
    ]);
};
Qong.add = function(indices, val) {
    var gates = [];
    for (var i = 0; i < val; i++) {
        gates.push(Qong.inc(indices));
    }
    return qseq(gates);
};
Qong.sub = function(indices, val) {
    var gates = [];
    for (var i = 0; i < val; i++) {
        gates.push(Qong.dec(indices));
    }
    return qseq(gates);
};

// x++ if r else x--
Qong.ts = function(x,r) {
    return qseq([
        qctrl(r, Qong.inc(x)),
        qflip(r),
        qctrl(r, Qong.dec(x)),
        qflip(r)
    ]);
};

// flip d if y is at the top or bottom
Qong.bnc = function(gridHeight, y, d) {
    let botMask0 = Qong.xorBits(y, gridHeight - 1, "0");
    let botMask1 = Qong.xorBits(y, gridHeight - 1, "1");
    
    return qseq([
        qflip(y),
        qctrl(y, qflip(d)),
        botMask1,
        qctrl(y, qflip(d)),
        botMask0
    ]);
    /*
    return qseq([
        qctrl(y, qflip(d)),
        qflip(y),
        qctrl(y, qflip(d)),
        qflip(y)
    ]);
    */
};

// move up/down but bounce off of walls
Qong.vstep = function(gridHeight, y, d) {
    return qseq([
        Qong.bnc(gridHeight, y, d),
        Qong.ts(y,d)
    ]);
};

// (x, y) -> (x=y, y)
Qong.eq2 = function(x,y) {
    return qseq([
        qctrl(y, qflip(x)),
        qflip(x),
    ]);
};
Qong.eq2_inv = function(x,y) {
    return qseq([
        qflip(x),
        qctrl(y, qflip(x))
    ]);
};

// flip b if x < p 
Qong.less = function(x,p,b) {
    let n = x.length;
    if (n == 1) {
        return qseq([
            qflip(x),
            qctrl(x.concat(p), qflip(b)),
            qflip(x)
        ]);
    }
    else {
        let xmsb = x.slice(n-1);
        let pmsb = p.slice(n-1);
        return qseq([
            Qong.eq2(xmsb, pmsb),
            qctrl(xmsb, Qong.less(x.slice(0,n-1), p.slice(0,n-1), b)),
            Qong.eq2_inv(xmsb, pmsb),
            Qong.less(xmsb, pmsb, b)
        ]);
    }
};
// flip r if next to the paddle about to reflect
Qong.hbnc = function(paddleSide, paddleLength, x, y, paddleMin, r) {
    /*
    let n = x.length;
    var xFlipGates = [];
    for (var i = 0; i < n; i++) {
        if (paddleSideX[i] == "0") {
            xFlipGates.push(qflip(x.slice(i,i+1)));
        }
    }
    let xMask = qseq(xFlipGates);
    */
    let xMask = Qong.xorBits(x, paddleSide, "0");
    return qseq([
        xMask,
        qctrl(x, Qong.less(y, paddleMin, r)),
        Qong.add(paddleMin, paddleLength-1),
        qctrl(x, Qong.less(paddleMin, y, r)),
        Qong.sub(paddleMin, paddleLength-1),
        qctrl(x, qflip(r)),
        xMask
    ]);
};
// horizontal motion, bouncing off paddles
Qong.hstep = function(gridWidth, paddleDist, paddleLength, x, y, paddleMinLeft, paddleMinRight, r) {
    let n = x.length;
    //let N = Math.pow(2,n);
    
    let paddleSideLeft = paddleDist + 1;
    let paddleSideRight = gridWidth-1 - paddleDist - 1;
    
    return qseq([
        Qong.hbnc(paddleSideLeft, paddleLength, x, y, paddleMinLeft, r),
        Qong.hbnc(paddleSideRight, paddleLength, x, y, paddleMinRight, r),
        Qong.ts(x, r)
    ]);
};

// move paddle in response to buttons, bounce off walls
// TODO: verify correct
/*
goal:
    p,d,u    p,d,b
    0,0,0 => 0,0,0
    0,0,1 => 0,1,0
    0,1,0 => 1,1,1
    
    1,0,0 => 1,0,0
    1,0,1 => 0,0,1
    1,1,0 => 2,1,1
    
    2,0,0 => 2,0,0
    2,0,1 => 1,0,1
    2,1,0 => 2,1,0
    
current:
    // p=pos, d=down, u=up
    if (down) {
        flip(up);
        if (pos == 2) {
            flip(up);
        }
    }
    if (pos == 0) {
        swap(up,down);
    }
    if (up) {
        ts(pos,down);
    }
*/
Qong.swap2 = function(x,y) {
    return qseq([
        qctrl(x, qflip(y)),
        qctrl(y, qflip(x)),
        qctrl(x, qflip(y)),
    ]);
};
Qong.joypad = function(gridHeight, paddleLength, paddleMin, down, up) {
    /*
    let m = paddleMin.length;
    let M = Math.pow(2,m);
    
    let maxBits = intBits(M - paddleLength, m);
    var maxFlip0 = [];
    var maxFlip1 = [];
    for (var i = 0; i < m; i++) {
        if (maxBits[i] == "0") {
            maxFlip0.push(qflip(paddleMin.slice(i,i+1)));
        }
        else {
            maxFlip1.push(qflip(paddleMin.slice(i,i+1)));
        }
    }
    */
    let maxFlip0 = Qong.xorBits(paddleMin, gridHeight - paddleLength, "0");
    let maxFlip1 = Qong.xorBits(paddleMin, gridHeight - paddleLength, "1");
    
    return qseq([
        maxFlip0,
        qctrl(down, qseq([
            qflip(up),
            qctrl(paddleMin, qflip(up))
        ])),
        maxFlip1,
        qctrl(paddleMin, Qong.swap2(down, up)),
        qflip(paddleMin),
        qctrl(up, Qong.ts(paddleMin, down))
    ]);
    /*
    // WRONG
    return qseq([
        qflip(paddleMin),
        qflip(buttonDir),
        qctrl(paddleMin.concat(buttonDir), qflip(buttonPressed)),
        qflip(buttonDir),
        qseq(maxFlip1),
        qctrl(paddleMin.concat(buttonPressed), qflip(buttonDir)),
        qseq(maxFlip0),
        qctrl(buttonPressed, Qong.ts(paddleMin, buttonDir))
    ]);
    */
};

Qong.obs = function(u, v, x, y, d, gate) {
    /*
    let n = x.length;
    let m = y.length;
    let uBits = intBits(u, n);
    let vBits = intBits(v, m);
    var uGates = [];
    var vGates = [];
    for (var i = 0; i < n; i++) {
        if (uBits[i] == "0") {
            uGates.push(qflip(x.slice(i,i+1)));
        }
    }
    for (var i = 0; i < m; i++) {
        if (vBits[i] == "0") {
            vGates.push(qflip(y.slice(i,i+1)));
        }
    }
    */
    let uGates = Qong.xorBits(x, u, "0");
    let vGates = Qong.xorBits(y, v, "0");
    
    return qseq([
        uGates,
        vGates,
        qctrl(x.concat(y), (gate == "H") ? qHadamard(d) : qYGate(d)),
        uGates,
        vGates,
    ]);
};

Qong.allObs = function(obsCoords, x, y, r, d) {
    var gates = [];
    for (var desc of obsCoords) {
        gates.push(Qong.obs(desc.x, desc.y,
                            x, y, (desc.dir == "v") ? d : r,
                            desc.gate));
    }
    return qseq(gates);
};

Qong.fullStep = function(gridHeight, gridWidth, paddleDist, paddleLength, paddleMinLeft, paddleMinRight,
                         x, y, r, d, buttonBits,
                         obsCoords) {
    
    return qseq([
        Qong.allObs(obsCoords, x, y, r, d),
        Qong.hstep(gridWidth, paddleDist, paddleLength, x, y, paddleMinLeft, paddleMinRight, r),
        Qong.vstep(gridHeight, y, d),
        Qong.joypad(gridHeight, paddleLength, paddleMinLeft,  buttonBits.slice(0,1), buttonBits.slice(1,2)),
        Qong.joypad(gridHeight, paddleLength, paddleMinRight, buttonBits.slice(2,3), buttonBits.slice(3,4))
    ]);
};

Qong.paddleAI = function(paddleLength, y, paddleMin, down, up) {
    let paddleHalfway = Math.floor(paddleLength / 2);
    
    return qseq([
        Qong.add(paddleMin, paddleHalfway),
        Qong.less(y, paddleMin, up),
        Qong.less(paddleMin, y, down),
        Qong.sub(paddleMin, paddleHalfway)
    ]);
};

Qong.fullStepAI = function(gridHeight, gridWidth, paddleDist, paddleLength, paddleMinLeft, paddleMinRight,
                         x, y, r, d, buttonBits,
                         obsCoords) {
    
    return qseq([
        Qong.allObs(obsCoords, x, y, r, d),
        Qong.hstep(gridWidth, paddleDist, paddleLength, x, y, paddleMinLeft, paddleMinRight, r),
        Qong.vstep(gridHeight, y, d),
        Qong.paddleAI(paddleLength, y, paddleMinLeft, buttonBits.slice(0,1), buttonBits.slice(1,2)),
        Qong.joypad(gridHeight, paddleLength, paddleMinLeft,  buttonBits.slice(0,1), buttonBits.slice(1,2)),
        Qong.joypad(gridHeight, paddleLength, paddleMinRight, buttonBits.slice(2,3), buttonBits.slice(3,4))
    ]);
};

////////////////////////////////////////////////////////////////////////////////

var ENGINE = {};

var numPlayers = 2;

// main menu state
ENGINE.Menu = {
    smoothing: false,
    
    render: function(dt) {
        this.app.layer.clear("#fff");
        
        this.app.layer.fillStyle("#000")
                      .fillRect(0,0, this.app.width/2, this.app.height);
        
        this.app.layer.font("64px Verdana")
                      .textAlign("left")
                      .fillText("ng", this.app.width/2 - 2, 64)
                      .fillStyle("#fff")
                      .textAlign("right")
                      .fillText("Qo", this.app.width/2 - 2, 64);
                      //.fillText("Players = " + numPlayers + "(press 1 or 2)", 32, 128)
                      //.fillText("Click to start!", 32, 196);
        
        this.app.layer.fillStyle("#fff")
                      .textAlign("left")
                      .fillText("Click", 32, 128+32)
                      .fillText("for 2P", 32, 196+32);
        
        this.app.layer.fillStyle("#000")
                      .textAlign("right")
                      .fillText("Click", this.app.width - 32, 128+32)
                      .fillText("for 1P", this.app.width - 32, 196+32);
        
        this.app.layer.textAlign("left");
    },
    /*
    keydown: function(event) {
        if (event.key == "1") {
            numPlayers = 1;
        }
        if (event.key == "2") {
            numPlayers = 2;
        }
    },
    */
    pointerdown: function(event) {
        if (event.x < this.app.layer.width / 2) {
            numPlayers = 2;
        }
        else {
            numPlayers = 1;
        }
        this.app.setState(ENGINE.Game)
    },
};

ENGINE.Pause = {
    smoothing: false,
    
    render: function(dt) {
        this.app.layer.clear("#60d")
                      .fillStyle("#fff")
                      .font("32px Verdana")
                      .fillText("Click to continue", 32, 128);
    },
    
    pointerdown: function(event) {
        this.app.setState(ENGINE.Game)
    },
};

////////////////////////////////////////////////////////////////////////////////

let gridHeightBits = 3;
let gridHeight = Math.pow(2, gridHeightBits) - 1;
let gridWidthBits = 5;
let gridWidth  = Math.pow(2, gridWidthBits);

let itoa = function(r) {
    var indices = [];
    for (var i = r.start; i < r.end; i++) {
        indices.push(i);
    }
    return indices;
}
// gameplay state
ENGINE.Game = {
    smoothing: false,
    
    // on first time state enter
    create: function() {
        // game output state
        this.recalculateGrid();
        // game input state
        this.leftDownPressed = false;
        this.leftUpPressed = false;
        this.rightDownPressed = false;
        this.rightUpPressed = false;
        this.mouseDown = false;
        // empty state
        this.state = qstate();
        // determine structure of the state
        this.paddleLength = 2;
        this.paddleGap = 1;
        
        this.structure = [
            { name: "x", len: gridWidthBits },
            { name: "y", len: gridHeightBits },
            { name: "paddleMinLeft", len: gridHeightBits },
            { name: "paddleMinRight", len: gridHeightBits },
            { name: "r", len: 1 },
            { name: "d", len: 1 },
            { name: "buttonBits", len: 4 },
            //{ name: "win", len: 1 },
        ];
        this.ranges = {};
        var pos = 0;
        for (var i = 0; i < this.structure.length; i++) {
            var r = { index: i, start: pos, end: pos + this.structure[i].len }
            pos = r.end;
            this.ranges[this.structure[i].name] = r;
        }
        // create a circuit that acts on this structure
        let _this = this;
        let indicesOf = function(name) {
            let r = _this.ranges[name];
            return itoa(r);
        };
        
        this.obsCoords = [
            { x: Math.floor(gridWidth*2/7),   y: Math.floor(gridHeight/2)+1, dir: "v", gate: "H" },
            { x: Math.floor(gridWidth*3/7),   y: Math.floor(gridHeight/2)-1, dir: "h", gate: "Y" },
            { x: Math.floor(gridWidth*4/7),   y: Math.floor(gridHeight/2)+1, dir: "h", gate: "H" },
            { x: Math.floor(gridWidth*5/7),   y: Math.floor(gridHeight/2)-1, dir: "v", gate: "Y" },
        ];
        
        if (numPlayers == 1) {
            this.leftAI = true;
            this.circuit = Qong.fullStepAI(
                gridHeight,
                gridWidth,
                this.paddleGap,
                this.paddleLength,
                indicesOf("paddleMinLeft"),
                indicesOf("paddleMinRight"),
                indicesOf("x"),
                indicesOf("y"),
                indicesOf("r"),
                indicesOf("d"),
                indicesOf("buttonBits"),
                this.obsCoords
            );
        }
        else if (numPlayers == 2) {
            this.leftAI = false;
            this.circuit = Qong.fullStep(
                gridHeight,
                gridWidth,
                this.paddleGap,
                this.paddleLength,
                indicesOf("paddleMinLeft"),
                indicesOf("paddleMinRight"),
                indicesOf("x"),
                indicesOf("y"),
                indicesOf("r"),
                indicesOf("d"),
                indicesOf("buttonBits"),
                this.obsCoords
            );
        }
        // initialize to new game state
        this.newGameState();
        
        if (false) {
            // print the circuit
            let tLog = new Array(this.state.n);
            for (var i = 0; i < this.state.n; i++) {
                tLog[i] = [];
            }
            this.circuit(circuitLogger(this.state.n, function(arr) {
                for (var i = 0; i < arr.length; i++) {
                    tLog[i].push(arr[i]);
                }
            }));
            
            console.log("---");
            for (var line of tLog) {
                console.log("[" + line.join("") + "]");
            }
            console.log("---");
        }
    },
    
    newGameState: function() {
        this.stop = false;
        
        this.timeSinceLastStep = 0;
        this.timePerStep = 1/8.0;
        this.lastStepDT = 0;
        
        this.probContinue = 1;
        this.leftProbLose = 0;
        this.rightProbLose = 0;
        
        this.leftTotalProbLost = 0;
        this.rightTotalProbLost = 0;
        
        let initState = "".concat(
            intBits(Math.floor(gridHeight/2), this.structure[0].len),
            intBits(Math.floor(gridHeight/2), this.structure[1].len),
            intBits(Math.floor(gridWidth/4)+1, this.structure[2].len),
            intBits(Math.floor(0), this.structure[3].len),
            "1",
            "0",
            "0000"
        );
        this.state.resetTo(initState, cOne());
    },
    
    // on entering this state
    enter: function() {
        
    },
    
    leave: function() { },

    // keyboard events
    keydown: function(event) {
        if (event.key == "s") {
            this.leftDownPressed = true;
        }
        if (event.key == "w") {
            this.leftUpPressed = true;
        }
        if (event.key == "down") {
            this.rightDownPressed = true;
        }
        if (event.key == "up") {
            this.rightUpPressed = true;
        }
        if (event.key == "i") {
            this.newGameState();
        }
    },
    keyup: function(event) {
        if (event.key == "s") {
            this.leftDownPressed = false;
        }
        if (event.key == "w") {
            this.leftUpPressed = false;
        }
        if (event.key == "down") {
            this.rightDownPressed = false;
        }
        if (event.key == "up") {
            this.rightUpPressed = false;
        }
    },
    
    pointerCheck: function() {
        let hspace = this.gstep * (gridWidth / 4);
        let vmid   = this.gridTop + this.gstep * (gridHeight / 2);
        
        this.leftUpPointed = false;
        this.leftDownPointed = false;
        this.rightUpPointed = false;
        this.rightDownPointed = false;
        
        for (var id in this.app.pointers) {
            let event = this.app.pointers[id];
            
            if (! event.pressed) continue;
            
            if (event.x - this.gridLeft <= hspace) {
                // left player
                if (event.y < vmid) {
                    // up
                    this.leftUpPointed = true;
                }
                else {
                    // down
                    this.leftDownPointed = true;
                }
            }
            else if (this.gridRight - event.x <= hspace) {
                // right player
                if (event.y < vmid) {
                    // up
                    this.rightUpPointed = true;
                }
                else {
                    // down
                    this.rightDownPointed = true;
                }
            }
        }
    },
    
    // pointers (mouse and touches)
    pointerdown: function(event) {
        //this.app.setState(ENGINE.Pause);
        if (this.stop) {
            if (   this.restartButton.x <= event.x
                && event.x <= this.restartButton.x + this.restartButton.width
                && this.restartButton.y <= event.y
                && event.y <= this.restartButton.y + this.restartButton.height) {
                // restart
                this.newGameState();
            }
        }
    },
    //pointerup: function(event) { },
    //pointermove: function(event) { },
    
    // mouse trap
    mousedown: function(event) { this.mouseDown = true; },
    mouseup: function(event) { this.mouseDown = false; },
    //mousemove: function(event) { },
    
    // finger trap - ouch
    //touchstart: function(event) { },
    //touchend: function(event) { },
    //touchmove: function(event) { },
    
    // gamepad events
    //gamepaddown: function(event) { },
    //gamepadup: function(event) { },
    //gamepadmove: function(event) { },
    getButtonBits: function() {
        // key events are recorded as they happen, but
        // must read the current state of any pointers
        this.pointerCheck();
        
        // apply inputs (simulating a joystick)
        let leftDown  =    ! this.leftAI
                        &&  (this.leftDownPressed  || this.leftDownPointed)
                        && !(this.leftUpPressed    || this.leftUpPointed);
        let leftUp    =    ! this.leftAI
                        &&  (this.leftUpPressed    || this.leftUpPointed)
                        && !(this.leftDownPressed  || this.leftDownPointed);
        let rightDown =     (this.rightDownPressed || this.rightDownPointed)
                        && !(this.rightUpPressed   || this.rightUpPointed);
        let rightUp   =     (this.rightUpPressed   || this.rightUpPointed)
                        && !(this.rightDownPressed || this.rightDownPointed);
        var newButtonBits = [
            leftDown  ? "1" : "0",
            leftUp    ? "1" : "0",
            rightDown ? "1" : "0",
            rightUp   ? "1" : "0"
        ];
        return newButtonBits;
    },
    
    step: function(dt) {
        if (this.stop) {
            return;
        }
        if (true) {
            this.timeSinceLastStep += dt;
            if (this.timeSinceLastStep < this.timePerStep) {
                return;
            }
            this.timeSinceLastStep = 0;
            this.lastStepDT = this.timePerStep;
        }
        else {
            this.lastStepDT = dt;
        }
        let newButtonBits = this.getButtonBits();
        for (var i = 0; i < newButtonBits.length; i++) {
            this.state.forceSet(this.ranges.buttonBits.start + i, newButtonBits[i]);
        }
        
        // compute circuit
        this.circuit(this.state);
        
        // determine if game is over
        let _this = this;
        let curProbLeftAlive = this.state.postSelect(function(tensor) {
            for (var i = _this.ranges.x.start; i < _this.ranges.x.end; i++) {
                if (tensor[i] == "1") return true;
            }
            return false;
        });
        this.leftTotalProbLost += (1 - curProbLeftAlive);
        let leftProbLoseDelta = (1 - curProbLeftAlive) * this.probContinue;
        this.leftProbLose += leftProbLoseDelta;
        this.probContinue -= leftProbLoseDelta;
        
        //if (feq(this.probContinue, 0)) {
        if (this.probContinue <= 0) {
            this.stop = true;
            this.probContinue = 0;
            return;
        }
        
        let curProbRightAlive = this.state.postSelect(function(tensor) {
            for (var i = _this.ranges.x.start; i < _this.ranges.x.end; i++) {
                if (tensor[i] == "0") return true;
            }
            return false;
        });
        this.rightTotalProbLost += (1 - curProbRightAlive);
        let rightProbLoseDelta = (1 - curProbRightAlive) * this.probContinue;
        this.rightProbLose += rightProbLoseDelta;
        this.probContinue -= rightProbLoseDelta;
        
        //if (feq(this.probContinue, 0)) {
        if (this.probContinue <= 0) {
            this.probContinue = 0;
            this.stop = true;
        }
    },
    
    recalculateGrid: function() {
        this.textHeight = 16 * 5;
        let appLeft = 0;
        let appTop = this.textHeight;
        let appRight  = this.app.width - appLeft;
        let appBottom = this.app.height - appTop;
        let siNaive = Math.floor((appRight - appLeft) / gridWidth);
        let sjNaive = Math.floor((appBottom - appTop) / gridHeight);
        this.gstep = Math.min(siNaive, sjNaive);
        
        // align to center
        let gridSizeX = this.gstep * gridWidth;
        let gridSizeY = this.gstep * gridHeight;
        
        let gridPadX = (appRight - gridSizeX) / 2;
        let gridPadY = (appBottom - gridSizeY) / 2;
        
        this.gridLeft = gridPadX;
        this.gridTop = gridPadY;
        this.gridRight  = gridPadX + gridSizeX;
        this.gridBottom = gridPadY + gridSizeY;
        
        this.restartButton = {
            x: this.app.width / 4 + 8,
            y: 8,
            width: this.app.width/2 - 16,
            height: this.gridTop - 16
        };
    },
    
    render: function(dt) {
        let layer = this.app.layer;
        let _this = this;
        
        // clear canvas
        layer.clear("#ffffff");
        
        // determine grid location
        this.recalculateGrid();
        let gridTop = this.gridTop;
        let gridLeft = this.gridLeft;
        let gridRight = this.gridRight;
        let gridBottom = this.gridBottom;
        let gstep = this.gstep;
        let textHeight = this.textHeight;
        
        // render a grid
        if (true) {
            layer.lineWidth(1).strokeStyle("#ff0000");
            for (var i = 0; i <= gridWidth; i++) {
                layer.strokeLine(gridLeft + i * gstep, gridTop, gridLeft + i * gstep, gridBottom);
            }
            for (var j = 0; j <= gridHeight; j++) {
                layer.strokeLine(gridLeft, gridTop + j * gstep, gridRight, gridTop + j * gstep);
            }
        }
        
        // render ball
        if (true) {
            let xRange = this.ranges.x;
            let yRange = this.ranges.y;
            let dBit = this.ranges.d.start;
            let rBit = this.ranges.r.start;
            
            layer.strokeStyle("#000")
                 .fillStyle("#000");
            
            // down-left arrow
            let arrow = [
                [gstep * -0.50, gstep * -0.50],
                [gstep *  0.25, gstep * -0.50],
                [gstep * -0.50, gstep *  0.25]
            ];
            
            this.state.forEach(function(tensor, amp) {
                let prob = amp.mag();
                let x = intFromBits(tensor.slice(xRange.start, xRange.end));
                let y = intFromBits(tensor.slice(yRange.start, yRange.end));
                let d = tensor.charAt(dBit);
                let r = tensor.charAt(rBit);
                
                layer.save()
                     .translate(gridLeft + gstep * (x + 0.5), gridTop + gstep * (y + 0.5));
                
                if (d == "1" && r == "0") {
                    // up-left => rotate 270
                    layer.rotate(-Math.PI/2);
                }
                else if (d == "1" && r == "1") {
                    // up-right, rotate 180
                    layer.rotate(Math.PI);
                }
                else if (d == "0" && r == "1") {
                    // down-right, rotate 90
                    layer.rotate(Math.PI/2);
                }
                //else if (d == "0" && r == "0") {
                    // down-left, don't rotate
                //}
                
                layer.a(prob)
                     .fillPolygon(arrow);
                
                layer.restore();
            });
            
            layer.a(1.0);
        }
        
        // render paddles
        if (true) {
            layer.fillStyle("#000");
            
            let lRange = this.ranges.paddleMinLeft;
            let rRange = this.ranges.paddleMinRight;
            
            let leftProbs = new Array(gridHeight);
            let rightProbs = new Array(gridHeight);
            for (var j = 0; j < gridHeight; j++) {
                leftProbs[j]  = 0.0;
                rightProbs[j] = 0.0;
            }
            
            this.state.forEach(function(tensor, amp) {
                let prob = amp.mag();
                let lPos = intFromBits(tensor.slice(lRange.start, lRange.end));
                let rPos = intFromBits(tensor.slice(rRange.start, rRange.end));
                
                for (var j = 0; j < _this.paddleLength; j++) {
                    leftProbs[lPos + j] += prob;
                    rightProbs[rPos + j] += prob;
                }
            });
            
            for (var j = 0; j < gridHeight; j++) {
                if (leftProbs[j] > 0) {
                    layer.a(Math.min(1.0, leftProbs[j]))
                     .fillRect(gridLeft + gstep * (_this.paddleGap + 0.5), gridTop + gstep * j,
                               gstep/2, gstep);
                }
                if (rightProbs[j] > 0) {
                    layer.a(Math.min(1.0, rightProbs[j]))
                     .fillRect(gridLeft + gstep * (gridWidth - _this.paddleGap - 1), gridTop + gstep * j,
                               gstep/2, gstep);
                }
            }
            
            layer.a(1.0);
        }
        
        // render obsCoords
        if (true) {
            let alpha = 0.5;
            
            let linesH = [
                [-0.3*gstep, -0.4*gstep,-0.3*gstep, 0.4*gstep],
                [-0.3*gstep,          0, 0.3*gstep,         0],
                [ 0.3*gstep, -0.4*gstep, 0.3*gstep, 0.4*gstep],
            ];
            let linesY = [
                [-0.3*gstep, -0.4*gstep,   0,          0],
                [ 0.3*gstep, -0.4*gstep,   0,          0],
                [         0,          0,   0,  0.4*gstep],
            ];
            
            layer.font("" + gstep + "px Verdana")
                 .textAlign("center")
                 .a(1.0);
            
            let fontMargin = Math.floor((gstep - layer.fontHeight()) / 2);
            
            for (var desc of this.obsCoords) {
                let color = (desc.dir == "v") ? "#00ff00" : "#0000ff";
                
                layer.strokeStyle("#000")
                     .strokeRect(gridLeft + gstep * desc.x, gridTop + gstep * desc.y,
                                 gstep, gstep);
                
                layer.strokeStyle(color);
                
                layer.save()
                     .translate(gridLeft + gstep * (desc.x + 0.5), gridTop + gstep * (desc.y + 0.5));
                
                if (desc.dir == "h") {
                    layer.rotate(Math.PI/2);
                }
                let lines = (desc.gate == "H") ? linesH : linesY;
                
                for (var line of lines) {
                    layer.strokeLine(line[0], line[1], line[2], line[3]);
                }
                
                layer.restore();
            }
        }
        
        // render prob lose text
        if (true) {
            layer.strokeStyle("#000")
                 .strokeLine(gridLeft, gridBottom, gridRight, gridBottom);
            
            layer.fillStyle("#000")
                 .font("1em Verdana")
                 .textAlign("left")
                 .fillText("Rate: " + Math.floor(1 / this.lastStepDT),
                           this.app.width / 4, gridBottom + 16 * 2)
                 .fillText("Left lost " + percent(this.leftTotalProbLost) +
                           "% in " + percent(this.leftProbLose) + "% of worlds",
                           this.app.width / 4, gridBottom + 16 * 3)
                 .fillText("Right lost " + percent(this.rightTotalProbLost) +
                           "% in " + percent(this.rightProbLose) + "% of worlds",
                           this.app.width / 4, gridBottom + 16 * 4)
                .fillText("Playing in " + percent(this.probContinue) + "% of worlds",
                           this.app.width / 4, gridBottom + 16 * 5)
                ;
        }
        
        // render try again button
        if (this.stop) {
            layer.save()
                 .strokeStyle("#000")
                 .lineWidth(3)
                 .strokeRect(this.restartButton.x, this.restartButton.y, this.restartButton.width, this.restartButton.height)
                 .fillStyle("#ddd")
                 .fillRect(this.restartButton.x, this.restartButton.y, this.restartButton.width, this.restartButton.height)
                 .restore();
            
            layer.fillStyle("#000")
                 .font("1em Verdana")
                 .textAlign("left")
                 .fillText("Touch here or press I to reinitialize!",
                           this.restartButton.x + 16, this.restartButton.y + this.restartButton.height/2)
        }
    }
};

////////////////////////////////////////////////////////////////////////////////

// disables transition animation
//PLAYGROUND.Transitions.plugin = false;

// entrance point
var app = playground({
    smoothing: false, // no interpolation
    scale: 1,         // force scaling
    
    // states related events (called only for application)
    //createstate: function() { },
    //enterstate: function() { },
    //leavestate: function() { },
    
    // before create
    //preload: function() { },
    
    // main loader
    create: function() {
        
    },
    
    // after create
    ready: function() {
        this.setState(ENGINE.Menu)
    }
});
