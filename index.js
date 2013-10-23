/*
 * plast - the Plotly AST thing
 *
 * Daisy chain plast with pipes to get it done!
 *
 *
 * >> plast -g -d 0 -f index.js | plast -w "wrapper" -f cmd.js > outfile.js
 *
 * Ben Postlethwaite
 * Oct 22nd
 */

var falafel = require('falafel');

module.exports = plast;

function plast() {

    var self = {}


    function getfuncs(src, depth) {
        var fl = [];
        falafel(src, function (node) {
            if (node.type === "FunctionDeclaration")
                if (getDepth(node) === depth)
                    fl.push(node.id.name)
        });
        return fl;
    }

    function wrapfuncs(src, funclist, modname) {

        return falafel(src, function (node) {
                   if (node.type === "CallExpression")
                       if (funclist.indexOf(node.callee.name) > -1)
                           node.update( modwrap(node, modname) )
               });
    }

    function modwrap(node, modname) {
        var src = node.source();
        var reg = new RegExp(node.callee.name);
        src = src.replace(reg, modname + "." + node.callee.name);
        return src;
    }


    function getDepth(node, depth) {

        if (node.type === "Program") return -1;

        if (!depth) depth = 0;
        var np = node.parent;

        if (np.type !== "Program")
            if (np.type ==="FunctionDeclaration") {
                depth += 1
                return getDepth(np, depth);
            } else
                return getDepth(np, depth);
        else return depth;
    }

    self.getfuncs = getfuncs;
    self.wrapfuncs = wrapfuncs;

    return self;
}