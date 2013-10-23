/*
 * Command line interface to plast
 *
 * Ben Postlethwaite
 */


var plast = require("../")()
  , fs = require('fs')
  , argv = require('optimist')
           .demand(["f"])
           .argv;

run(argv)

function run(argv) {

    getStdin( function (stdindata) {

        var depth = argv.d ? argv.d : 0
          , modname = argv.w
          , src = fs.readFileSync(argv.f, {encoding:'utf-8'})
          , funclist = [];


        /*
         * g option, hunt global functions.
         * Depth default is 0, set it on the command line
         */

        if (argv.g) {


            funclist = getfuncs(src, depth)

            /*
             * If user doesn't want to wrap, print and call it a day
             */
            if (!argv.w) {
                funclist.forEach( function (f) {console.log(f)} );
                process.exit()
            }
        }

        /*
         * If we got no data from global hunting lets see if the user
         * is piping in function names to wrap.
         */
        if (argv.w) {
            funclist = argv.g ? funclist : stdindata;
            var newsrc = wrapfuncs(src, funclist, modname);
            console.log(newsrc);
        }

    })
}

function getStdin ( cb ) {
    /*
     * First lets check for stdin, and grab it all
     */
    var chunks = '';
    if (!process.stdin.isTTY) {

        process.stdin.resume();
        process.stdin.setEncoding('utf8');

        process.stdin.on('data', function(chunk) {
            chunks += chunk;
        });

        process.stdin.on('end', function () {
            chunks.replace( /\n/g, " " ).split( " " )
            cb(chunks)
        })
    }
    else cb(chunks);

}
