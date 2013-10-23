function fncfitter() {
function fitplotdata(){
        var td=Tabs.get();
        if( td.data === undefined ){
            Plotly.Lib.notifier('Upload some data to make a graph!');
            return;
        }

        // make popover with trend line options
        var fitlink=$(td).find('.fitlink')[0];
        newPopover(td,fitlink,'fitlink',fitBoxContent);

        var id = ".context-id-" + td.id; // get popup specific to this tab

        InitFitTabs(td, $(id).find(".fitopts"), 'rerunfit(this);');

        // special bootstrap styled <select> dropdown
        var sp = $('.selectpicker');
        if( sp.length !== 0 ) sp.selectpicker();

        // create fit, if fit doesn't exist yet
        if(!plotfitexists(td)) { rerunfit(); }
    }

    // checks whether a fit already exists in the plot by looking for trace names
    // preceded by "fit_"
    function plotfitexists(td){
        var fitexists = false;
        $(td.data).each(function(i,v){ if(v.name.substr(0,4) == 'fit_') fitexists = true; });
        return fitexists;
    }

    // creates a trend line and adds an annotation
    // xold - original x array (mysql dates instead of unix timestamps, used for plotting)
    // xnew - cleaned x array (used for fit calculation)
    // y - y data, used to calculate fit
    // order - polynomial order, used to calculate fit
    function createplotfit(td, xnew, fit, params, rsquared){
        // Ben note: I am putting in xnew into the plot now... Is there a reason not too?
        // make a linear fit the data
        /*console.log('FIT DATA');
          console.log(xnew,y,order);*/
        // colorbrewer Set1 set
        var anno = {}
          , fitdata = {}
          , fit_colors = [ 'rgb(228, 26, 28)', 'rgb(55, 126, 184)', 'rgb(77, 175, 74)',
                           'rgb(52, 78, 163)', 'rgb(255, 127, 0)', 'rgb(66, 86, 40)',
                           'rgb(166, 86, 40)', 'rgb(247, 129, 191)', 'rgb(53, 153, 153)' ]
          , data_index = Number($('select.selectfitdata').val())
          , tracename = td.data[data_index].name
          , numfits = 0

        // how many fits are thre?

        $(td.data).each(function(i,v){ if(v.name.substr(0,4)=='fit_') numfits++; });
        if( numfits > 8 ) numfits = Math.floor(Math.random()*8);

        anno['text'] = write_fit_eqn(params) +
            '<br><br>R<sup>2</sup> = ' + rsquared.toString().substr(0,6);

        var xfloor = Math.min.apply(null, xnew)
          , yfloor = Math.min.apply(null, fit);
        var xmiddle = ((Math.max.apply(null,xnew) - xfloor) / 2) + xfloor;
        var direction = 1;

        if (numfits % 2 === 0) direction = -1;
        anno['x'] = xmiddle + ( xmiddle * 0.3 * numfits * direction );
        // offset subsequent annotations so not
        // right on top of each other if multiple

        if( td.layout.xaxis.type == 'date' ) anno['x'] = unix_to_mysql(anno['x']);

        anno['y'] = ((Math.max.apply(null, fit) - yfloor) / 2) + yfloor;
        anno['arrowwidth'] = 3;
        anno['font'] = { size:16 };
        anno['arrowcolor'] = fit_colors[numfits+1];
        anno['tracename'] = 'fit_' + tracename;

        // make the annotation
        Plotly.Annotations.draw(td,undefined,anno,'add');

        fitdata['line'] = {color: fit_colors[numfits+1], width: 6};
        fitdata['opacity'] = 0.5;
        fitdata['name'] = 'fit_' + td.data[data_index].name;
        fitdata['x'] = xnew;
        fitdata['y'] = fit;

        // make the fit line!
        Plotly.plot(td,[fitdata]);
    }

    // remove fit for selected trace
    // identify fits by names starting with "fit_"
    // remove only one at a time
    function removeplotfits(td){

        td.changed = true;

        var data_index = Number($('select.selectfitdata').val());
        var tracename = td.data[data_index].name;

        $(td.data).each(function(i,v){
            // ~ TRACE ~
            if(v.name.substr(0,4) == 'fit_'  && v.name.substr(4) == tracename ){
                td.data.splice(i,1); // REMOVE FROM gd.data
                td.calcdata.splice(i,1); // REMOVE FROM gd.calcdata
                // ~ ANNOTATION ~
                $(td.layout.annotations).each(function(j,u){
                    if( 'tracename' in u ){
                        if( u.tracename == 'fit_' + tracename ){
                            td.layout.annotations.splice(j,1); // REMOVE FROM gd.layout.annotations
                        }
                    }
                });
            }
        });
        Plotly.plot(td);
    }

    // This function generates the html the fills the trend link popover
    function fitBoxContent(popover){

        var gd = popover[0].gd
          , pc = popover.find('.popover-content>div')
          , tracedd = '<select class="selectpicker selectfitdata">' // html for trace dropdown
          , pt = popover.find('.trace-dropdown').html('Make Linear or Polynomial Fits');


        $(gd.data).each(function(i,v){
            if( v.name.substr(0,4) != 'fit_' ) tracedd+='<option value="'+i+'">'+v.name+'</option>';
        });
        tracedd += '</select><div class="sublabel" style="margin-top:-10px;">Select Data to Fit</div><br>';

        var html =
            '<div style="border-right:1px solid #444;padding-right:10px;">'+
            // Left-Hand Side <div>
            '<div style="margin-top:12px;">'+tracedd+'</div>'+
            '<div style="display:inline-block;margin-right:10px;">'+
            '<a class="btn btn-success" '+
            'onclick="rerunfit();">'+
            '<i class="icon-plus"></i> Fit</a>'+
            '</div>'+
            '<div style="display:inline-block;">'+
            '<a class="btn btn-danger" '+
            'onclick="removeplotfits(Tabs.get());">'+
            '<i class="icon-remove"></i> Remove Fit</a>'+
            '</div>'+
            '<div class="popover-footer">'+
            '<b style="color:#3a87ad;">Pro Tip:</b>'+
            ' You can also make fits with the '+
            '<a onclick="graphToGrid(\'Fit\');$(\'div.popover\').remove();" '+
            'style="color:#3a87ad;cursor:pointer!important;text-decoration:none;">'+
            '<b>Plotly grid</b></a>'+
            '</div>'+
            '</div>'+
            // Right-Hand Side <span>
            '<span class="fitopts" style="float:right;">'+
            '</span>'+
            '</div>';

        pc.html( html );

        // special bootstrap styled <select> dropdown
        var sp = $('.selectpicker');
        if( sp.length !== 0 ) sp.selectpicker();
    }

    // This function looks for the selected trace in gd.data with a name starting with fit_
    // It then removes it from gd.data, and makes a new fit from the params in the fit popover
    function rerunfit(el){

        addMessage('Generating Fit...');
        var td = Tabs.get();
        removeplotfits(td);
        var data_index = Number($('select.selectfitdata').val());
        if( isNaN( data_index ) ) data_index = 0;
        var order=Number($('body>.popover').find('.porder').val());
        if( isNaN( order ) ) order = 1;

        var x = td.data[data_index].x
          , y = td.data[data_index].y
          , xnew = mysql_to_unix(x);

        var data = applyFitParameters(td, xnew, y, x);

        var xfit = data.xfit
          , fit = data.fit
          , p = data.fitdata

        createplotfit(td, xfit, fit, p.params, p.corr)

        delMessage('Generating Fit...');
    }


    function applyFitParameters(td, x, y, x_old){


        // Make copies
        var x_cp = x.slice(0)
          , y_cp = y.slice(0)
          , xchanged = false    // hack for signalling to grid if we need to show a new x column

        // Strip out leading and trailing headers and garbage: Mutates passed in arrays (pass by ref)
        stripNaNs(x_cp, y_cp)

        // Get the fit parameters from the dom elements
        var fitp = getFitParameters(td, x_cp, x_old);

        var mn = fitp.mn
          , mx = fitp.mx
          , polyorder = fitp.polyorder
          , nintrp = fitp.nintrp;

        // Sort arrays - makes subRanges, interpolation and extrapolation much easier.
        // We need to let the Grid know if this thing was sorted so it can show a new x column
        // I would like to have another work around, but this will do for now... so compare against
        // a tmp array
        var z = d3.zip(x_cp, y_cp).sort( function (a, b) {return a[0] - b[0]} );
        var x_tmp = [];
        z.forEach( function (zip, idx) {
            x_tmp[idx] = zip[0];
            y_cp[idx] = zip[1];
        })


        if (!arrays_equal(x_tmp, x_cp))
            xchanged = true

        x_cp = x_tmp;

        var xo = isNaN(mn) ? 0 : findInsertionIndex(x_cp, mn);
        var xf = isNaN(mx) ? (x_cp.length - 1) : findInsertionIndex(x_cp, mx);

        // take subarrays
        x_cp = x_cp.splice(xo, xf - xo + 1);
        y_cp = y_cp.splice(xo, xf - xo + 1);

        // fit the data
        var p = polyfit(x_cp, y_cp, polyorder);

        var fit = p.fit,
            xfit = x_cp.slice(0);

        /*
         * If mn and mx from subarray are outside domain (this function checks this internally)
         * then assign copied arrays to variables. Pass in fit function, since this is what we
         * are appending to, not the actual function y.  returns {'xnew': xnew, 'ynew': ynew}
         */
        var xtrp = extrapolate(xfit, fit, [mn, mx], p.params);
        xfit = xtrp.xnew;
        fit = xtrp.ynew;
        if (xtrp.extrapd) {
            console.log("performing extrapolation")
            xchanged = true
        }

        var itrp = interpolate(xfit, fit, p.params, nintrp);
        xfit = itrp.xnew;
        fit = itrp.ynew;
        if (itrp.interpd) {
            console.log("performing interpolation")
            xchanged = true
        }

        return {
            x: x_cp
          , y: y_cp
          , xfit: xfit
          , fit: fit
          , fitdata: p
          , xo: xo
          , xf: xf
          , xchanged: xchanged
        };
    }

    function getFitParameters(td, x, xold) {
        // Parse min and max if they are there: Note this grabs first min and max
        // IF user has other plot or grid tabs open with min and max, how does it know
        // which min and max to grab?

        var id = ".context-id-" + td.id;

        var mn=$(id).find('.xmin').val()
          , mx=$(id).find('.xmax').val()
          , polyorder = clean_num($(id).find('.porder').val())
          , nintrp;

        if (td.tabtype === "grid")
            nintrp = clean_num($(id).find('.interp').val())

        else if (td.tabtype === "plot")
            nintrp = clean_num($(id).find('#interp-slider').slider("option", "value"))

        else
            console.log("could not determine tab type")

        // if xmin and xmax are datetime strings, convert them to a unix timestamp
        // TODO - convert this to
        if( mn.match(/^\d\d\d\d-(\d)?\d-(\d)?\d \d\d:\d\d:\d\d$/g) !== null ){
            mn=Date.parse(mn)/1000 - Date.parse(xold[0])/1000;
        }
        if( mx.match(/^\d\d\d\d-(\d)?\d-(\d)?\d \d\d:\d\d:\d\d$/g) !== null ){
            mx=Date.parse(mx)/1000 - Date.parse(xold[0])/1000;
        }

        // check for exponential notation
        if(typeof(mn)==='string'){ if(mn.indexOf('e')>0 || mn.indexOf('E')>0) mn=Number(mn); }
        if(typeof(mx)==='string'){ if(mx.indexOf('e')>0 || mx.indexOf('E')>0) mx=Number(mx); }

        // take out any accidental non-numeric characters in min and max
        mn=clean_num(mn);
        mx=clean_num(mx);

        // if xmin or xmax are empty strings, set them to 1st number of array ends
        if( mn==='' || mn===undefined ) mn = x[0];
        if( mx==='' || mx===undefined ) mx = x[x.length - 1];


        // switch min and max if min is greater than max
        if(!isNaN(mn) && !isNaN(mx)){
            if(mn>mx){
                var mn_cp=mn;
                mn=mx;
                mx=mn_cp;
            }
        }


        return {
            "mn" : mn
          , "mx" : mx
          , "polyorder": polyorder
          , "nintrp" : nintrp
        };

    }



    function interpolateXRange(x, Nps) {
        /*
         * Interpolate Nps values between the elements of x.
         * Returns an array of array. Each array are the interpolated
         * values between each x elements. So that there are x.length - 1 arrays
         * returned.
         */

        var i, new_x = []
          , nx = x.length
          , xrange = x[nx - 1] - x[0]
          , interps = [];

        /*
         * Relatively space interpolation points between x's
         * for case where x is not uniform
     */
    var npts = []
    for (i = 0; i < nx - 1; i++)
        npts[i] = Math.floor( ((x[i+1] - x[i]) / xrange) * Nps);

    for (i = 0; i < nx - 1; i++) {
        interps = interps.concat(numeric.linspace(x[i], x[i+1], npts[i] + 2))
        // Pop off the ends of inside x points to avoid doubles
        if (i < nx - 2)
            interps.pop();
    }

        return interps;
    }

    function PolySolver(params) {
        return function (x) {
            return params.reduce( function(prev, curr, j) {
                       return prev + curr * Math.pow(x, j);
                   });
        };
    }


    function interpolate(x, fit, params, Nps) {
        /*
         * Wrapper around xpolate that interpolates values
         */

        if (isNaN(Nps) || Nps === '' || Nps < 1)
            return {'xnew': x, 'ynew': fit, 'interpd': false};

        //var xnew = x.slice(0);

        /*
         * Sort xnew and ynew
         */
        //xnew.sort( function (a, b) {return a - b} )
        var i
          , xinterps = interpolateXRange(x, Nps)

        if (xinterps.length === x.length)
            return {'xnew': x, 'ynew': fit, 'interpd': false};

        var solver = PolySolver(params);
        var yinterps = xinterps.map(solver);

        return {'xnew': xinterps, 'ynew': yinterps, 'interpd': true};
    }


    function extrapolate(x, y, xtraps, params) {
        /*
         * Performs extrapolation for a given array of new points
         *
         * y = c0 + c1*x + c2*x^2 + ... cN*x^N
         *
         * Using coefficients from polynomial fit.
         * It returns new x array and y array containing xpolated points and
         * old points. If you pass in xtrap points inside current x domain it
         * returns the interpolated points. It inserts points in correct positions
         * in copied array.
         * It has an optional yold argument, that if supplied, is extended and padded
         * along with x and yfit and returned as ypad.
         */
        if (!Array.isArray(xtraps))
            xtraps = [xtraps]
        if (!Array.isArray(params))
            xtraps = [params]

        xtraps = xtraps.filter( function (xt) {
                     if (!isNaN(xt) && (xt < x[0] || xt > x[x.length - 1]))
                         return true
                     else return false
                 })

        if (xtraps.length === 0)
            return {'xnew': x, 'ynew': y, 'extrapd': false};


        var i;
        var solver = PolySolver(params);
        var ytraps = xtraps.map(solver);

        /*
         * Insert extrap/interp points in the right order into new x and y arrays
         */
        var xnew = x.slice(0);
        var ynew = y.slice(0);

        for (i = 0; i < xtraps.length; i++)  {
            var idx = findInsertionIndex(xnew, xtraps[i]);
            xnew.splice(idx, 0, xtraps[i]);
            ynew.splice(idx, 0, ytraps[i]);
        }

        // console.log("xnew in xpolate", xnew)
        // console.log("ynew in xpolate", ynew)

        return {'xnew': xnew, 'ynew': ynew, 'extrapd': true};

    }

    function findInsertionIndex(sortedArr, val) {
        /*
         * Find insertion point for a numeric val
         * Can easily generalize by making a compare function for arbitrary data types, but
         * using this implementation for efficiency
         * @param sortedArr The sorted array
         * @param val The value for which to find an insertion point (index) in the array
         */

        var low = 0, high = sortedArr.length;
        /*
         * quick precheck for extrapolated points
         */
        if (val < sortedArr[0])
            return low;
        if (val > sortedArr[high])
            return high;
        /*
         * Perform binary search for the index to insert at
         */
        var mid = -1, c = 0;
        while(low < high)   {
            mid = (low + high) >> 1; // fast Math.floor(low+high)/2
            if(sortedArr[mid] < val)   {
                low = mid + 1;
            }else if(sortedArr[mid] > val) {
                high = mid;
            }else {
                return mid;
            }
        }
        return low;
    }



    function InitInterpSlider (td) {

        var id = ".context-id-" + td.id;

        var nmin = 0//fitdata.x.length
          , nmax = 100
          , npoints = nmin
          , valstr = "Smooth with "
          , valstr2 = " extra points (optional)"
          , options = {
              value: npoints
            , min: nmin
            , max: nmax
            , step: 10
            , slide: function(event, ui) {
                  npoints = ui.value;
                  $(slider).slider( "option", "value", npoints);
                  $(display).text(valstr + npoints + valstr2);
                  rerunfit();
              }
          }
          , slider  = $("<div>")
                      .attr("id", "interp-slider")
                      .css("margin-left", "10px")
                      .css("margin-right", "10px")
                      .css("font-size", "small")
                      .slider(options)

          , display = $('<div>').attr('id',"interp-slider-value")
                      .addClass("sublabel")
                      .css("margin-left", "10px")
                      .text(valstr + npoints + valstr2)


        $(slider).find(".ui-slider-handle")
        .css("width","0.8em");

        $(id).find(".fit-tab-content").append(slider, display);


    }

    function InitFitTabs (td, elem, onchange, fitTabName) {


        var id = ".context-id-" + td.id;
        var tabname = ''

        $(elem).html(fitopts(td, onchange, fitTabName))

        $(id).find(".fit-tab")
        .click( function (e) {

            tabname = e.target.innerHTML
            InitFitTabs(td, elem, onchange, tabname)
        })


        if (td.tabtype === "plot")
            InitInterpSlider(td);
    }

    function functionfitter (options) {

        if (!options) options = {};
        var log = options.log || false

        var self = {};

        /*
         * Definitions for automating help
         */
        var funcdefs = {
            "max":"","min":"","abs":"","sqrt":"","exp":"",
            "ln":"","log10":"","log2":"","power":"","pow":"",
            "sin":"","cos":"","tan":"","cot":"","sec":"",
            "csc":"","atan":"","asin":"","acos":"","acot":"",
            "asec":"","acsc":"","sinh":"","cosh":"","tanh":"",
            "coth":"","sech":"","csch":"","asinh":"","acosh":"",
            "atanh":"","acoth":"","asech":"","acsch":"","sum":"",
            "fact":"","gamma":"","chisq":"","norm":"","gauss":"",
            "studt":"","pi":"","statcom":"","inverse":"","agauss":"",
            "anorm":"","aerf":"","achisq":"","astudt":"","afishf":""
        };

        var constdefs = {"pi":""};



        /*
         * Function and variable definitions
         */
        var PI=Math.PI, PID2=PI/2, Deg=180/PI, SQRT2=Math.sqrt(2);

        function MAX(x1,x2) { return Math.max(x1,x2); }
        function MIN(x1,x2) { return Math.min(x1,x2); }
        function ABS(x) { return Math.abs(x); }
        function SQRT(x) { return Math.sqrt(x); }
        function EXP(x) { return Math.exp(x); }
        function LN(x) { return Math.log(x); }
        function LOG10(x) { return LN(x)/Math.LN10; }
        function LOG2(x) { return LN(x)/Math.LN2; }
        function POWER(x,y) { return Math.pow(x,y); }
        function POW(x,y) { return Math.pow(x,y); }
        function SIN(x) { return Math.sin(x); }
        function COS(x) { return Math.cos(x); }
        function TAN(x) { return SIN(x)/COS(x); }
        function COT(x) { return COS(x)/SIN(x); }
        function SEC(x) { return 1/COS(x); }
        function CSC(x) { return 1/SIN(x); }
        function ATAN(x) { return Math.atan(x); }
        function ASIN(x) { return Math.asin(x); }
        function ACOS(x) { return Math.acos(x); }
        function ACOT(x) { return ATAN(1/x); }
        function ASEC(x) { return ACOS(1/x); }
        function ACSC(x) { return ASIN(1/x); }
        function SINH(x) { return (EXP(x)-EXP(-x))/2; }
        function COSH(x) { return (EXP(x)+EXP(-x))/2; }
        function TANH(x) { return SINH(x)/COSH(x); }
        function COTH(x) { return 1/TANH(x); }
        function SECH(x) { return 1/COSH(x); }
        function CSCH(x) { return 1/SINH(x); }
        function ASINH(x) { return LN(x+SQRT(x*x+1)); }
        function ACOSH(x) { return LN(x+SQRT(x*x-1)); }
        function ATANH(x) { return 0.5*LN((1+x)/(1-x)); }
        function ACOTH(x) { return 0.5*LN((x+1)/(x-1)); }
        function ASECH(x) { return LN(1/x+SQRT(1/(x*x)+1)); }
        function ACSCH(x) { return LN(1/x+SQRT(1/(x*x)-1)); }
        function SUM(arr) {
	    return arr.reduce(function(prev, curr) { return prev + curr; });
        }

        function FACT(n) { // extention of factorial to all real numbers
            if(n===0 | n==1) { return 1; }
            if(n<0) { return FACT(n+1)/(n+1); }
            if(n>1) { return n*FACT(n-1); }
            var r= (n<0.5) ? n : 1-n;
            r = 1 / (1 +
                     r*( 0.577215664819072 +
                         r*(-0.655878067489187 +
                            r*(-0.042002698827786 +
                               r*(0.166538990722800 +
                                  r*(-0.042197630554869 +
                                     r*(-0.009634403818022 +
                                        r*(0.007285315490429 +
                                           r*(-0.001331461501875)))))))));
            if( n > 0.5 ) { r = n*(1-n)*PI / (r*SIN(PI*n)); }
            return r;
        }
        function GAMMA(n) { return FACT(n-1); }

        function lpadzero(value, padding) { // courtesy http://samuelmullen.com/2012/03/left-pad-zeroes-in-javascript/
            var zeroes = "0";
            for (var i = 0; i < padding; i++) { zeroes += "0"; }
            return (zeroes + value).slice(padding * -1);
        }

        function CHISQ(x,n) {
            if(x>1000 | n>1000) {
                var q=NORM((POWER(x/n,1/3)+2/(9*n)-1)/SQRT(2/(9*n)))/2;
                return (x>n) ? q : 1-q;
            }
            var p=EXP(-0.5*x);
            if((n%2)==1) { p=p*SQRT(2*x/PI); }
            var k=n;
            while(k>=2){
                p=p*x/k;
                k=k-2;
            }
            var t=p, a=n;
            while(t>1e-15*p){
                a+=2;
                t*=x/a;
                p+=t;
            }
            return 1-p;
        }

        function NORM(z) {
            var q=z*z;
            if(ABS(z)>7) { return (1-1/q+3/(q*q))*EXP(-q/2)/(ABS(z)*SQRT(PID2)); }
            else { return CHISQ(q,1); }
        }

        function GAUSS(z) { return ((z<0) ? ((z<-10) ? 0 : CHISQ(z*z,1)/2 ) : ((z>10) ? 1 : 1-CHISQ(z*z,1)/2)); }
        function ERF(z) { return ((z<0) ? (2*GAUSS(SQRT2*z)-1) : (1-2*GAUSS(-SQRT2*z))); }

        function STUDT(t,n) {
            t=ABS(t);
            var w=t/SQRT(n), th=ATAN(w);
            if(n==1) { return 1-th/PID2; }
            var sth=SIN(th), cth=COS(th);
            if((n%2)==1) { return 1-(th+sth*cth*STATCOM(cth*cth,2,n-3,-1))/PID2; }
            else { return 1-sth*STATCOM(cth*cth,1,n-3,-1); }
        }

        function FISHF(f,n1,n2) {
            var x=n2/(n1*f+n2);
            if((n1%2)===0) { return STATCOM(1-x,n2,n1+n2-4,n2-2)*POWER(x,n2/2); }
            if((n2%2)===0) { return 1-STATCOM(x,n1,n1+n2-4,n1-2)*POWER(1-x,n1/2); }
            var th=ATAN(SQRT(n1*f/n2)), a=th/PID2, sth=SIN(th), cth=COS(th);
            if(n2>1) a=a+sth*cth*STATCOM(cth*cth,2,n2-3,-1)/PID2;
            if(n1==1) { return 1-a; }
            var c=4*STATCOM(sth*sth,n2+1,n1+n2-4,n2-2)*sth*POWER(cth,n2)/PI;
            if(n2==1) { return 1-a+c/2; }
            var k=2;
            while(k<=(n2-1)/2){
                c*=k/(k-0.5);
                k++;
            }
            return 1-a+c;
        }
        function STATCOM(q,i,j,b) {
            var zz=1, z=zz, k=i;
            while(k<=j){
                zz*=q*k/(k-b);
                z+=zz;
                k+=2;
            }
            return z;
        }

        function INVERSE(f,p){ // general inversion of f(x) for x>=0
            var v=0.5, dv=0.25, x;
            while(dv>1e-15){
                x=1/v-1;
                v+= (f(x)>p) ? -dv : dv;
                dv/=2;
            }
            return x;
        }


        function AGAUSS(p) {
            if(p>0.5) { return SQRT(ACHISQ(2*(1-p),1)); }
            else { return -SQRT(ACHISQ((2*p,1))); }
        }

        function ANORM(p) { return INVERSE(NORM,p); }
        function AERF(p) { return AGAUSS(p/2+0.5)/SQRT2; }
        function ACHISQ(p,n) { return INVERSE(function(a){return CHISQ(a,n);},p);}
        function ASTUDT(p,n) { return INVERSE(function(a){return STUDT(a,n);},p);}
        function AFISHF(p,n1,n2) { return INVERSE(function(a){return FISHF(a,n1,n2);},p);}

        // helper functions
        function Fmt(x) {
            var v;
            if(ABS(x)<0.00005) { x=0; }
            v = ' '+(x+((x>0) ? 0.00005 : -0.00005));
            v = v.substring(0,v.indexOf('.')+5);
            return v.substring(v.length-10,v.length);
        }

        function vFmt(x) {
            var v;
            if(ABS(x)<0.0000005) { x=0; }
            v = ' '+(x+((x>0) ? 0.0000005 : -0.0000005));
            v = v.substring(0,v.indexOf('.')+7);
            return v.substring(v.length-18,v.length);
        }

        function createArray() {
            this.length = createArray.arguments.length;
            for (var i = 0; i < this.length; i++) { this[i] = createArray.arguments[i]; }
        }

        function ix(j,k,p) { return j*(p+1)+k; }

        function cropNum(v,vmin,vmax,vdefault){
            if(typeof v != 'number') { return vdefault; }
            else { return MAX(vmin,MIN(vmax,v)); }
        }

        // --------------------------------------- Iterate -------------------------------------
        //
        // ***** input ******************
        // data: array of points
        //   each point is array of x1,..,xn, then y; or y, weight; or y1,..,yn
        // params: initial parameter values, to be referred to in the eqn as A-H, or P1-P8
        // npar: how many of these to vary (requires held params at the end)
        // nvar: number of independent variables (default 1)
        // pctile: percent of points to put below the curve, if using least abs scaling (default 50)
        // centered: use centered derivative method (slower but more accurate) instead of one-sided (default false)
        // transform: array of strings containing code to transform each variable (Y is first, then X1-X8)
        //   all variables are 'V', and missing vals go to just 'V'
        // funct: string of the function to fit, using the form of the reference syntax, described at the top of this file
        // errorform: what kind of y error to use. default "1"
        //   expected values are:
        //      "1": all errors=1
        //      "Y": errors=Y values, good for exponential-type data
        //      "Sqrt(Y)": poisson noise, good for count data
        //      "w": looks for an extra column in data with weight values
        //      "Rep": a bunch of Y values are on the line, get spread of these values
        // leastabs: fit to least absolute error, rather than least squares
        // frelax: fractional relaxation: default 1, use a smaller number to converge more slowly but robustly
        // cov: covariance matrix from a previous iteration, if there is one.
        //
        // ***** output ******************
        // params: [{val, err, pval},...]
        // fit: array of Y_calc values
        // res: array of Y_in-Y_calc
        // conf_lo: array of 95% lower confidence vals
        // conf_hi: array of 95% upper confidence vals
        // rms: average rms deviation of data (vs. error per errorform)
        // cov: covariance matrix (make sure to save this for future iterations)
        function Iterate(input) {
                         var xArr = new createArray(0,0,0,0,0,0,0,0,0);
                         var Par = input.params.slice(0);
                         var SEP = new createArray(1,1,1,1,1,1,1,1);
                         var Der = new createArray(0,0,0,0,0,0,0,0,0);
                         var Arr = new createArray(0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0);
                         var Cov;
                         if('cov' in input && input.cov) { Cov=input.cov; }
                         else { Cov = new createArray(0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0 ); }

                         var cccSW = 0, cccAv = 0, cccSD = 0;

                         var NL = unescape("%0D") + unescape("%0A"); // CR+LF
                         var i=0, j=0, k=0, l=0, m=0;

                         var nPts = input.data.length; // IN: # data points
                         var nPar = ('npar' in input) ? input.npar : input.params.length; // IN: # fit parameters
                         var nVar = cropNum(input.nvar,1,8,1); // IN: # independent variables
                         var dgfr = nPts - nPar;
                         var St95 = ASTUDT( 0.05 , dgfr );
                         var Pctile = cropNum(input.pctile/100,0,1,0.5); // IN: percent of data to try to put below the curve (only if using least abs)

                         // IN: Parameter initial values in array Par, aliased to A-H, P1-P8 for the sake of eval(Y)

                         // IN: evaluate the function
                         var func = makefunc(input.funct);

                         var SSq = 0;

                         // make the output object. include the input data for convenience
                         var out = {'params':[], 'fit':[], 'res':[], 'conf_lo':[], 'conf_hi':[], 'rms':0, 'cov':[]};

                         var X,yc,vSEy,w,nY,sY,sYY,Rep,yT,
                         Save,Del,ccSW,ccAv,ccSD,GenR2,GenR,FRelax;

                         for (j = 0; j<nPar*(nPar+1); j++) { Arr[j] = 0; }
                         // big loop over data points
                         for (i = 0; i<nPts; i++) {

                             // grab X from the data string
                             X = input.data[i][0]

                             yc = func(X, Par);

                             // grab Y from the data string
                             var Y = input.data[i][nVar];

                             // IN: form of std. error.
                             // expected values are:
                             //  equal: "1": all errors=1
                             //  relative: "Y": errors prop. to Y values, good for exponential-type data
                             //  counts: "Sqrt(Y)": poisson noise on count data
                             //  data: "w": looks for an extra column in data with weight values
                             //  replicates: "Rep": a bunch of Y values are on the line, get spread of these values
                             vSEy = input.errorform;
                             if(vSEy===undefined) { vSEy='1'; }

                             if(vSEy=="w") { w = input.data[i][nVar+1]; }
                             if(vSEy=="Rep") {
                                 nY = 1;
                                 sY = Y;
                                 sYY = Y*Y;
                                 for(j=nVar+1; j<input.data[i].length; j++){
                                     Y = input.data[i][j];
                                 nY = nY + 1; sY = sY + Y; sYY = sYY + Y*Y;
                                 }
                                 Y = sY/nY;
                                 Rep = SQRT(ABS(sYY/nY - Y*Y)/nY);
                             }
                             var yo = Y;
                             w = eval(vSEy);
                             if(w===0) { w = 0.001; }



                             // IN: least absolute fitting, rather than least squares
                             if(input.leastabs) {
                                 w *= MAX( SQRT(ABS(Y-yc)*((Y<yc) ? Pctile : (1-Pctile))), 0.001*Y );
                             }
                             cccSW += 1 / (w * w);
                             cccAv += yc / (w * w);
                             cccSD += (Y - ccAv) * (Y - ccAv) / (w * w);

                             // Calculate the derivatives
                             for (j=0; j<nPar; j++) {
                                 Save = Par[j];
                                 Del = (Save===0) ? 0.0001 : Save/1000;
                                 Par[j] = Save + Del;

                                 var ycDecr, ycIncr = func(X, Par);

                                 if(input.centered) {
                                     Par[j] = Save - Del;
                                     ycDecr = func(X, Par);
                                     Der[j] = (ycIncr - ycDecr) / (2 * Del * w);
                                 }
                                 else Der[j] = (ycIncr - yc) / (Del * w);
                                 Par[j] = Save;
                             }
                             Der[nPar] = (Y - yc) / w;

                             SSq += Der[nPar]*Der[nPar];
                             var SEest = 0;
                             for (j=0; j<nPar; j++) {
                                 for (k=0; k<=nPar; k++)
                                     Arr[ix(j,k,nPar)] += Der[j] * Der[k];
                                 SEest += Cov[ix(j,j,nPar)] * Der[j] * Der[j];
                                 for (k=j+1; k<nPar; k++)
                                     SEest += 2 * Cov[ix(j,k,nPar)] * Der[j] * Der[k];
                             }
                             SEest=w*SQRT(SEest);
                             var yco=yc, ycl=yc-St95*SEest, ych=yc+St95*SEest;

                             out.fit.push(yco);
                             out.res.push(yo-yco);
                             out.conf_lo.push(ycl);
                             out.conf_hi.push(ych);
                         }

                         ccSW = cccSW;
                         ccAv = cccAv / ccSW;
                         ccSD = cccSD / ccSW;
                         GenR2 = (ccSD-(SSq/ccSW))/ccSD;
                         GenR = SQRT(GenR2);
                                 if (log)
                                     console.log("Corr. Coeff. = " + vFmt(GenR) + "; r*r = " + vFmt(GenR2));

                         var RMS = SQRT(SSq/MAX(1,dgfr));
                             if (log)
                                 console.log("RMS Error = " + vFmt(RMS) + "; d.f = " + dgfr + "; SSq = " + vFmt(SSq));

                         var AICc, AIC = nPts * LN(SSq / nPts) + 2*(nPar+1);
                         if(nPts>=(nPar+2)) { AICc = AIC + (2 * (nPar+1) * ((nPar+1) + 1)) / (nPts - (nPar+1) - 1); }
                         else { AICc = AIC; }
                         if (log)
                             console.log("AIC = " + vFmt(AIC) + "; AIC(corrected) = " + vFmt(AICc));

                         for (i=0; i<nPar; i++) {
                             var s = Arr[ix(i,i,nPar)];
                             Arr[ix(i,i,nPar)] = 1;
                             for (k=0; k<=nPar; k++) { Arr[ix(i,k,nPar)] = Arr[ix(i,k,nPar)] / s; }
                             for (j=0; j<nPar; j++) {
                                 if (i!=j) {
                                     s = Arr[ix(j,i,nPar)];
                                     Arr[ix(j,i,nPar)] = 0;
                                     for (k=0; k<=nPar; k++) {
                                         Arr[ix(j,k,nPar)] -= s * Arr[ix(i,k,nPar)];
                                     }
                                 }
                             }
                         }

                         // IN: relaxation factor (set to <1 for more robust but slower convergence
                         FRelax = cropNum(input.frelax,0,1,1);
                         if (log)
                             console.log(NL+"Parameter Estimates...");
                         for( i=0; i<nPar; i++) {
                             Par[i] = Par[i] + FRelax * Arr[ix(i,nPar,nPar)];
                             SEP[i] = RMS * SQRT(Arr[ix(i,i,nPar)]);
                             if (log)
                                 console.log("p"+(i+1)+"="+vFmt(Par[i])+" +/- "+vFmt(SEP[i])+"; p="+Fmt(STUDT(Par[i]/SEP[i],dgfr)));
                         }
                         if (log)
                             console.log(NL+"Covariance Matrix Terms and Error-Correlations...");
                         var v;
                         for (j=0; j<nPar; j++) {
                             for (k=j; k<nPar; k++) {
                                 Cov[ix(j,k,nPar)] = Arr[ix(j,k,nPar)] * RMS * RMS;
                                 v = Arr[ix(j,k,nPar)]/SQRT(Arr[ix(j,j,nPar)]*Arr[ix(k,k,nPar)]);
                                 var o = ("B(" + (j+1) + "," + (k+1) + ")=");
                                 if(j!=k) o += ("B(" + (k+1) + "," + (j+1) + ")=");
                                 v = Fmt(v);
                                 o += " " + Cov[ix(j,k,nPar)] + "; r=" + v.substring(v.length-7,v.length);
                                 if (log)
                                     console.log(o);
                             }
                         }

            for(i=0;i<8;i++) out.params[i]={'val':Par[i], 'err':SEP[i], 'pval':vFmt(STUDT(Par[i]/SEP[i],dgfr))};
            out.rms=RMS;
            out.cov=Cov;

            return out;
        }


        function makefunc (funct) {

            var func
              , parinds = {
                  "A": "P[0]"
                , "B": "P[1]"
                , "C": "P[2]"
                , "D": "P[3]"
                , "E": "P[4]"
                , "F": "P[5]"
                , "G": "P[6]"
                , "H": "P[7]"
                , "T": "X"
              }
              , fstr;

            funct = funct.toUpperCase()

            Object.keys(parinds).forEach( function (key) {
                var reg = new RegExp("(^|[^A-Z])"+key+"([^A-Z]|$)", "g")
                funct = funct.replace(reg, "$1"+parinds[key]+"$2")
            })

            fstr = "function func (X, P) {" +
                "return " + funct + ";" +
                "}";

            eval(fstr)

            return func;
        }




        self.Iterate = Iterate;
        self.makefunc = makefunc;

        return self
    }
}