// A SCRIPT FOR MAKING JUPYTER NOTEBOOK TABLES OF CONTENTS
//
// Put this script in your directory, and then (in Julia) run a cell that has
// 
//     macro javascript_str(s) display("text/javascript", s); end
//
//     javascript"""
//     $.getScript('make_table_of_contents.js')
//     """
//
// Carlos Brody got this file from 
//
//        https://github.com/kmahelona/ipython_notebook_goodies
//
// And then modified it slightly to remove the roman numerals.
//
//   Carlos Brody   4-Sep-2016


// Converts integer to roman numeral
function romanize(num) {
    var lookup = {M:1000,CM:900,D:500,CD:400,C:100,XC:90,L:50,XL:40,X:10,IX:9,V:5,IV:4,I:1},
	roman = '',
	    i;
	for ( i in lookup ) {
	    while ( num >= lookup[i] ) {
		roman += i;
		num -= lookup[i];
	    }
	}
	return roman;
 }

// Builds a <ul> Table of Contents from all <headers> in DOM
function createTOC(){
    var toc = "";
    var level = 0;
    var levels = {}
    $('#toc').html('');

    $(":header").each(function(i){
	    if (this.id=='tocheading'){return;}
        
	    titleText = this.innerHTML;
	    openLevel = this.tagName[1];

	    if (levels[openLevel]){
		levels[openLevel] += 1;
	    } else{
		levels[openLevel] = 1;
	    }

	    if (openLevel > level) {
		toc += (new Array(openLevel - level + 1)).join('<ul class="toc">');
	    } else if (openLevel < level) {
		toc += (new Array(level - openLevel + 1)).join("</ul>");
		for (i=level;i>openLevel;i--){levels[i]=0;}
	    }

	    level = parseInt(openLevel);


	    if (this.id==''){this.id = this.innerHTML.replace(/ /g,"-")}
	    var anchor = this.id;

        // The line below adds a roman numeral at the start of each table of contents line.
	    // toc += '<li><a href="#' + anchor + '">' +  romanize(levels[openLevel].toString()) + '. ' + titleText + '</a></li>'
        // The line below just does the header as is.
        toc += '<li><a href="#' + anchor + '">' +  titleText + '</a></li>';
        
	});

    
    if (level) {
	toc += (new Array(level + 1)).join("</ul>");
    }

 
    $('#toc').append(toc);

};

// Executes the createToc function
setTimeout(function(){createTOC();},100);

// Rebuild to TOC every minute
setInterval(function(){createTOC();},60000);
