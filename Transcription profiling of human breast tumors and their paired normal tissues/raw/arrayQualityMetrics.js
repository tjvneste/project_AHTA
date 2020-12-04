// (C) Wolfgang Huber 2010-2011

// Script parameters - these are set up by R in the function 'writeReport' when copying the 
//   template for this script from arrayQualityMetrics/inst/scripts into the report.

var highlightInitial = [ false, true, true, true, false, true, false, true, true, false, false, false, false, false, false, true, false, false, true, false, false, true, true, true, true, false, true, true, true, true, true, true, true, true, true, true, false, true, true, true, true, true, true, true, true, true, true, true, true, true, false, false, false, false, true, true, true, true, true, true, true, true, true, false, true, true, true, true, true, false, false, false, false, false, true, true, true, false, true, true, true, true, false, true, false, true ];
var arrayMetadata    = [ [ "1", "GSM398074.CEL", "1", "05/23/06 14:12:04" ], [ "2", "GSM398075.CEL", "2", "05/23/06 14:22:12" ], [ "3", "GSM398076.CEL", "3", "02/15/06 11:59:17" ], [ "4", "GSM398077.CEL", "4", "02/15/06 12:09:47" ], [ "5", "GSM398078.CEL", "5", "02/15/06 12:20:06" ], [ "6", "GSM398079.CEL", "6", "02/15/06 12:31:19" ], [ "7", "GSM398080.CEL", "7", "03/24/06 15:52:39" ], [ "8", "GSM398081.CEL", "8", "03/24/06 16:02:56" ], [ "9", "GSM398082.CEL", "9", "03/24/06 12:27:22" ], [ "10", "GSM398083.CEL", "10", "03/24/06 12:37:30" ], [ "11", "GSM398084.CEL", "11", "03/30/06 11:33:19" ], [ "12", "GSM398085.CEL", "12", "03/30/06 11:43:29" ], [ "13", "GSM398086.CEL", "13", "03/24/06 12:47:52" ], [ "14", "GSM398087.CEL", "14", "03/24/06 13:04:03" ], [ "15", "GSM398088.CEL", "15", "08/08/06 15:23:28" ], [ "16", "GSM398089.CEL", "16", "08/08/06 15:34:13" ], [ "17", "GSM398090.CEL", "17", "05/23/06 14:32:30" ], [ "18", "GSM398091.CEL", "18", "05/23/06 13:07:46" ], [ "19", "GSM398092.CEL", "19", "05/23/06 13:18:03" ], [ "20", "GSM398093.CEL", "20", "05/23/06 13:29:48" ], [ "21", "GSM398094.CEL", "21", "05/23/06 11:51:31" ], [ "22", "GSM398095.CEL", "22", "05/23/06 12:02:31" ], [ "23", "GSM398096.CEL", "23", "10/12/06 13:10:04" ], [ "24", "GSM398097.CEL", "24", "10/12/06 12:10:33" ], [ "25", "GSM398098.CEL", "25", "11/22/06 13:17:54" ], [ "26", "GSM398099.CEL", "26", "10/12/06 13:00:02" ], [ "27", "GSM398100.CEL", "27", "11/22/06 13:38:28" ], [ "28", "GSM398101.CEL", "28", "10/12/06 12:20:49" ], [ "29", "GSM398102.CEL", "29", "10/12/06 12:00:22" ], [ "30", "GSM398103.CEL", "30", "11/22/06 13:28:17" ], [ "31", "GSM398104.CEL", "31", "04/03/07 12:42:35" ], [ "32", "GSM398105.CEL", "32", "11/22/06 15:41:48" ], [ "33", "GSM398106.CEL", "33", "04/03/07 12:53:38" ], [ "34", "GSM398107.CEL", "34", "12/28/06 14:53:31" ], [ "35", "GSM398108.CEL", "35", "01/23/07 12:38:12" ], [ "36", "GSM398109.CEL", "36", "01/23/07 13:53:45" ], [ "37", "GSM398110.CEL", "37", "04/03/07 13:05:21" ], [ "38", "GSM398111.CEL", "38", "11/22/06 15:52:04" ], [ "39", "GSM398112.CEL", "39", "12/28/06 15:17:48" ], [ "40", "GSM398113.CEL", "40", "12/28/06 16:06:04" ], [ "41", "GSM398114.CEL", "41", "12/28/06 15:28:22" ], [ "42", "GSM398115.CEL", "42", "12/28/06 15:48:25" ], [ "43", "GSM398116.CEL", "43", "01/23/07 12:49:09" ], [ "44", "GSM398117.CEL", "44", "01/23/07 12:59:22" ], [ "45", "GSM398118.CEL", "45", "01/23/07 14:51:19" ], [ "46", "GSM398119.CEL", "46", "01/23/07 15:01:35" ], [ "47", "GSM398120.CEL", "47", "12/28/06 12:32:04" ], [ "48", "GSM398121.CEL", "48", "12/28/06 12:43:08" ], [ "49", "GSM398122.CEL", "49", "12/28/06 12:56:27" ], [ "50", "GSM398123.CEL", "50", "12/28/06 13:06:41" ], [ "51", "GSM398124.CEL", "51", "03/06/07 14:35:15" ], [ "52", "GSM398125.CEL", "52", "03/06/07 13:03:41" ], [ "53", "GSM398126.CEL", "53", "03/06/07 15:44:38" ], [ "54", "GSM398127.CEL", "54", "03/06/07 13:35:15" ], [ "55", "GSM398128.CEL", "55", "03/06/07 13:14:20" ], [ "56", "GSM398129.CEL", "56", "03/06/07 15:23:46" ], [ "57", "GSM398130.CEL", "57", "03/06/07 13:24:50" ], [ "58", "GSM398131.CEL", "58", "04/03/07 13:15:38" ], [ "59", "GSM398132.CEL", "59", "04/03/07 14:31:39" ], [ "60", "GSM398133.CEL", "60", "04/03/07 14:20:50" ], [ "61", "GSM398134.CEL", "61", "04/03/07 16:35:57" ], [ "62", "GSM398135.CEL", "62", "04/03/07 16:56:30" ], [ "63", "GSM398136.CEL", "63", "04/03/07 15:22:30" ], [ "64", "GSM398137.CEL", "64", "04/03/07 15:36:12" ], [ "65", "GSM398138.CEL", "65", "04/17/07 14:17:52" ], [ "66", "GSM398139.CEL", "66", "04/17/07 14:07:38" ], [ "67", "GSM398140.CEL", "67", "04/03/07 15:12:08" ], [ "68", "GSM398141.CEL", "68", "04/03/07 16:46:13" ], [ "69", "GSM398142.CEL", "69", "04/17/07 12:11:43" ], [ "70", "GSM398143.CEL", "70", "04/17/07 12:22:13" ], [ "71", "GSM398144.CEL", "71", "04/17/07 12:32:33" ], [ "72", "GSM398145.CEL", "72", "04/17/07 12:42:59" ], [ "73", "GSM398146.CEL", "73", "04/17/07 14:28:13" ], [ "74", "GSM398147.CEL", "74", "04/17/07 13:57:03" ], [ "75", "GSM398148.CEL", "75", "04/17/07 16:12:57" ], [ "76", "GSM398149.CEL", "76", "04/17/07 14:39:01" ], [ "77", "GSM398150.CEL", "77", "04/17/07 16:02:40" ], [ "78", "GSM398151.CEL", "78", "04/17/07 15:00:14" ], [ "79", "GSM398152.CEL", "79", "05/24/07 12:35:04" ], [ "80", "GSM398153.CEL", "80", "05/24/07 12:45:37" ], [ "81", "GSM398154.CEL", "81", "05/24/07 13:06:17" ], [ "82", "GSM398155.CEL", "82", "05/24/07 13:43:45" ], [ "83", "GSM398156.CEL", "83", "05/24/07 14:05:48" ], [ "84", "GSM398157.CEL", "84", "05/24/07 14:54:48" ], [ "85", "GSM398158.CEL", "85", "05/24/07 14:16:19" ], [ "86", "GSM398159.CEL", "86", "05/24/07 15:16:51" ] ];
var svgObjectNames   = [ "pca", "dens", "dig" ];

var cssText = ["stroke-width:1; stroke-opacity:0.4",
               "stroke-width:3; stroke-opacity:1" ];

// Global variables - these are set up below by 'reportinit'
var tables;             // array of all the associated ('tooltips') tables on the page
var checkboxes;         // the checkboxes
var ssrules;


function reportinit() 
{
 
    var a, i, status;

    /*--------find checkboxes and set them to start values------*/
    checkboxes = document.getElementsByName("ReportObjectCheckBoxes");
    if(checkboxes.length != highlightInitial.length)
	throw new Error("checkboxes.length=" + checkboxes.length + "  !=  "
                        + " highlightInitial.length="+ highlightInitial.length);
    
    /*--------find associated tables and cache their locations------*/
    tables = new Array(svgObjectNames.length);
    for(i=0; i<tables.length; i++) 
    {
        tables[i] = safeGetElementById("Tab:"+svgObjectNames[i]);
    }

    /*------- style sheet rules ---------*/
    var ss = document.styleSheets[0];
    ssrules = ss.cssRules ? ss.cssRules : ss.rules; 

    /*------- checkboxes[a] is (expected to be) of class HTMLInputElement ---*/
    for(a=0; a<checkboxes.length; a++)
    {
	checkboxes[a].checked = highlightInitial[a];
        status = checkboxes[a].checked; 
        setReportObj(a+1, status, false);
    }

}


function safeGetElementById(id)
{
    res = document.getElementById(id);
    if(res == null)
        throw new Error("Id '"+ id + "' not found.");
    return(res)
}

/*------------------------------------------------------------
   Highlighting of Report Objects 
 ---------------------------------------------------------------*/
function setReportObj(reportObjId, status, doTable)
{
    var i, j, plotObjIds, selector;

    if(doTable) {
	for(i=0; i<svgObjectNames.length; i++) {
	    showTipTable(i, reportObjId);
	} 
    }

    /* This works in Chrome 10, ssrules will be null; we use getElementsByClassName and loop over them */
    if(ssrules == null) {
	elements = document.getElementsByClassName("aqm" + reportObjId); 
	for(i=0; i<elements.length; i++) {
	    elements[i].style.cssText = cssText[0+status];
	}
    } else {
    /* This works in Firefox 4 */
    for(i=0; i<ssrules.length; i++) {
        if (ssrules[i].selectorText == (".aqm" + reportObjId)) {
		ssrules[i].style.cssText = cssText[0+status];
		break;
	    }
	}
    }

}

/*------------------------------------------------------------
   Display of the Metadata Table
  ------------------------------------------------------------*/
function showTipTable(tableIndex, reportObjId)
{
    var rows = tables[tableIndex].rows;
    var a = reportObjId - 1;

    if(rows.length != arrayMetadata[a].length)
	throw new Error("rows.length=" + rows.length+"  !=  arrayMetadata[array].length=" + arrayMetadata[a].length);

    for(i=0; i<rows.length; i++) 
 	rows[i].cells[1].innerHTML = arrayMetadata[a][i];
}

function hideTipTable(tableIndex)
{
    var rows = tables[tableIndex].rows;

    for(i=0; i<rows.length; i++) 
 	rows[i].cells[1].innerHTML = "";
}


/*------------------------------------------------------------
  From module 'name' (e.g. 'density'), find numeric index in the 
  'svgObjectNames' array.
  ------------------------------------------------------------*/
function getIndexFromName(name) 
{
    var i;
    for(i=0; i<svgObjectNames.length; i++)
        if(svgObjectNames[i] == name)
	    return i;

    throw new Error("Did not find '" + name + "'.");
}


/*------------------------------------------------------------
  SVG plot object callbacks
  ------------------------------------------------------------*/
function plotObjRespond(what, reportObjId, name)
{

    var a, i, status;

    switch(what) {
    case "show":
	i = getIndexFromName(name);
	showTipTable(i, reportObjId);
	break;
    case "hide":
	i = getIndexFromName(name);
	hideTipTable(i);
	break;
    case "click":
        a = reportObjId - 1;
	status = !checkboxes[a].checked;
	checkboxes[a].checked = status;
	setReportObj(reportObjId, status, true);
	break;
    default:
	throw new Error("Invalid 'what': "+what)
    }
}

/*------------------------------------------------------------
  checkboxes 'onchange' event
------------------------------------------------------------*/
function checkboxEvent(reportObjId)
{
    var a = reportObjId - 1;
    var status = checkboxes[a].checked;
    setReportObj(reportObjId, status, true);
}


/*------------------------------------------------------------
  toggle visibility
------------------------------------------------------------*/
function toggle(id){
  var head = safeGetElementById(id + "-h");
  var body = safeGetElementById(id + "-b");
  var hdtxt = head.innerHTML;
  var dsp;
  switch(body.style.display){
    case 'none':
      dsp = 'block';
      hdtxt = '-' + hdtxt.substr(1);
      break;
    case 'block':
      dsp = 'none';
      hdtxt = '+' + hdtxt.substr(1);
      break;
  }  
  body.style.display = dsp;
  head.innerHTML = hdtxt;
}
