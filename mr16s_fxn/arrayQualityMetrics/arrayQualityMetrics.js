// (C) Wolfgang Huber 2010-2011

// Script parameters - these are set up by R in the function 'writeReport' when copying the 
//   template for this script from arrayQualityMetrics/inst/scripts into the report.

var highlightInitial = [ false, false, true, false, false, false, false, false, false, false, true, false, false, false, false, false, false, true, false, true, false, true, false, false, false, false, false, false, false, false, false, false, false, false, false, false, true, true, false, false, false, false, false, false, false, true, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, true, true, false, false, false, false, false, false ];
var arrayMetadata    = [ [ "1", "A1", "0.718330549031588", "A", "1", "Moorea NE Inner", "Moorea", "MNW", "in", "MNWI" ], [ "2", "A10", "1.1219273970582", "A", "10", "Moorea SE Inner", "Moorea", "MSE", "in", "MSEI" ], [ "3", "A11", "0.421688910120062", "A", "11", "Moorea SE Outer", "Moorea", "MSE", "out", "MSEO" ], [ "4", "A12", "0.521649165373198", "A", "12", "Moorea SE Outer", "Moorea", "MSE", "out", "MSEO" ], [ "5", "A2", "1.07328100543593", "A", "2", "Moorea NE Inner", "Moorea", "MNW", "in", "MNWI" ], [ "6", "A3", "0.776915913402586", "A", "3", "Moorea NE Outer", "Moorea", "MNW", "out", "MNWO" ], [ "7", "A5", "1.19529699899849", "A", "5", "Tahiti Inner", "Tahiti", "TNW", "in", "TI" ], [ "8", "A6", "1.743088419951", "A", "6", "Tahiti Inner", "Tahiti", "TNW", "in", "TI" ], [ "9", "A7", "0.780187691505056", "A", "7", "Tahiti Outer", "Tahiti", "TNW", "out", "TO" ], [ "10", "A9", "0.147229063146162", "A", "9", "Moorea SE Inner", "Moorea", "MSE", "in", "MSEI" ], [ "11", "B1", "3.30530915578713", "B", "1", "Moorea NE Inner", "Moorea", "MNW", "in", "MNWI" ], [ "12", "B10", "1.81634141093491", "B", "10", "Moorea SE Inner", "Moorea", "MSE", "in", "MSEI" ], [ "13", "B11", "0.776309555820615", "B", "11", "Moorea SE Outer", "Moorea", "MSE", "out", "MSEO" ], [ "14", "B12", "0.266937153786741", "B", "12", "Moorea SE Outer", "Moorea", "MSE", "out", "MSEO" ], [ "15", "B2", "2.22503814585589", "B", "2", "Moorea NE Inner", "Moorea", "MNW", "in", "MNWI" ], [ "16", "B3", "1.1916351531922", "B", "3", "Moorea NE Outer", "Moorea", "MNW", "out", "MNWO" ], [ "17", "B4", "1.44599927773479", "B", "4", "Moorea NE Outer", "Moorea", "MNW", "out", "MNWO" ], [ "18", "B6", "1.54710134963843", "B", "6", "Tahiti Inner", "Tahiti", "TNW", "in", "TI" ], [ "19", "B7", "2.3313746107888", "B", "7", "Tahiti Outer", "Tahiti", "TNW", "out", "TO" ], [ "20", "B9", "3.3035973669727", "B", "9", "Moorea SE Inner", "Moorea", "MSE", "in", "MSEI" ], [ "21", "C1", "0.228829585554471", "C", "1", "Moorea NE Inner", "Moorea", "MNW", "in", "MNWI" ], [ "22", "C10", "1.09529147134485", "C", "10", "Moorea SE Inner", "Moorea", "MSE", "in", "MSEI" ], [ "23", "C11", "0.635326235925568", "C", "11", "Moorea SE Outer", "Moorea", "MSE", "out", "MSEO" ], [ "24", "C3", "0.823437774503989", "C", "3", "Moorea NE Outer", "Moorea", "MNW", "out", "MNWO" ], [ "25", "C4", "1.52635769848673", "C", "4", "Moorea NE Outer", "Moorea", "MNW", "out", "MNWO" ], [ "26", "C5", "1.74131705007696", "C", "5", "Tahiti Inner", "Tahiti", "TNW", "in", "TI" ], [ "27", "C6", "0.338022208971007", "C", "6", "Tahiti Inner", "Tahiti", "TNW", "in", "TI" ], [ "28", "C7", "1.62513575551175", "C", "7", "Tahiti Outer", "Tahiti", "TNW", "out", "TO" ], [ "29", "C8", "1.57317956269512", "C", "8", "Tahiti Outer", "Tahiti", "TNW", "out", "TO" ], [ "30", "C9", "0.774553234841931", "C", "9", "Moorea SE Inner", "Moorea", "MSE", "in", "MSEI" ], [ "31", "D1", "2.18901114189168", "D", "1", "Moorea NE Inner", "Moorea", "MNW", "in", "MNWI" ], [ "32", "D10", "0.868889022565817", "D", "10", "Moorea SE Inner", "Moorea", "MSE", "in", "MSEI" ], [ "33", "D11", "1.3762003963164", "D", "11", "Moorea SE Outer", "Moorea", "MSE", "out", "MSEO" ], [ "34", "D12", "1.89329709711487", "D", "12", "Moorea SE Outer", "Moorea", "MSE", "out", "MSEO" ], [ "35", "D3", "1.20135829140029", "D", "3", "Moorea NE Outer", "Moorea", "MNW", "out", "MNWO" ], [ "36", "D5", "1.6156088775713", "D", "5", "Tahiti Inner", "Tahiti", "TNW", "in", "TI" ], [ "37", "D6", "1.32047797818325", "D", "6", "Tahiti Inner", "Tahiti", "TNW", "in", "TI" ], [ "38", "D7", "3.27231205644663", "D", "7", "Tahiti Outer", "Tahiti", "TNW", "out", "TO" ], [ "39", "D8", "1.80363459444482", "D", "8", "Tahiti Outer", "Tahiti", "TNW", "out", "TO" ], [ "40", "D9", "1.37263563740889", "D", "9", "Moorea SE Inner", "Moorea", "MSE", "in", "MSEI" ], [ "41", "E1", "1.43870607837353", "E", "1", "Moorea NE Inner", "Moorea", "MNW", "in", "MNWI" ], [ "42", "E10", "0.830380452378072", "E", "10", "Moorea SE Inner", "Moorea", "MSE", "in", "MSEI" ], [ "43", "E11", "0.32267364326428", "E", "11", "Moorea SE Outer", "Moorea", "MSE", "out", "MSEO" ], [ "44", "E12", "0.673228854565288", "E", "12", "Moorea SE Outer", "Moorea", "MSE", "out", "MSEO" ], [ "45", "E2", "1.95843161583446", "E", "2", "Moorea NE Inner", "Moorea", "MNW", "in", "MNWI" ], [ "46", "E3", "1.11620814950709", "E", "3", "Moorea NE Outer", "Moorea", "MNW", "out", "MNWO" ], [ "47", "E4", "1.44696035306325", "E", "4", "Moorea NE Outer", "Moorea", "MNW", "out", "MNWO" ], [ "48", "E5", "0.395992027626821", "E", "5", "Tahiti Inner", "Tahiti", "TNW", "in", "TI" ], [ "49", "E6", "0.47285104590063", "E", "6", "Tahiti Inner", "Tahiti", "TNW", "in", "TI" ], [ "50", "E7", "1.48572745545635", "E", "7", "Tahiti Outer", "Tahiti", "TNW", "out", "TO" ], [ "51", "E9", "0.696732953817553", "E", "9", "Moorea SE Inner", "Moorea", "MSE", "in", "MSEI" ], [ "52", "F1", "2.10098353974591", "F", "1", "Moorea NE Inner", "Moorea", "MNW", "in", "MNWI" ], [ "53", "F10", "0.461632450926662", "F", "10", "Moorea SE Inner", "Moorea", "MSE", "in", "MSEI" ], [ "54", "F11", "0.255168455243733", "F", "11", "Moorea SE Outer", "Moorea", "MSE", "out", "MSEO" ], [ "55", "F12", "0.625161984441419", "F", "12", "Moorea SE Outer", "Moorea", "MSE", "out", "MSEO" ], [ "56", "F2", "2.17199941248312", "F", "2", "Moorea NE Inner", "Moorea", "MNW", "in", "MNWI" ], [ "57", "F3", "1.45812686354146", "F", "3", "Moorea NE Outer", "Moorea", "MNW", "out", "MNWO" ], [ "58", "F4", "0.617918054238025", "F", "4", "Moorea NE Outer", "Moorea", "MNW", "out", "MNWO" ], [ "59", "F5", "1.29351922579161", "F", "5", "Tahiti Inner", "Tahiti", "TNW", "in", "TI" ], [ "60", "F6", "1.41577557857935", "F", "6", "Tahiti Inner", "Tahiti", "TNW", "in", "TI" ], [ "61", "F7", "1.58061981465079", "F", "7", "Tahiti Outer", "Tahiti", "TNW", "out", "TO" ], [ "62", "F8", "1.59959909452249", "F", "8", "Tahiti Outer", "Tahiti", "TNW", "out", "TO" ], [ "63", "G1", "1.17097846308102", "G", "1", "Moorea NE Inner", "Moorea", "MNW", "in", "MNWI" ], [ "64", "G10", "0.367873101357328", "G", "10", "Moorea SE Inner", "Moorea", "MSE", "in", "MSEI" ], [ "65", "G11", "0.669578593001964", "G", "11", "Moorea SE Outer", "Moorea", "MSE", "out", "MSEO" ], [ "66", "G12", "1.02549015319176", "G", "12", "Moorea SE Outer", "Moorea", "MSE", "out", "MSEO" ], [ "67", "G2", "0.883425666205832", "G", "2", "Moorea NE Inner", "Moorea", "MNW", "in", "MNWI" ], [ "68", "G3", "0.856537448476863", "G", "3", "Moorea NE Outer", "Moorea", "MNW", "out", "MNWO" ], [ "69", "G4", "0.679872808973325", "G", "4", "Moorea NE Outer", "Moorea", "MNW", "out", "MNWO" ], [ "70", "G5", "1.41868157547912", "G", "5", "Tahiti Inner", "Tahiti", "TNW", "in", "TI" ], [ "71", "G6", "0.777646564318586", "G", "6", "Tahiti Inner", "Tahiti", "TNW", "in", "TI" ], [ "72", "G7", "1.00586792019449", "G", "7", "Tahiti Outer", "Tahiti", "TNW", "out", "TO" ], [ "73", "G8", "1.82264881456698", "G", "8", "Tahiti Outer", "Tahiti", "TNW", "out", "TO" ], [ "74", "G9", "0.508590203783608", "G", "9", "Moorea SE Inner", "Moorea", "MSE", "in", "MSEI" ], [ "75", "H1", "0.77886133903717", "H", "1", "Moorea NE Inner", "Moorea", "MNW", "in", "MNWI" ], [ "76", "H10", "1.30231327492551", "H", "10", "Moorea SE Inner", "Moorea", "MSE", "in", "MSEI" ], [ "77", "H11", "1.04989434843661", "H", "11", "Moorea SE Outer", "Moorea", "MSE", "out", "MSEO" ], [ "78", "H12", "1.29472878477699", "H", "12", "Moorea SE Outer", "Moorea", "MSE", "out", "MSEO" ], [ "79", "H2", "0.456537863408214", "H", "2", "Moorea NE Inner", "Moorea", "MNW", "in", "MNWI" ], [ "80", "H3", "0.842300542245631", "H", "3", "Moorea NE Outer", "Moorea", "MNW", "out", "MNWO" ], [ "81", "H6", "1.19042250494552", "H", "6", "Tahiti Inner", "Tahiti", "TNW", "in", "TI" ], [ "82", "H7", "0.644897092245675", "H", "7", "Tahiti Outer", "Tahiti", "TNW", "out", "TO" ], [ "83", "H8", "1.8542105348027", "H", "8", "Tahiti Outer", "Tahiti", "TNW", "out", "TO" ], [ "84", "H9", "1.35341967113754", "H", "9", "Moorea SE Inner", "Moorea", "MSE", "in", "MSEI" ] ];
var svgObjectNames   = [ "pca", "dens" ];

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
