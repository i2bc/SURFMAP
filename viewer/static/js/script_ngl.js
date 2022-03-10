// document.addEventListener( "DOMContentLoaded", function() {
//     var stage = new NGL.Stage( "viewport" );
//     stage.loadFile( "1g3n_A.pdb", { defaultRepresentation: true } );
// } );


var stage = new NGL.Stage( "viewport", { backgroundColor: "white" });
var structureComponent;
var cartoonRepr;
var stickRepr;

var isSeqClicked = false;

$('#pdb-file').change(function () {
  file = this.files[0];
   
  stage.loadFile(file).then(function (o) {
    structureComponent = o;


    structureComponent.setSelection();
    cartoonRepr = structureComponent.addRepresentation("cartoon", {
      colorScheme: "bfactor",
    });

    stickRepr = structureComponent.addRepresentation("ball+stick", {
      sele: "none",
    });
    
    structureComponent.autoView();

    createSequenceTable(structureComponent.structure);

    $('#sequence-container table tbody').on('click', function(event) {

      if ( event.target.className.indexOf('res-') == -1 ) return
  
      className = '.' + event.target.className.split(' ')[0]
      cells = $(className)
      // $("#sequence-container td").removeClass('hovered')
      cells.toggleClass('hovered')

      rects = $('rect')
      surfmapCells = $(`rect${className}`)

      if (surfmapCells.length == 0) return

      if (isSeqClicked == false) {
        rects.attr({
          "fill-opacity": 0.45,
        })

        surfmapCells.attr({
          "fill-opacity": 1,
          "stroke": "black",
          "stroke-width": 0.85
        })

        sele = className.split('_').slice(1).join(':')

        stickRepr.setSelection(sele);
        structureComponent.autoView(sele, 1000);
        isSeqClicked = true;
      }

      else {
        rects.attr({
          "fill-opacity": 1,
        })

        surfmapCells.attr({
          "stroke": "gray",
          "stroke-width": 0.25
        })
        isSeqClicked = false;

      }

  
      // console.log(event.target)
  
    });
  
    console.log(o);
  });
    
});


get2DArray = function(nglStruct, ncol) {

  /*  
  Returns a 2D array given ncol and an NGL structure object.
  The 2D arrays is focused on residues of the structure.
  
  [
    [
      res1 = {
        resno: 1,
        resname: "GLN",
        rescode: "Q",
        chain: "A"

      },
      res2 = {
        ...
      },
      ...
    ],
    ...
  ]

  
  Examples:
  For data = 'abcdefgh' and resIndexes undefined
  -> Array([['abc'], ['def'], ['gh'],])  
  */

  if ( typeof(nglStruct) === 'undefined' ) return;
  ncol = ncol ? ncol: 40

  residues = []
  nglStruct.eachResidue(function(res) {
    residue = {
      resno: res.resno,
      resname: res.resname,
      rescode: res.getResname1(),
      chain: res.chainname
    }
    residues.push(residue)
  });

  rowArray = []
  sequence2dArray = []

  for (var i=0; i<residues.length; i++) {
    if ( i%ncol == 0 && i != 0 ) {      
      sequence2dArray.push(rowArray);
      rowArray = [];
    }
    rowArray.push(residues[i]);    
  }

  if (rowArray.length > 0) {
    sequence2dArray.push(rowArray);
  }    

  return sequence2dArray
}


createSequenceTable = function(nglStruct, params) {
  if ( typeof(nglStruct) === 'undefined' ) return;

  if ( typeof(params) === 'undefined' ) {
    params = {
      ncols: 40
    }
  }

  // get protein sequence in a 2d array (1st dim will represent rows, 2nd cols)
  sequence2DArray = get2DArray(nglStruct, params.ncols);

  sequenceContainer = $("#sequence-container");

  var table = $('<table>');
  var tbody = $('<tbody>');

  sequence2DArray.forEach(function(rowData) {
    var rowIndex = $('<tr>').addClass('res-index');
    var row = $('<tr>')

    rowData.forEach(function(cellData) {
      className = `res-${cellData.resname}_${cellData.resno}_${cellData.chain}`

      cellIndex = $('<td>').text(cellData.resno);
      cellIndex.addClass(className);
      rowIndex.append(cellIndex);

      cellResname = $('<td>').text(cellData.rescode);
      cellResname.addClass(className);
      row.append(cellResname);
    });

    tbody.append(rowIndex);
    tbody.append(row);
  });

  table.append(tbody);
  sequenceContainer.append(table);

}
