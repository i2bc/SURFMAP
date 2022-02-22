// document.addEventListener( "DOMContentLoaded", function() {
//     var stage = new NGL.Stage( "viewport" );
//     stage.loadFile( "1g3n_A.pdb", { defaultRepresentation: true } );
// } );

var stage = new NGL.Stage( "viewport" );

$(function () {
    
    stage.loadFile("1g3n_A.pdb").then(function (o) {
        o.setSelection();
        o.addRepresentation("surface", {colorScheme: "bfactor"});
        o.autoView();
        console.log(o);
      });

      stage.setParameters({
        backgroundColor: "white"
      })
    
});

