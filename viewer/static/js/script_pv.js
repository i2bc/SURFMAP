// override the default options with something less restrictive.
var options = {
    width: 600,
    height: 600,
    antialias: true,
    quality : 'medium'
  };

var viewer = pv.Viewer(document.getElementById('viewer'), options);
  
viewer.on('viewerReady', function() {

    pv.io.fetchPdb('1ake.pdb', function(structure) {
        viewer.cartoon('protein', structure);
        viewer.autoZoom();
    
    });


});
