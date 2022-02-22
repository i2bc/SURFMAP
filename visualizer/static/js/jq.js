$(function () {
    
    $('#svg-container').on('click', function(event) {
        if (event.target.tagName == 'rect') {
            console.log(event.target);
        }
    });

});