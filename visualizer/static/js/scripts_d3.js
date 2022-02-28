
$('input[type=file]').change(function () {
    file = this.files[0];

	if(file){
		var reader = new FileReader();
		reader.onload = function(e) {
            surfmap(reader.result)
		}
		reader.readAsText(file);
	}

});

var _surfmap = function(data) {
    console.log(data);
}


var surfmap = function(data) {

    // create tooltip element and add to document body
    var tooltip = document.createElement("div")
    Object.assign(tooltip.style, {
        display: "none",
        position: "fixed",
        zIndex: 10,
        pointerEvents: "none",
        backgroundColor: "rgba( 0, 0, 0, 0.6 )",
        color: "lightgrey",
        padding: "8px",
        fontFamily: "sans-serif"
    })
    document.body.appendChild(tooltip)


    dataset = new Promise(( resolve, reject ) => {
        resolve(d3.tsvParse(data));
    });

    console.log(dataset)

    var scale_factor = 2;
    var width = 360*scale_factor;
    var height = 180*scale_factor;
    var pixel_size = 5; 

    //Create SVG element
    var svg = d3.select("#svg-container ")
        .append("svg")
        .attr("width", width)
        .attr("height", height);

    dataset.then(function(data) {
        console.log(data)
        var minVal = d3.min(data.map(e => parseFloat(e.svalue)))
        var maxVal = d3.max(data.map(function(ele) {
            if (ele.svalue != "Inf") {
                return parseFloat(ele.svalue)
            }
        }));

        var palette_color = ["#5E4FA2", "#3288BD", "#66C2A5", "#ABDDA4", "#E6F598", "#FFFFBF", "#FEE08B", "#FDAE61", "#F46D43", "#D53E4F", "#9E0142"]

        // var colors = d3.scaleQuantize()
        //     .domain([minVal ,maxVal])
        //     .range(palette_color);

        var colors = d3.scaleLinear()
            .domain([minVal, 0, maxVal])
            .range(["blue", "white", "red"]);


        var rects = svg.selectAll(".rects")
            .data(data)
            .enter()
            .append("rect")
            .attr("y", d => (d.ord-pixel_size)*scale_factor)
            .attr("x", d => (d.absc-pixel_size)*scale_factor)
            .attr("height", pixel_size*scale_factor)
            .attr("width", pixel_size*scale_factor)
            .attr("fill", function(d) {
                if (d.svalue != 'Inf') return colors(d.svalue)
                else return 'white';
            })
            .attr("stroke", "gray");



        rects.append("title")
            .text(function (d) {
                if (d.residues != "NA") {
                    resname = d.residues.split('_')[0];
                    resnb = d.residues.split('_')[1];
                    chain = d.residues.split('_')[2];
                    return `${resname} ${resnb} ${chain}`
                }
            });

        rects.on('click', function(event) {
            console.log(this.__data__);

            residues = this.__data__.residues.split(',')
            if (residues.length == 1 && residues[0] == 'NA') { return }

            tooltip_table = `
            <table id='tooltip'>
                <thead>
                    <tr>
                        <th>resname</th>
                        <th>resnb</th>
                        <th>chain</th>
                    </tr>                    
                </thead>
                <tbody>
            `
            table_tr = '';

            residues.forEach(function(d) {
                resname = d.split('_')[0];
                resnb = d.split('_')[1];
                chain = d.split('_')[2];

                table_tr = table_tr + `
                <tr>
                    <td>${resname}</td>
                    <td>${resnb}</td>
                    <td>${chain}</td>
                </tr>
                `
            });

            tooltip_table = tooltip_table + table_tr + `                  
                </tbody>
            </table>
            `

            tooltip.innerHTML = tooltip_table

            tooltip.style.bottom = window.innerHeight - event.clientY + 3 + "px"
            tooltip.style.left = event.clientX + 3 + "px"
            if (tooltip.style.display == 'none') {
                tooltip.style.display = "block" 
            }
            else {
                tooltip.style.display = "none" 
            }


        })
        // console.log(minVal, maxVal)
    });
};
