
$('#matrix-file').change(function () {
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


mapColors = {
    'stickiness': {
        'minVal': -1.273,
        'maxVal': 1.273,
        'scale': ["royalblue", "white", "darkgreen"]
    }
}

var surfmap = function(data) {

    d3.select("svg").remove()

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

    // console.log(dataset)

    var scale_factor = 1.8;
    var width = 360*scale_factor + 80;
    var height = 180*scale_factor;
    var pixel_size = 5;

    console.log(width);

    //Create SVG element
    var svg = d3.select("#svg-container")
        .append("svg")
        .attr("width", width)
        .attr("height", height);

    var minVal = mapColors.stickiness.minVal
    var maxVal = mapColors.stickiness.maxVal

    var colors = d3.scaleLinear()
        .domain([minVal, 0, maxVal])
        .range(mapColors.stickiness.scale);


    dataset.then(function(data) {

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
            .attr("stroke", "gray")
            .attr("stroke-width", 0.25)
            .attr("class", function(d) {
                residues = d.residues.split(',')
                if (residues.length == 1 && residues[0] == 'NA') return

                className = []
                residues.forEach(function(d) {
                    d = d.trim();
                    resname = d.split('_')[0];
                    resnb = d.split('_')[1];
                    chain = d.split('_')[2];

                    className.push(`res-${resname}_${resnb}_${chain}`)
                });

                return className.join(' ')
                });

        // add the legend now
        var legendFullHeight = height;
        var legendFullWidth = 50;

        var legendMargin = { top: 20, bottom: 20, left: 5, right: 20 };

        // use same margins as main plot
        var legendWidth = legendFullWidth - legendMargin.left - legendMargin.right;
        var legendHeight = legendFullHeight - legendMargin.top - legendMargin.bottom;

        var legendSvg = svg.append("g")
            .attr("transform", "translate(" + 670 + "," + 0 + ")")

        // append gradient bar
        var gradient = legendSvg.append('defs')
            .append('linearGradient')
            .attr('id', 'gradient')
            .attr('x1', '0%') // bottom
            .attr('y1', '100%')
            .attr('x2', '0%') // to top
            .attr('y2', '0%')


        gradient.append('stop')
            .attr('offset', '0%')
            .attr('stop-color', mapColors.stickiness.scale[0])
            .attr('stop-opacity', 1)
        gradient.append('stop')
            .attr('offset', '50%')
            .attr('stop-color', mapColors.stickiness.scale[1])
            .attr('stop-opacity', 1)
        gradient.append('stop')
            .attr('offset', '100%')
            .attr('stop-color', mapColors.stickiness.scale[2])
            .attr('stop-opacity', 1)

        legendSvg.append('rect')
            .attr('x1', 0)
            .attr('y1', 0)
            .attr('width', legendWidth)
            .attr('height', height)
            .style('fill', 'url(#gradient)')
            .attr("stroke", "black")
            .attr("stroke-width", 0.15)

        // create a scale and axis for the legend
        var legendScale = d3.scaleLinear()
            .domain([mapColors.stickiness.minVal, mapColors.stickiness.maxVal])
            .range([height, 0]);

        // var colors = d3.scaleLinear()
        //     .domain([minVal, 0, maxVal])
        //     .range(mapColors.stickiness.scale);

        var legendAxis = d3.axisRight(legendScale)
                

        // var legendAxis = d3.svg.axis()
        //     .scale(legendScale)
        //     .orient("right")
        //     .tickValues(d3.range(mapColors.stickiness.minVal, mapColors.stickiness.maxVal))
        //     .tickFormat(d3.format("d"));

        legendSvg.append("g")
            .attr("class", "legend axis")
            .attr("transform", "translate(" + legendWidth + ", 0)")
            .call(legendAxis);            

        rects.on('mouseover', function(event) {
            residues = this.__data__.residues.split(',')
            if (residues.length == 1 && residues[0] == 'NA') return

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
        });

        rects.on('click', function(event) {
            residues = this.__data__.residues.split(',')
            if (residues.length == 1 && residues[0] == 'NA') return

            // resIndex = residues.split('_')[1]
            // console.log(residues)
            console.log(this)

            this.classList.forEach(function(d) {
                className = `.${d}`
                sequenceRes = $(`table td${className}`)
                // console.log(sequenceRes)

                sequenceRes.toggleClass('hovered')

        })

        });

    });
};
