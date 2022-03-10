
var fileName = null;

$('#matrix-file').change(function () {
    file = this.files[0];

	if(file){
        fileName = file.name

		var reader = new FileReader();
		reader.onload = function(e) {
            initSurfmap(reader.result)
		}
		reader.readAsText(file);
	}

});



mapColors = {
    'stickiness': {
        'minVal': -1.273,
        'maxVal': 1.273,
        'scale': ["royalblue", "white", "darkgreen"],
    },
    'electrostatics': {
        'minVal': null,
        'maxVal': null,
        'scale': ["red", "white", "blue"],
    },
    'kyte_doolittle': {
        'minVal': -4.5,
        'maxVal': 4.5,
        'scale': ["#5f9ea0", "#faf4e0", "#8b4726"],
        // 'scale': ["#5f9ea0", "#7ac5cd", "#faf4e0", "#cd8500", "#8b4726"],
    },
    'wimley_white': {
        'minVal': 2.23,
        'maxVal': 6.1,
        'scale': ["#5F9EA0", "#faf4e0", "#8B4726"]
        // 'scale': ["#5F9EA0", "#7AC5CD", "#faf4e0", "#CD8500", "#8B4726"]
    },
    'circular_variance': {
        'minVal': 0,
        'maxVal': 1,
        'scale': ["black", "white", "blue"]
    },
    'bfactor': {
        'minVal': null,
        'maxVal': null,
        'scale': ["#CD661D", "white", "#00008B"]
    },
    'discreteBfactor': {
        'minVal': null,
        'maxVal': null,
        'scale': [
            "#FFFFFF", "#EE2C2C", "#6495ED", "#FFB90F",
            "#7A378B", "#00CD00", "#87CEFF", "#CD853F",
            "#EE82EE", "#006400", "#FF7F24", "#4D4D4D",
            "#FA8072", "#CD6889", "#27408B"
        ]
    },
}
  
var scaleColorMap;


initSurfmap = function(data) {
    $("#scale-controller").empty()
    $("svg").remove()
    $("#tooltip-table").remove()
    $("#tooltip").remove()
    
    var mapProp = null;
    Object.keys(mapColors).forEach(function(el) {
        if (fileName.indexOf(el) != -1) {
            mapProp = el
        }
    });

    if (mapProp != null) {
        scaleColorMap = mapProp;
        surfmap(data);
    }
    else {
        var labelSel = $("<label>").appendTo("#scale-controller");
        labelSel.attr('for', 'scale-select').text('Please select the property corresponding to your matrix file:')
        labelSel.attr('class', 'label-select-scale')
    
        var sel = $('<select>').appendTo("#scale-controller");
        sel.attr('id', 'scale-select')
        sel.append($("<option>").attr('value', 'empty').text('-- Select a property --'));

        Object.keys(mapColors).forEach(function(d) {
            sel.append($("<option>").attr('value',d).text(d));
        });

        // attach event listener to control
        $('#scale-select').on('change', function(e) {
            if (this.value != 'empty') {
                scaleColorMap = this.value;
                surfmap(data);
            }
        });


    }
}

var mapColorProp;

var surfmap = function(data) {
    $("#scale-controller").empty()
    $("svg").remove()
    $("#tooltip-table").remove()

    // create tooltip element and add to document body
    var tooltip = document.createElement("div")
    tooltip.id = 'tooltip'
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

    var scale_factor = 1.8;
    var width = 360*scale_factor + 80;
    var height = 180*scale_factor;
    var pixel_size = 5;


    //Create SVG element
    var svg = d3.select("#svg-container")
        .append("svg")
        .attr("width", width)
        .attr("height", height);
        
    mapColorProp = mapColors[`${scaleColorMap}`];

    dataset.then(function(data) {
        if (mapColorProp.minVal != null) {
            var minVal = mapColorProp.minVal
            var maxVal = mapColorProp.maxVal
        }
        else {
            var minVal = d3.min(data.map(e => parseFloat(e.svalue)))
            var maxVal = d3.max(data.map(function(ele) {
                if (ele.svalue != "Inf") {
                    return parseFloat(ele.svalue)
                }
            }));
            mapColorProp.minVal = minVal;
            mapColorProp.maxVal = maxVal;
        }
    
        var colors = d3.scaleLinear()
            .domain([minVal, (minVal+maxVal)/2, maxVal])
            .range(mapColorProp.scale);
            

        var rects = svg.selectAll(".rects")
            .data(data)
            .enter()
            .append("rect")
            .attr("y", d => (180-d.ord-pixel_size)*scale_factor)
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

        var legendMargin = { top: 5, bottom: 5, left: 5, right: 20 };

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
            .attr('stop-color', mapColorProp.scale[0])
            .attr('stop-opacity', 1)
        gradient.append('stop')
            .attr('offset', '50%')
            .attr('stop-color', mapColorProp.scale[1])
            .attr('stop-opacity', 1)
        gradient.append('stop')
            .attr('offset', '100%')
            .attr('stop-color', mapColorProp.scale[2])
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
            .domain([mapColorProp.minVal, mapColorProp.maxVal])
            .range([height, 0]);

        var legendAxis = d3.axisRight(legendScale)                

        legendSvg.append("g")
            .attr("class", "legend axis")
            .attr("transform", "translate(" + legendWidth + ", 0)")
            .call(legendAxis);            

        rects.on('mouseover', function(event) {
            residues = this.__data__.residues.split(',')
            if (residues.length == 1 && residues[0] == 'NA') return

            value = this.__data__.svalue

            tooltip_table = `
            <table id='tooltip-table'>
                <thead>
                    <tr>
                        <th>resname</th>
                        <th>resnb</th>
                        <th>chain</th>
                        <th>value</th>
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
                    <td>${value}</td>
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


            this.classList.forEach(function(d) {
                className = `.${d}`
                sequenceRes = $(`table td${className}`)
                // console.log(sequenceRes)

                sequenceRes.toggleClass('hovered')

        })

        });

    });
};
