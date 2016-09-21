'use strict';

angular.module('ResCtxVis.services', [
        'ngResource'
    ])

    .factory('AlignmentJob', ['$resource', function ($resource) {
        return $resource('/alignments/:alignmentRunId.json', {}, {
            get: {
                method: 'GET',
                headers: {'Accept': 'application/json'}
            },
            create: {
                method: 'POST',
                url: '/alignments/ ',
                transformRequest: function(obj) {
                        var str = [];
                        for(var p in obj) {
                            str.push(encodeURIComponent(p) + "=" + encodeURIComponent(obj[p]));
                        }
                        return str.join("&");
                    },
                headers: {
                    'Accept': 'application/json',
                    'Content-Type': 'application/x-www-form-urlencoded'
                }
            }
        });
    }])

    .factory('ContextResourceFactory', ['$resource', function ($resource) {
        return function (kind, alignmentRunId) {
            var resourceBase = '';

            if(alignmentRunId !== null) {
                resourceBase += '/alignments/' + alignmentRunId;
            }

            resourceBase += '/contexts/' + kind;

            return $resource(resourceBase + '/:identifier.json', {}, {
                get: { method: 'GET', headers: {'Accept': 'application/json'} }
            });
        };
    }])

    .factory('ChainContextFactory', ['ContextResourceFactory', function (ContextResourceFactory) {
        return function(alignmentRunId) {
            return ContextResourceFactory('sterics', alignmentRunId === undefined ? null : alignmentRunId);
        };
    }])

    .factory('GridPointContextFactory', ['ContextResourceFactory', function (ContextResourceFactory) {
        return function(kind, alignmentRunId) {
            return ContextResourceFactory(kind, alignmentRunId === undefined ? null : alignmentRunId);
        };
    }])

    .factory('Alignment', ['$resource', function($resource ) {
        return $resource('/alignments/:alignmentRunId/alignment.json', {}, {
            get: { method: 'GET', headers: {'Accept': 'application/json'} }
        });
    }])

    .factory('SparseHeatMap', [function() {
        return function (data, options) {
            var thumbnail = [],
                xOffset = 0,
                arrayIndex = 0,
                histogram = data.histogram,
                azmuthBins = options.azmuthBins || 12,
                zenithBins = options.zenithBins || 6,
                numRadialBins = options.numRadialBins || 72,
                thumbnailSize = options.thumbnailSize || 4,
                colorMap = options.colorMap || [
                        "#FFF", "#DDD", "#BBB", "#999",
                        "#777", "#555", "#555", "#333",
                        "#333", "#333", "#111", "#111",
                        "#111", "#000", "#000", "#000"
                    ],
                boundarySize = options.boundarySize || 0,
                getColor = function (val) {
                    return colorMap[Math.min(val, colorMap.length -1)];
                };
            var value, color;
            for(var i = 0; i < numRadialBins; i++) {
                for(var j = 0; j < zenithBins; j++) {
                    for(var k = 0; k < azmuthBins; k++) {
                        if(histogram[arrayIndex] !== undefined) {
                            value = histogram[arrayIndex].length || histogram[arrayIndex];
                            color = getColor(value);
                            thumbnail.push({
                                x: k * thumbnailSize + xOffset,
                                y: j * thumbnailSize,
                                width: thumbnailSize,
                                height: thumbnailSize,
                                color: color
                            });
                        }
                        arrayIndex++;
                    }
                }
                xOffset += azmuthBins * thumbnailSize + boundarySize;
            }
            return thumbnail;
        };
    }])

    .factory('DenseHeatMap', [function() {
        return function(data, options) {
            var heatmap = [],
                arrayIndex = 0,
                histogram = data.histogram,
                azmuthBins = options.azmuthBins || 12,
                zenithBins = options.zenithBins || 6,
                numRadialBins = options.numRadialBins || 5,
                colorMap = options.colorMap || [
                        "#FFF", "#DDD", "#BBB", "#999",
                        "#777", "#555", "#555", "#333",
                        "#333", "#333", "#111", "#111",
                        "#111", "#000", "#000", "#000"
                    ],
                getColor = function (val) {
                    return colorMap[Math.min(val, colorMap.length -1)];
                },
                histogram = data.histogram
            var value, color;
            for(var i = 0; i < numRadialBins; i++) {
                var radialBin = [];
                for(var j = 0; j < zenithBins; j++) {
                    var zenithBin = [];
                    for(var k = 0; k < azmuthBins; k++) {
                        var value = histogram[arrayIndex],
                            color = getColor(value);
                        zenithBin.push({
                            value: value || 0,
                            index: arrayIndex,
                            binIndex: [i, j, k],
                            color: color
                        });
                        arrayIndex++;
                    }
                    radialBin.push(zenithBin);
                }
                heatmap.push(radialBin);
            }
            return heatmap;
        };
    }])

    .factory('ResidueContextHeatMap', ['DenseHeatMap', function (DenseHeatMap) {
        return function (chainContext, index) {
            var options = {
                    oneRadialSize: chainContext.bin_data.num_bins_per_radial_unit,
                    numRadialBins: chainContext.bin_data.max_radial_bin,
                    zenithBins: chainContext.bin_data.bin_divisor,
                    azmuthBins: 2 * chainContext.bin_data.bin_divisor,
                    oneRowBins: chainContext.bin_data.num_bins / chainContext.bin_data.bin_divisor,
                    colorMap: [
                        "#FFF", "#DDD", "#BBB", "#999",
                        "#777", "#555", "#555", "#333",
                        "#333", "#333", "#111", "#111",
                        "#111", "#000", "#000", "#000"
                    ]
                },
                context = chainContext.contexts[index];
            return DenseHeatMap(context, options);
        }
    }])

    .factory('JsMolChainColoring', ['Alignment', function(Alignment) {
        return function(jsmol, $scope) {
            var ident = $scope.identifier;
            $scope.residueColoring = {};
            var setAlignmentColors = function (jsmol, alignment) {
                var command = [],
                    colors = [
                        "#B22222", "#B0FFB0", "#FFC0C8", "#FFFF80", "#FFC0FF",
                        "#B0F0F0", "#FFD070", "#F08080", "#F5DEB3", "#00BFFF",
                        "#CD5C5C", "#66CDAA", "#9ACD32", "#EE82EE", "#00CED1",
                        "#00FF7F", "#3CB371", "#00008B", "#BDB76B", "#006400",
                        "#800000", "#808000", "#800080", "#008080", "#B8860B",
                        "#C0D0FF"
                    ],
                    blocks = alignment[ident]
                            || alignment[ident.toUpperCase()]
                            || alignment[ident.toLowerCase()]
                            || [];
                for(var i = 0; i < blocks.length; i++) {
                    var color = colors[i % colors.length],
                        block = blocks[i];
                    for(var j = 0; j < block.indices.length; j++) {
                        $scope.residueColoring[block.indices[j]] = color;
                        command.push("COLOR { " + block.indices[j] + ".CA } " + '"' + color + '"');
                    }
                }
                var jmolScript = command.join(";");
                Jmol.script(jsmol, jmolScript);
                return jmolScript;
            };

            $scope.setAlignmentColors = setAlignmentColors;

            Alignment.get({alignmentRunId: $scope.alignmentRunId}, function(results) {
                $scope.alignmentData = results;
                jsmol = jsmol || $scope.jsmol;

                if(jsmol) {
                    $scope.setAlignmentColors(jsmol, results);
                }
            });
        };
    }])

    .factory('JsMolSphericalHistogram', [function () {
        var binPrefix = 'ResCtx_Bin_';
        var makeBins = function (context, options) {
            var numBins = options.numBins,
                numRadialBins = options.numRadialBins,
                vertexTemplates = options.templates.vertices,
                faces = options.templates.faces,
                frequencyColors = options.frequencyColors === undefined ? true : options.frequencyColors,
                includeEmpty = options.includeEmpty == undefined ? false: options.includeEmpty,
                getColor = function (val) {
                    return val * 10;
                },
                numRadialShells = Math.floor(numBins / numRadialBins),
                histogram = context.histogram,
                center = context.coords,
                vertices = [],
                planes = [],
                colors = [],
                lookup = {},
                empty = [];
            if(frequencyColors) {
                var scale = 0;
                for(var hI = 0; hI < numBins; hI++) {
                    scale += histogram[hI] === undefined ? 0 : histogram[hI];
                }
                scale = 10 * 255 / scale;
                getColor = function (val) {
                    return val !== undefined ? 255 - val * scale : 255;
                };
            }
            for(var vI = 0; vI < vertexTemplates.length; vI++) {
                vertices.push([
                    vertexTemplates[vI][0] + center[0],
                    vertexTemplates[vI][1] + center[1],
                    vertexTemplates[vI][2] + center[2]
                ]);
            }
            for(var hfI = 0; hfI < numBins; hfI++) {
                planes.push([
                    vertices[faces[hfI][0]],
                    vertices[faces[hfI][1]],
                    vertices[faces[hfI][2]],
                    vertices[faces[hfI][3]],
                ]);
                colors.push(getColor(histogram[hfI]));
                empty.push(histogram[hfI] == undefined || histogram[hfI] == 0);
                lookup[hfI] = planes.length- 1;
            }
            return {
                center: center,
                planes: planes,
                colors: colors,
                lookup: lookup,
                empty: empty
            };
        };

        var drawBins = function (applet, bins, options) {
            var includeEmpty = options.includeEmpty,
                oneRadialSize = options.oneRadialSize,
                script = [
                "delete $" + binPrefix + "*",
                "draw " + binPrefix + "self DIAMETER 25 POINTS [ { " + bins.center.join(' ') + " } ] COLOR red "
            ];
            for(var i = 0; i < bins.planes.length; i++) {
                var binName = binPrefix + i,
                    iPrevious = i - oneRadialSize
                if(includeEmpty || !bins.empty[i]) {
                    if(iPrevious >= 0) {
                        script.push(
                            "draw " + binName
                            + " polygon 8 "
                            + "{" + bins.planes[i][0].join(' ') + "} "
                            + "{" + bins.planes[i][1].join(' ') + "} "
                            + "{" + bins.planes[i][2].join(' ') + "} "
                            + "{" + bins.planes[i][3].join(' ') + "} "
                            + "{" + bins.planes[iPrevious][0].join(' ') + "} "
                            + "{" + bins.planes[iPrevious][1].join(' ') + "} "
                            + "{" + bins.planes[iPrevious][2].join(' ') + "} "
                            + "{" + bins.planes[iPrevious][3].join(' ') + "} "
                            + " 12 "
                            + "[ 0 1 2 5 ] [ 1 3 2 2 ] " // Front
                            + "[ 4 5 6 5 ] [ 5 7 6 2 ] " // Back
                            + "[ 0 2 4 5 ] [ 2 6 4 2 ] " // Left
                            + "[ 1 3 5 5 ] [ 3 7 5 2 ] " // Right
                            + "[ 0 1 4 5 ] [ 1 5 4 2 ] " // Bottom
                            + "[ 2 3 6 5 ] [ 3 7 6 2 ] " // Top
                            + ( bins.empty[i] ? " NOFILL " : " FILL ") + " "
                        );
                    } else {
                        script.push(
                            "draw " + binName
                            + " polygon 5 "
                            + "{" + bins.planes[i][0].join(' ') + "} "
                            + "{" + bins.planes[i][1].join(' ') + "} "
                            + "{" + bins.planes[i][2].join(' ') + "} "
                            + "{" + bins.planes[i][3].join(' ') + "} "
                            + "{" + bins.center.join(' ') + "} "
                            + " 6 " // No back needed
                            + "[ 0 1 2 5 ] [ 1 3 2 2 ] " // Front
                            + "[ 0 2 4 5 ] " // Left
                            + "[ 1 3 4 5 ] " // Right
                            + "[ 0 1 4 5 ] " // Bottom
                            + "[ 2 3 4 5 ] " // Top
                            + ( bins.empty[i] ? " NOFILL " : " FILL ") + "  "
                        );
                    }

                    script.push("COLOR $" + binName
                        + " TRANSLUCENT  .99 "
                        + "[" + bins.colors[i] + ", " + bins.colors[i] + ", " + bins.colors[i] + "]"
                    );
                }
            }
            Jmol.script(applet, script.join(';'));
        };

        var highlightBin = function (applet, bins, binIndex, clear, options) {
            var binName = binPrefix + 'highlight',
                i = bins.lookup[binIndex],
                numRadialBins = options.numRadialBins,
                oneRadialSize = options.oneRadialSize,
                iPrevious = i - oneRadialSize;
            if(clear === true) {
                Jmol.script(applet, 'DELETE $' + binName);
            } else {
                if(iPrevious >= 0) {
                    Jmol.script(applet,
                        "draw " + binName
                        + " polygon 8 "
                        + "{" + bins.planes[i][0].join(' ') + "} "
                        + "{" + bins.planes[i][1].join(' ') + "} "
                        + "{" + bins.planes[i][2].join(' ') + "} "
                        + "{" + bins.planes[i][3].join(' ') + "} "
                        + "{" + bins.planes[iPrevious][0].join(' ') + "} "
                        + "{" + bins.planes[iPrevious][1].join(' ') + "} "
                        + "{" + bins.planes[iPrevious][2].join(' ') + "} "
                        + "{" + bins.planes[iPrevious][3].join(' ') + "} "
                        + " 12 "
                        + "[ 0 1 2 5 ] [ 1 3 2 2 ] " // Front
                        + "[ 4 5 6 5 ] [ 5 7 6 2 ] " // Back
                        + "[ 0 2 4 5 ] [ 2 6 4 2 ] " // Left
                        + "[ 1 3 5 5 ] [ 3 7 5 2 ] " // Right
                        + "[ 0 1 4 5 ] [ 1 5 4 2 ] " // Bottom
                        + "[ 2 3 6 5 ] [ 3 7 6 2 ] " // Top
                        + " FILL MESH "
                    );
                } else {
                    Jmol.script(applet,
                        "draw " + binName
                        + " polygon 5 "
                        + "{" + bins.planes[i][0].join(' ') + "} "
                        + "{" + bins.planes[i][1].join(' ') + "} "
                        + "{" + bins.planes[i][2].join(' ') + "} "
                        + "{" + bins.planes[i][3].join(' ') + "} "
                        + "{" + bins.center.join(' ') + "} "
                        + " 6 " // No back needed
                        + "[ 0 1 2 5 ] [ 1 3 2 2 ] " // Front
                        + "[ 0 2 4 5 ] " // Left
                        + "[ 1 3 4 5 ] " // Right
                        + "[ 0 1 4 5 ] " // Bottom
                        + "[ 2 3 4 5 ] " // Top
                        + ( bins.empty[i] ? " NOFILL " : " FILL ") + " MESH"
                    );
                }
                Jmol.script(applet,
                    'COLOR $' + binName + " TRANSLUCENT .4 [ 169, 68, 66 ] [ 169, 68, 66 ] "
                );
            }
        }

        var drawBinMembers = function (applet, points, clear) {
            var memberPrefix = binPrefix + "members_";
            if(clear) {
                Jmol.script(applet, 'DELETE $' + memberPrefix + "*");
            } else {
                var script = []
                for(var i = 0; i < points.length; i++) {
                    script.push("draw " + memberPrefix + i
                        + " DIAMETER 20 POINTS [ { " + points[i].join(' ') + " } ] COLOR [ 245, 231, 158 ] "
                    );
                }
                Jmol.script(applet, script.join(';'));
            }
        };

        return function (applet, chainContext, index) {
            var options = {
                    oneRadialSize: chainContext.bin_data.num_bins_per_radial_unit,
                    numRadialBins: chainContext.bin_data.max_radial_bin,
                    numBins: chainContext.bin_data.num_bins,
                    templates: chainContext.templates,
                    frequencyColors: true,
                    includeEmpty: false
                },
                context = chainContext.contexts[index],
                bins = makeBins(context, options);
            drawBins(applet, bins, options);
            return {
                bins: bins,
                highlight: function (binIndex, clear) {
                    return highlightBin(applet, bins, binIndex, clear, options);
                },
                center: function () {
                    Jmol.script(applet, 'CENTER { ' + bins.center.join(' ') + ' }');
                    Jmol.script(applet, 'MOVETO 1 BACK');
                },
                drawMembers: function (points, clear) {
                    return drawBinMembers(applet, points, options);
                },
                setActiveBin: function (binIndex) {
                    if(binIndex === undefined || binIndex == null) {
                        return;
                     }
                    var members = [],
                        memberIndexes = context.inverse[binIndex] == undefined ? [] : context.inverse[binIndex];
                    for(var i = 0; i < memberIndexes.length; i++) {
                        members.push(chainContext.contexts[memberIndexes[i]].coords);
                    }
                    highlightBin(applet, bins, binIndex, true, options);
                    highlightBin(applet, bins, binIndex, false, options);
                    drawBinMembers(applet, [], true);
                    drawBinMembers(applet, members);
                },
                clearActiveBin: function (binIndex) {
                    highlightBin(applet, bins, binIndex, true, options);
                    drawBinMembers(applet, [], true);
                }
            };
        };

    }])

    .factory('ChainContextHeatMap', ['SparseHeatMap', function(SparseHeatMap) {
        return function(chainContext, index) {
            var options = {
                    oneRadialSize: chainContext.bin_data.num_bins_per_radial_unit,
                    numRadialBins: chainContext.bin_data.max_radial_bin,
                    numBins: chainContext.bin_data.num_bins,
                    zenithBins: chainContext.bin_data.bin_divisor,
                    azmuthBins: 2 * chainContext.bin_data.bin_divisor,
                    oneRowBins: chainContext.bin_data.num_bins / chainContext.bin_data.bin_divisor,
                    boundarySize: 1,
                    thumbnailSize: 4,
                    colorMap: [
                        "#FFF", "#DDD", "#BBB", "#999",
                        "#777", "#555", "#555", "#333",
                        "#333", "#333", "#111", "#111",
                        "#111", "#000", "#000", "#000"
                    ]
                };
            var thumbnails = [];
            for(var thumbIndex = 0; thumbIndex < chainContext.contexts.length; thumbIndex++) {
                thumbnails.push(SparseHeatMap(chainContext.contexts[thumbIndex], options));
            }

            var boundaries = [];
            for(var i = 1; i < options.numRadialBins; i++) {
                boundaries.push({
                    x: i * options.thumbnailSize * options.azmuthBins + (i-1) * options.boundarySize,
                    y: 0,
                    width: options.boundarySize,
                    height: options.zenithBins * options.thumbnailSize,
                    color: "#000"
                });
            }

            return {
                heatmaps: thumbnails,
                boundaries: boundaries,
                dimensions: {
                    width: options.thumbnailSize * options.oneRowBins + options.boundarySize * boundaries.length,
                    height: options.thumbnailSize * options.zenithBins
                }
            };
        };
    }])

    .factory('ContextSvgPaths', [function () {
        var scale = 25;
        return function (context) {
            var zenithBins = [],
                zenithBoundaries = [0].concat(context.bin_data.zenith_boundaries).concat(Math.PI),
                azimuthBoundaries = context.bin_data.azimuth_boundaries;
            for(var zI = 0; zI < zenithBoundaries.length; zI++) {
                var azimuthBins = [];
                for(var aI = 0; aI < azimuthBoundaries.length; aI++) {
                    var a = azimuthBoundaries[aI],
                        z = zenithBoundaries[zI],
                        r = scale * z + 10 + zI;
                    azimuthBins.push({
                        a: Math.round(a, 0),
                        z: Math.round(z, 0),
                        r: Math.round(r, 0),
                        x: Math.round(r * Math.cos(a), 0),
                        y: Math.round(r * Math.sin(a), 0)
                    });
                }
                zenithBins.push(azimuthBins);
            }

            var bins = [],
                index = 0;
            for(var i = 0; i < zenithBins.length-1; i++) {
                for(var j = 0; j < zenithBins[i].length; j++) {
                    var i1 = i+1,
                        j1 = (j+1) % azimuthBoundaries.length,
                        b00 = zenithBins[i][j],
                        b01 = zenithBins[i][j1],
                        b10 = zenithBins[i1][j],
                        b11 = zenithBins[i1][j1],
                        cx1 = (b00.x + b01.x) / 2,
                        cx2 = (b10.x + b11.x) / 2,
                        cy1 = (b00.y + b01.y) / 2,
                        cy2 = (b10.y + b11.y) / 2;
                    bins.push({
                        index: index,
                        azimuth: b00.a,
                        zenith: b00.z,
                        azimuthIndex: j,
                        zenithIndex: i,
                        grad: {
                            x1: cx1,
                            x2: cx2,
                            y1: cy1,
                            y2: cy2,
                        },
                        path: [
                            "M" + b00.x + "," + b00.y,
                            "L" + b10.x + "," + b10.y,
                            "A" + b10.r + "," + b11.r + " 0 0 1 " + b11.x + "," + b11.y,
                            "L" + b01.x + "," + b01.y,
                            "A" + b01.r + "," + b00.r + " 0 0 0 " + b00.x + "," + b00.y,
                            "Z"
                        ].join(' ')
                    });
                    index++;
                }
            }

            return bins;
        };
    }])
;
