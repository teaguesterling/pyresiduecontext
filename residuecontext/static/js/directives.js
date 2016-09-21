'use strict';
(function() {
    var directives = angular.module('ResCtxVis.directives', [])

    .filter('round', [function () {
        return function (input,places) {
            if(places === undefined) {
                return Math.round(input);
            } else {
                var scale = Math.pow(10, places)
                return Math.round(scale * input) / scale;
            }
        };
    }])

    .filter('min', [function () {
        return function (input,min) {
            return Math.min(input, min);
        };
    }])

    .filter('joinBy', [function () {
        return function (input,delimiter) {
            return (input || []).join(delimiter || ',');
        };
    }])

    .directive('ngJsMol', ['$window', function ($window) {
        var JSMOL_TPL_NUMBER = 1;
        var JSMOL_LINK_NUMBER = 1;

        return {
            restrict: 'A',
            scope: false,
            template: function (elem, attr) {
                return '<div id="' + attr.id + (JSMOL_TPL_NUMBER++) + '-jsmol">Loading...</div>';
            },
            link: function(scope, elem, attrs) {
                var structurePath;
                if(scope.alignmentRunId) {
                    structurePath = "/alignments/" + scope.alignmentRunId + "/pdb/" + attrs.structure + ".pdb";
                } else {
                    structurePath = "/pdb/" + attrs.structure + ".pdb";
                }

                var script = [
                        "load " + structurePath,
                        //"delete not backbone or not :" + attrs.chain,
                        "display protein",
                        "display :" + attrs.chain,
                        "backbone only",
                        "backbone 40",
                        //"sticks",
                        "set antialiasDisplay ON",
                        "set antialiasTranslucent ON",
                        "moveto 1 BACK"
                    ].join(';'),
                    options = {
                        width: "100%",
                        height: 500,
                        debug: false,
                        use: "JAVA HTML5 WEBGL",   // JAVA HTML5 WEBGL are all options
                        jarPath: "/static/lib/jsmol/java",
                        j2sPath: "/static/lib/jsmol/j2s", // this needs to point to where the j2s directory is.
                        jarFile: "JmolAppletSigned0.jar",
                        disableJ2SLoadMonitor: true,
                        disableInitialConsole: true,
                        addSelectionOptions: false,
                        readyFunction: function (applet) {
                            scope.jmol = applet;
                        },
                        isSigned: true,
                        disableJ2SLoadMonitor: false,
                        disableInitialConsole: false,
                        allowJavaScript: true,
                        coverTitle: "JSMOl initializing",
                        zIndexBase: 10,
                        script: script
                    };
                var id = attrs.id + (JSMOL_LINK_NUMBER++);
                scope.jsmol = Jmol.getApplet(id, options);
                scope.__proto__.jsmol = scope.jsmol;  // HACK!
                elem.find('#' + id + '-jsmol').html(Jmol.getAppletHtml(scope.jsmol));

                if(scope.setAlignmentColors && scope.alignmentData) {
                   scope.setAlignmentColors(scope.jsmol, scope.alignmentData);
                }
            }
        };
    }])

   .directive('ngSvg', function () {
      return {
        restrict: 'E',
        scope: true,
        transclude: true,
        replace: true,
        template: '<div ng-transclude></div>',
        controller: 'ngSvgController'
      };
    })

   .controller('ngSvgController', ['$scope', '$element', function ($scope, $element) {
       var self = this;
       // yes, belongs in link, but easier here
       $element.svg({onLoad: function (svg) {
         self.svg = svg;
         self.paths = svg.defs('paths');
       }});
   }])

    .directive('ngSvg', ['$compile', function($compile) {
        return {
            restrict: 'A',
            link: {
                post: function (scope, element, attrs) {

                }
            }
        }
    }])

    .directive('ngPolarBin', ['$compile', function($compile) {
        return {
            link: function (scope, element, attrs) {
                var index = parseInt(attrs['ngPolarBin']),
                    item = scope.active.item,
                    value = attrs['ngPolarBinValue'];
                element.click(function(e) {
                    scope.toggleFocus(index);
                    e.stopPropagation();
                });
                value = value === undefined ? 0 : value;
                if(typeof value === 'string' && value && value[0] == '#') {
                    element.attr('fill', 'url(' + scope.absUrl + value + ')');
                } else if(value > 0) {
                    var color = color = scope.getContextBinColor(value);
                    element.attr('fill', 'rgb(' + ([color, color, color].join(',')) + ')');
                }
            }
        };
    }])

    .directive('ngPolarBinGrad', ['$compile', function($compile) {
        return {
            link: function (scope, element, attrs) {
                var index = parseInt(attrs['ngPolarBin']),
                    item = scope.active.item,
                    value = item.histogram[index] === undefined ? 0 : item.histogram[index],
                    color = scope.getContextBinColor(value);
                element.click(function(e) {
                    scope.toggleFocus(index);
                    e.stopPropagation();
                });
                if(value > 0) {
                    element.css({
                        fill: 'rgb(' + ([color, color, color].join(',')) + ')'
                    });
                }
                //attrs.$set('title', 'boo!');
            }
        };
    }])

    .directive('ngHistogramBin', ['$compile', function($compile) {
        return {
            link: function (scope, element, attrs) {
                var index = parseInt(attrs['ngPolarBin']),
                    item = scope.active.item,
                    value = item.histogram[index] === undefined ? 0 : item.histogram[index],
                    color = scope.getContextBinColor(value);
                element.click(function(e) {
                    scope.toggleFocus(index);
                    e.stopPropagation();
                });
                if(value > 0) {
                    element.css({
                        fill: 'rgb(' + ([color, color, color].join(',')) + ')'
                    });
                }
                //attrs.$set('title', 'boo!');
            }
        };
    }])

    // Add special attributes to handle the mess of SVG attribute errors
    angular.forEach(['x', 'y', 'width', 'height', 'd', 'offset', 'color', 'stop-opacity'], function (attrName) {
        var directiveName = 'ng' + attrName[0].toUpperCase() + attrName.slice(1);
        directives.directive(directiveName, function () {
            return function (scope, element, attrs) {
                attrs.$observe(directiveName, function (value) {
                    attrs.$set(attrName, value);
                });
            };
        });
    });


})();

