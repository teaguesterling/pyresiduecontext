'use strict';

angular.module('ResCtxVis.controllers', [
    'ResCtxVis.services'
])
    .controller('PdbSelector', ['$scope', '$location', function ($scope, $location) {
        $scope.viewContext = function(pdbid1, chain1, pdbid2, chain2) {
            if(!pdbid2 || !chain2) {
                $location.path("/view/" + pdbid1 + chain1);
            } else {
                $location.path("/compare/" + pdbid1 + chain1 + "-" + pdbid2 + chain2);
            }
        };
    }])

    .controller('View2', [
        '$scope',
        '$route',
        function($scope, $route) {
            $scope.identifier1 = $route.current.params.identifier1;
            $scope.identifier2 = $route.current.params.identifier2;
        }
    ])

    .controller('ViewAligned', [
        '$scope',
        '$route',
        function($scope, $route) {
            $scope.identifier1 = $route.current.params.identifier1;
            $scope.identifier2 = $route.current.params.identifier2;
            $scope.alignmentRunId = $route.current.params.alignmentRunId;
        }
    ])

    .controller('ContextViewer', [
        '$scope',
        '$route',
        '$interval',
        'ChainContextFactory',
        'GridPointContextFactory',
        'JsMolChainColoring',
        'ChainContextHeatMap',
        'ResidueContextHeatMap',
        'JsMolSphericalHistogram',
        'ContextSvgPaths',
        function (
            $scope,
            $route,
            $interval,
            ChainContextFactory,
            GridPointContextFactory,
            JsMolChainColoring,
            ChainContextHeatMap,
            ResidueContextHeatMap,
            JsMolSphericalHistogram,
            ContextSvgPaths) {
        if($scope.identifier === undefined) {
            $scope.identifier = $route.current.params.identifier;
        }
        if($scope.alignmentRunId === undefined) {
            $scope.alignmentRunId = $route.current.params.alignmentRunId || null;
        }

        var ChainContext = ChainContextFactory($scope.alignmentRunId),
            GridPointContext = GridPointContextFactory($scope.alignmentRunId);

        $scope.pdbid = $scope.identifier.slice(0, 4);
        $scope.chain = $scope.identifier[4];
        $scope.jsmol = null;
        $scope.active = {
            index: -1,
            item: null
        };
        $scope.focus = null;
        $scope.highlights = [];
        $scope.animationTimeout = null;
        $scope.sphericalBinData = null;
        $scope.alignmentData = null;

        var $document = angular.element(document);
        $document.on('keydown', $scope.handleKeyPress);
        $scope.$on('$destory', function () {
            $document.off('keydown',$scope.handleKeyPress);
        });

        $scope.toggleAnimate = function () {
            if($scope.animationTimeout === null) {
                $scope.animationTimeout = $interval(function () {
                    $scope.stepActiveContext(1);
                }, 100);
            } else {
                $interval.cancel($scope.animationTimeout);
                $scope.animationTimeout = null;
            }
        };

        $scope.stepActiveContext = function (step) {
            step = step || 1;
            var newIndex = $scope.active.index + step;
            if(newIndex < 0) {
                newIndex = $scope.context.contexts.length-1;
            } else if(newIndex >= $scope.context.contexts.length) {
                newIndex = 0;
            }
            $scope.setActiveContext(newIndex);
        };

        $scope.setActiveContext = function (index) {
            var currentBin = $scope.focus;
            $scope.active = {
                index: index,
                item: $scope.context.contexts[index],
                heatmap: ResidueContextHeatMap($scope.context, index)
            };

            var scale = 0;
            for(var hI = 0; hI < $scope.context.bin_data.num_bins; hI++) {
                scale += $scope.active.item.histogram[hI] === undefined ? 0 : $scope.active.item.histogram[hI];
            }
            scale = 20 * 255 / scale;
            $scope.getContextBinColor = function (val) {
                return Math.max(0, val !== undefined ? Math.round(255 - val * scale,0) : 255);
            };

            $scope.gridHistograms = null;
            GridPointContext.get(
                {
                    identifier: $scope.identifier,
                    coords: $scope.active.item.coords.join(' ')
                }, function (data) {
                    $scope.gridHistograms = data.histograms;
                }
            );
            $scope.resetHighlights();
            $scope.sphericalBinData = JsMolSphericalHistogram($scope.jsmol, $scope.context, index);
            $scope.sphericalBinData.center();

            if(currentBin !== null && currentBin !== undefined) {
                $scope.setActiveBinHighlights(currentBin);
            }
        };

        $scope.changeContext = function(pdbid, chain) {
            $location.path("/view/" + pdbid + chain);
        };

        $scope.setActiveBinHighlights = function (binNumber) {
            if($scope.focus !== undefined) {
                var current = $scope.getContextBinColor($scope.active.item.histogram[$scope.focus]);
                angular.element("#polar-bin-" + $scope.pdbid + "-" + $scope.focus).css({
                    'fill': 'rgb(' + [current, current, current].join(',') + ')',
                    'stroke': '#ddd'
                });
            }
            $scope.focus = binNumber;
            $scope.highlights = $scope.active.item.inverse[binNumber] === undefined ? [] :
                $scope.active.item.inverse[binNumber];
            if($scope.sphericalBinData !== null) {
                $scope.sphericalBinData.setActiveBin(binNumber);
            }

            angular.element("#polar-bin-" + $scope.pdbid + "-" + binNumber).css({
                'fill': '#a94442',
                'stroke': '#000'
            });
            $scope.$apply();
        };

        $scope.resetHighlights = function (force) {
            if($scope.focus !== undefined && $scope.focus !== null) {
                var current = $scope.getContextBinColor($scope.active.item.histogram[$scope.focus]);
                angular.element("#polar-bin-" + $scope.pdbid + "-" + $scope.focus).css({
                    'fill': 'rgb(' + [current, current, current].join(',') + ')',
                    'stroke': '#ddd'
                });
            }
            if($scope.sphericalBinData !== null) {
                $scope.sphericalBinData.clearActiveBin($scope.focus);
            }
            $scope.highlights = [];
            $scope.focus = null;
        };

        $scope.toggleFocus = function (binIndex) {
            if(binIndex === null) {
                $scope.resetHighlights();
            } else {
                $scope.setActiveBinHighlights(binIndex);
            }
        };

        $scope.clientTiming = {};
        var start = Date.now();

        $scope.status = "Downloading Alignment Data";
        JsMolChainColoring($scope.jsmol, $scope);

        $scope.status = "Downloading PDB/Generating Residue Contexts";
        start = Date.now();
        ChainContext.get({identifier: $scope.identifier}, function (context) {
            $scope.clientTiming['download'] = (Date.now() - start) / 1000;
            $scope.clientTiming['transmission'] = ($scope.clientTiming['download'] -
                context.timing.download - context.timing.generate - context.timing.render
            );
            $scope.status = "Generating Residue Context Thumbprints";
            $scope.context = context;
            start = Date.now();
            $scope.thumbnails = ChainContextHeatMap($scope.context);
            $scope.polarBinsTemplate = ContextSvgPaths(context);
            $scope.clientTiming['drawing'] = (Date.now() - start) / 1000;
            $scope.status = "";


        });
    }])

    .controller('Drawer', ['$scope', 'SparseHeatMap', 'DenseHeatMap',
    function ($scope, SparseHeatMap, DenseHeatMap) {

    }])
;
