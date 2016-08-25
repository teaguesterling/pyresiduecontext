'use strict';

angular.module('ResCtxVis', [
        'ngRoute',
        'ngResource',
        'ResCtxVis.services',
        'ResCtxVis.controllers',
        'ResCtxVis.filters',
        'ResCtxVis.directives'
    ])

    .config(['$routeProvider', '$locationProvider', function ($routeProvider, $locationProvider) {
        $routeProvider
            .when('/', {
                templateUrl: '/static/tpl/input.html',
                controller: 'PdbSelector'
            })
            .when('/aligned/:alignmentRunId/:identifier1-:identifier2', {
                templateUrl: '/static/tpl/view-aligned.html',
                controller: 'ViewAligned'
            })
            .when('/compare/:identifier1-:identifier2', {
                templateUrl: '/static/tpl/view-2.html',
                controller: 'View2'
            })
            .when('/view/:identifier', {
                templateUrl: '/static/tpl/view-1.html',
                controller: 'ContextViewer'
            })
            .when('/draw/', {
                templateUrl: '/static/tpl/drawer.html',
                controller: 'Drawer'
            })
            .otherwise({
                redirectTo: '/'
            });

        $locationProvider.html5Mode(false);
    }])

    .run(['$rootScope', '$window', '$templateCache', function ($rootScope, $window, $templateCache) {
        $rootScope.$on('$viewContentLoaded', function() {
            $templateCache.removeAll();
        });
    }])
;
