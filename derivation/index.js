var equations = {
    'e1': '\\ddt{H} + \\ddx{}(UH) = 0',
    'e2':
        '\\begin{aligned} ' +
        'H		&\\equiv	\\zeta + h \\\\ ' +
        '\\zeta	&=		\\text{free surface departure from the geoid} \\\\ ' +
        'h 		&=		\\text{bathymetric depth (distance from the geoid to the bottom)} \\\\' +
        'u 		&=		\\text{vertically varying velocity in the x-direction} \\\\' +
        'U 		&=		\\frac{1}{H}\\int_{-h}^\\zeta u \\mathrm{d}z = \\text{depth-averaged velocity in the x-direction} \\\\' +
        '\\end{aligned}'
};

var macros = {
    '\\ddt': '\\frac{\\partial #1}{\\partial t}',
    '\\dddt': '\\frac{\\partial^2 #1}{\\partial t^2}',
    '\\ddx': '\\frac{\\partial #1}{\\partial x}',
    '\\dddx': '\\frac{\\partial^2 #1}{\\partial x^2}',
    '\\ddxt': '\\frac{\\partial #1}{\\partial x\\partial t}',
    '\\Jx': '\\tilde{J}_x'
};

function render_equations ( equations ) {

    for ( var e in equations ) {

        if ( equations.hasOwnProperty( e ) ) {

            d3.select('#' + e).call(render_equation, equations[e]);

        }

    }

}

function render_equation (selection, equation) {

    katex.render(equation, selection.node(), { displayMode: true, macros: macros } );

}

render_equations(equations);