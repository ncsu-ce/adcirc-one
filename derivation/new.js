var tags = {};
var equations = {};
var num_equations = 0;

// Render equations
d3.selectAll('.equation')
    .each(render_block);

d3.selectAll('p')
    .each(render_inline);
d3.selectAll('li')
    .each(render_inline);

// Generate tags
d3.selectAll('.equation')
    .each(generate_tag);

d3.selectAll('a.tag')
    .each(generate_tag);

// Reference tags
d3.selectAll('a')
    .each(reference_tag);


function reference_tag () {

    var a = d3.select(this);
    var tag = a.attr('data-tag');

    if (tag) {

        a.attr('href', '#' + tag)
            .style('white-space', 'nowrap');

        // Reference generic tag
        if (tag in tags) {
            a.text(tags[tag]);
        }

        // Reference tag for equation
        if (tag in equations) {
            a.text('(' + equations[tag] + ')');
        }

    }

}

function generate_tag () {

    var element = d3.select(this);

    // Generate tag for element with the .tag class
    // These will generally be # characters next to
    // section headers.
    if (element.classed('tag')) {
        var name = element.attr('name');
        if (name) {
            element.attr('href', '#' + name);
            tags[name] = name;
        }

    }

    // Generate tags for equations
    if (element.classed('equation')) {

        var tag = element.attr('data-tag');
        var num = num_equations + 1;
        if (tag) {

            element.insert('a', ':first-child')
                .classed('eq-tag', true)
                .attr('name', tag)
                .attr('href', '#' + tag)
                .text('(' + num + ')');

            equations[tag] = num;
            num_equations += 1;

        }

    }

}

function render_block () {

    var div = d3.select(this);
    var text = div.text().trim();
    var block = '\\begin{aligned}' + text + '\\end{aligned}';
    katex.render(block, div.node(), { displayMode: true });

}

function render_inline () {

    var selection = d3.select(this);
    var text = selection.html();
    var regex = /\\\((.*?)\\\)/g;

    text = text.replace(regex, function (something) {
        return katex.renderToString(something.slice(2, -2));
    });

    selection.html(text);

}

d3.selectAll('.inset')
    .each(function () {

        // Select the element
        var div = d3.select(this);

        // Create the new parent div
        var parent = d3.select(document.createElement('div'))
            .attr('class', 'expander');

        // Add the expander button
        var button = parent.append('span')
            .attr('class', 'fa-layers fa-fw fa-lg expander-button');

        // Add the icons to the expander button
        button.append('i')
            .attr('class', 'fas fa-circle background')
            .attr('data-fa-transform', 'grow-18')
            .style('filter', 'drop-shadow(0 0.25em 2px #dedede)');
        button.append('i')
            .attr('class', 'fas fa-minus-circle');

        // Add the container div immediately before this one
        this.parentNode.insertBefore(parent.node(), this);

        // Remove the div and place into the parent
        parent.append(function () {
            return div.remove().node();
        });

        // Hook up the button
        var collapsed = false;
        button.on('click', function () {

            var fold_time = 500;
            var height = null;

            if (!collapsed) {

                // Get the computed height of the div
                height = div.style('height');
                div.datum(height);

                // Collapse the div
                div.style('height', height)
                    .transition(fold_time)
                    .style('height', '0px');

                // Set the icon
                button.selectAll('.fa-minus-circle')
                    .classed('fa-minus-circle', false)
                    .classed('fa-plus-circle', true);

                // Remove the drop shadow
                button.selectAll('.background')
                    .style('filter', null);
                console.log(button);

            } else {

                // Get the height from before the div was collapsed
                height = div.datum();

                // Expand the div
                div.style('height', '0px')
                    .transition(fold_time)
                    .style('height', height);

                // Set the icon
                button.selectAll('.fa-plus-circle')
                    .classed('fa-plus-circle', false)
                    .classed('fa-minus-circle', true);

                button.selectAll('.background')
                    .style('filter', 'drop-shadow(0 0.25em 2px #dedede)')

            }

            collapsed = !collapsed;

        });


    });

hljs.initHighlightingOnLoad();