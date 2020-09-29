/* eslint-disable */
import * as d3 from "d3";
import {
  NewickTools
} from "newick";
export class PhyloHeatmap {
  constructor(element) {
    this.container = element;
    this.props = {
      width: 800,
      height: 400
    };
    var treeview = document.createElement("div");
    treeview.id = "ph_treeview";
    treeview.style.width = this.props.width / 5 + "px";
    treeview.style.height = this.props.height + "px";
    treeview.style.float="left";
    var heatmapview = document.createElement("div");
    heatmapview.id = "ph_heatmapview";
    heatmapview.style.width = (this.props.width-this.props.width / 5) + "px";
    heatmapview.style.height = this.props.height + "px";
    heatmapview.style.float="left";
    heatmapview.style.position="relative";
    this.container.appendChild(treeview);
    this.container.appendChild(heatmapview);


  }
  load(phylotree, heatmap) {
    this.heatmap = heatmap;


    var margin = {
      top: 30,
      right: 30,
      bottom: 60,
      left: 30
    };
    var width = 200 - margin.left - margin.right;
    var height = 400 - margin.top - margin.bottom;


    var newick = NewickTools.parse(phylotree);

    //console.log(newick);
    // declares a tree layout and assigns the size
    const treemap = d3.tree().size([height, width]);

    //  assigns the data to a hierarchy using parent-child relationships
    let nodes = d3.hierarchy(newick, d => d.branchset);
    //console.log(nodes);
    // maps the node data to the tree layout
    let dist = 20;

    var stack = [];
    nodes.height = 0;
    stack.push(nodes);
    var max_heigth = 0;
    var max_depth = 0;
    //travel to cal height of nodes
    let count = 0;
    var num_leaf = 0;
    this.node_leaf = [];
    this.arr_sample_from_tree = [];
    while (stack.length > 0) {
      var n = stack.pop();
      n.id = count;

      count++;
      if (n.children != undefined) {
        for (var i = 0; i < n.children.length; i++) {
          stack.push(n.children[i]);
          n.children[i].height = n.height + n.children[i].data.length;
        }
      } else {
        if (n.height > max_heigth) max_heigth = n.height;
        if (n.depth > max_depth) max_depth = n.depth;
        this.node_leaf.push(n);
        this.arr_sample_from_tree.push(
          n.data.name.replace("'", "").replace("'", "")
        );
        num_leaf++;
      }
    }

    var distance_per_depth = width / max_depth;
    var unit_distance = width / max_heigth;
    stack = [];
    stack.push(nodes);
    count = 0;
    nodes.ay = 0;
    while (stack.length > 0) {
      var n = stack.pop();

      if (n.children != undefined) {
        for (var i = 0; i < n.children.length; i++) {
          stack.push(n.children[i]);
        }
        //setup node by depth
        n.ay = n.depth * distance_per_depth;
      } else {
        n.ax = (count * height) / num_leaf;
        count++;
        //setup node by depth
        n.ay = max_depth * distance_per_depth;
      }
      //setup node by heigth
      //n.ay = n.height * unit_distance;
    }
    for (var d = max_depth; d >= 0; d--) {
      stack = [];
      stack.push(nodes);
      while (stack.length > 0) {
        var n = stack.pop();

        if (n.children != undefined) {
          for (var i = 0; i < n.children.length; i++) {
            stack.push(n.children[i]);
            if (n.depth == d) {
              n.ax = (n.children[0].ax + n.children[1].ax) / 2;
            }
          }
        }
      }
    }
    //console.log(nodes);
    //travel again, alter position of nodes

    //
    this.fnodes = treemap(nodes);

    //split amr and vir:

    this.vir_genes = [];
    this.amr_genes = [];
    for (var i = 0; i < heatmap.hits.length; i++) {
      if (heatmap.hits[i].type == "amr") {

        //add vir gene to list gene (need remove when result got list vir genes)
        if (!this.amr_genes.includes(this.heatmap.hits[i].gene)) {
          this.amr_genes.push(this.heatmap.hits[i].gene);
        }
      }
      if (heatmap.hits[i].type == "vir") {

        //add vir gene to list gene (need remove when result got list vir genes)
        if (!this.vir_genes.includes(this.heatmap.hits[i].gene)) {
          this.vir_genes.push(this.heatmap.hits[i].gene);
        }
      }

    }
    this.genes = [];
    this.genes = this.genes.concat(this.vir_genes).concat(this.amr_genes);
    this.data = heatmap.hits;
  }
  setOptions(options) {
    this.props.width = options.width;
    this.props.height = options.height;

  }
  draw() {
    var margin = {
      top: 30,
      right: 30,
      bottom: 60,
      left: 30
    };
    var width = 200 - margin.left - margin.right;
    var height = 400 - margin.top - margin.bottom;
    var height_rect = height / this.node_leaf.length / 2;

    var svg2 = d3
      .select("#ph_treeview")
      .append("svg")
      .attr("width", width + margin.left + margin.right)
      .attr("height", height + margin.top + margin.bottom);
    var g = svg2
      .append("g")
      .attr(
        "transform",
        "translate(" + margin.left + "," + (margin.top + height_rect) + ")"
      );

    // adds the links between the nodes

    const link = g
      .selectAll(".link")
      .data(this.fnodes.descendants().slice(1))
      .enter()
      .append("path")
      .attr("class", "link")
      .style("stroke", "black")
    .style("fill", "none")

      .attr("d", d => {
        return (
          "M" +
          d.ay +
          "," +
          d.ax +
          "L" +
          d.parent.ay +
          "," +
          d.ax +
          " " +
          d.parent.ay +
          "," +
          d.parent.ax
        );
      });

    // adds each node as a group
    const node = g
      .selectAll(".node")
      .data(this.fnodes.descendants())
      .enter()
      .append("g")
      .attr(
        "class",
        d => "node" + (d.children ? " node--internal" : " node--leaf")
      )
      .attr("transform", d => "translate(" + d.ay + "," + d.ax + ")");

    // adds the circle to the node
    const leftnode = g
      .selectAll(".node--leaf")
      .append("circle")
      .attr("r", 2)
      .style("stroke", "black")
      .style("fill", "black");
    //add connectors  in case draw phylogeny by length
    const connectors = g
      .selectAll(".connector")
      .data(this.node_leaf)
      .enter()
      .append("path")
      .attr("class", "connector")

      .attr("d", d => {
        return (
          "M" + d.ay + "," + d.ax + "L" + (width + margin.right) + "," + d.ax
        );
      });
    // append the svg object to the body of the page
    width = 800 - margin.left - margin.right;
    var svg = d3
      .select("#ph_heatmapview")
      .append("svg")
      .attr("width", width + margin.left + margin.right)
      .attr("height", height + margin.top + margin.bottom)
      .append("g")
      .attr("transform", "translate(" + (margin.left + 50) + "," + margin.top + ")");

    var x = d3
      .scaleBand()
      .range([0, width])
      .domain(this.genes)
      .padding(0.01);
    svg
      .append("g")
      .attr("transform", "translate(" + 0 + "," + height + ") ")
      .call(d3.axisBottom(x))
      .selectAll("text")
      .style("text-anchor", "end")
      .attr("dx", "-.8em")
      .attr("dy", ".25em")
      .attr("transform", "rotate(-65)");;

    // Build X scales and axis:
    var y = d3
      .scaleBand()
      .range([height, 0])
      .domain(this.arr_sample_from_tree)
      .padding(0.01);
    svg.append("g").call(d3.axisLeft(y));
    var amrColor = d3
      .scaleOrdinal()
      .domain(this.heatmap.list_class)
      .range([
        "gold",
        "peru",
        "chocolate",
        "yellow",
        "red",
        "pink",
        "#C83200",
        "#CD3700",
        "#FF6103",
        "#CC7722",
        "#FD6302",
        "#883000",
        "#FFBF00",
        "#CF9812A",
        "coral",
        "pumpkin",
        "tomato",
        "brown",
        "vermilion",
        "orange red",
        "orange",
        "crimson",
        "dark red",
        "hot pink",
        "smitten",
        "magenta",
        "indigo",
        "blue violet"
      ]);

    // create a tooltip
    var tooltip = d3
      .select("#ph_heatmapview")
      .append("div")
      .style("opacity", 0)
      .attr("class", "tooltip")
      .style("position","absolute")
      .style("background-color", "white")
      .style("border", "solid")
      .style("border-width", "2px")
      .style("border-radius", "5px")
      .style("padding", "5px");

    // Three function that change the tooltip when user hover / move / leave a cell
    var mouseover = function(d) {
      tooltip.style("opacity", 1);

      d3.select(this)
        .style("stroke", "black")
        .style("opacity", 1);

    };
    var mousemove = function(d) {
      tooltip
        .html(
          d.gene +
          " in " +
          d.sample +
          "<br>" +
          "Identity: " +
          d.identity +
          "<br>" + (d.type == "amr" ? d.class : d.product)

        )
        .style("left", d3.mouse(this)[0] + 70 + "px")
        .style("top", d3.mouse(this)[1] + "px");
    };
    var mouseleave = function(d) {
      tooltip.style("opacity", 0);
      d3.select(this).style("stroke", "white");
    };

    // add the squares
    console.log(this.amr_genes);
    var _amr_genes=this.amr_genes;
    var _data=this.data;
    svg
      .selectAll()
      .data(this.data, function(d) {
        return d.sample + ":" + d.gene;
      })
      .enter()
      .append("rect")
      .attr("x", function(d) {
        return x(d.gene);
      })
      .attr("y", function(d) {
        return y(d.sample);
      })
      .attr("width", x.bandwidth())
      .attr("height", y.bandwidth())

      .style("fill", function(d) {
        if (d.type == "amr" || _amr_genes.includes(d.gene)) {
          if (d.type == "amr")
            return amrColor(d.class);
          else {
            //find class for amr gene
            for (var i = 0; i < _data.length; i++) {
              if (_data[i].gene == d.gene && _data[i].type == "amr") {
                return amrColor(_data[i].class);
              }
            }
          }
        } else {

          return "lime";
        }

      })
      .style("stroke", "white")
      .on("mouseover", mouseover)
      .on("mousemove", mousemove)
      .on("mouseleave", mouseleave);
  }
}
export default PhyloHeatmap
