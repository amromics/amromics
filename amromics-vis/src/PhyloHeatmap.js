/* eslint-disable */
import * as d3 from "d3";
import {
  NewickTools
} from "newick";
export class PhyloHeatmap {
  constructor(element) {
    this.container = element;
    this.props = {
      width: this.container.clientWidth,
      height: 400,
      
    };
    var treeview = document.createElement("div");
    treeview.id = "ph_treeview";
    treeview.style.width = this.props.width / 5 + "px";
    treeview.style.float="left";
    var heatmapview = document.createElement("div");
    heatmapview.id = "ph_heatmapview";
    heatmapview.style.width = (this.props.width-this.props.width / 5) + "px";
    heatmapview.style.float="left";
    heatmapview.style.overflowX = "scroll";
    heatmapview.style.position="relative";
    this.container.appendChild(treeview);
    this.container.appendChild(heatmapview);
    this.active_names=[];
    this.cell_size=20;
  }
  load(phylotree, hits) {
    this.hits = hits;    
    var margin = {
      top: 30,
      right: 30,
      bottom: 60,
      left: 30
    };
    //var width_tree =  this.props.width / 5 - margin.left - margin.right-50;
    var newick = NewickTools.parse(phylotree);
     //estimate length of sample name, by average length plus 10
     var numchar=0;
     var listsamples=Object.keys(NewickTools.dfs(newick));
    
     for (var s in listsamples){
       numchar=numchar+listsamples[s].length;
      
     }
     var namelength=numchar/listsamples.length*5+10;
     var width_tree = this.props.width / 5 - margin.left - margin.right-namelength;
     //  assigns the data to a hierarchy using parent-child relationships
    let nodes = d3.hierarchy(newick, d => d.branchset);
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
          n.data.name.replace(/\'/g, "")
        );
        num_leaf++;
      }
    }
    var height=num_leaf*this.cell_size;
    this.container.style.height= (height + 100) + "px";
    var distance_per_depth = width_tree / max_depth;
    stack = [];
    stack.push(nodes);
    count = 0;
    nodes.ay = 0;
    //travel again, alter position of nodes
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
    
    const treemap = d3.tree().size([height, width_tree]);
    this.fnodes = treemap(nodes);
    //collect item for x axist
    this.genes = [];
    for (var i = 0; i < this.hits.length; i++) {

        if (!this.genes.includes(this.hits[i].gene)) {
          this.genes.push(this.hits[i].gene);
        }
      }


  }
  setOptions(options) {
    this.props.width = options.width;
    this.props.height = options.height;

  }
  drawTree(){
    var active_names=this.active_names;
    document.getElementById("ph_treeview").innerHTML = "";
    var margin = {
      top: 30,
      right: 30,
      bottom: 60,
      left: 30
    };
    var width = this.props.width / 5 ;
    //var height = 400 - margin.top - margin.bottom;
    var height = this.node_leaf.length*this.cell_size;
    var height_rect =this.cell_size;
    var svg2 = d3
      .select("#ph_treeview")
      .append("svg")
      .attr("width", width )
      .attr("height", this.cell_size*this.node_leaf.length + margin.top + margin.bottom);
    var g = svg2
      .append("g")
      .attr(
        "transform",
        "translate(" + margin.left + "," + (margin.top + height_rect/2) + ")"
      );

    // adds the links between the nodes
    this.graphic_tree=g;
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
    var active_names=this.active_names;

    //console.log(active_names);
    const node = g
      .selectAll(".node")
      .data(this.fnodes.descendants())
      .enter()
      .append("g")
      .attr(
        "class",
        d => "node" + (d.children ? " node--internal" : " node--leaf") +(active_names.includes(d.data.name.replace(/\'/g,''))?" node--active":"")
      )
      .attr("transform", d => "translate(" + d.ay + "," + d.ax + ")");

    // adds the circle to the node

    const leftnode = g
      .selectAll(".node--leaf")
      .append("circle")
      .attr("r", 2)
      .style("stroke", "black")
      .style("fill", "black");
    const leftnode_active = g
      .selectAll(".node--active")
      .append("circle")
      .attr("r", 5)
      .style("stroke", "blue")
      .style("fill", "none");
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
      var y_data=this.arr_sample_from_tree.slice();
      y_data.reverse();
      // Build X scales and axis:
      var y = d3
        .scaleBand()
        .range([height, 0])
        .domain(y_data)
        .padding(0.01);
      svg2
      .append("g")
      .attr(
        "transform",
        "translate(" +
        (this.props.width / 5) +
        "," +
        margin.top +
        ")"
      )
      .call(d3.axisLeft(y));
  }

  drawHeatmap(){
    var margin = {
      top: 30,
      right: 30,
      bottom: 100,
      left:0
    };
    var height = this.node_leaf.length*this.cell_size;
    var width = this.genes.length*this.cell_size ;
    document.getElementById("ph_heatmapview").innerHTML = "";
    var svg = d3
      .select("#ph_heatmapview")
      .append("svg")
      .attr("width", width )
      .attr("height", height + margin.top + margin.bottom)
      .append("g")
      .attr("transform", "translate(" + (margin.left ) + "," + margin.top + ")");
    var x = d3
      .scaleBand()
      .range([0, this.genes.length*this.cell_size])
      .domain(this.genes)
      .padding(0.01);
    var selected_samples=this.active_names;  
    var highlight_genes=new Set();
    for (var i=0;i<this.hits.length;i++){
      //console.log(this.hits[i].sample);
      if (selected_samples.includes(this.hits[i].sample)){
        highlight_genes.add(this.hits[i].gene);
      }
    }
    //console.log(highlight_genes);
    var tick_x=svg
      .append("g")
      .attr("transform", "translate(" + 0 + "," + height + ") ")
      .call(d3.axisBottom(x))
      .selectAll("text")
      .style("text-anchor", "end")
      .attr("dx", "-.8em")
      .attr("dy", ".25em")
      .attr("transform", "rotate(-65)")
      .style("font-weight", function(d) {
        //console.log(d);
        if(highlight_genes.has(d)){
         
          return "bold";
        }
        else{
          //console.log("no bold"+d);
          return "normal";
        }
        
      })
     ;
    var selected_samples=this.active_names;
    var y_data=this.arr_sample_from_tree.slice();
    y_data.reverse();
    // Build X scales and axis:
    var y = d3
      .scaleBand()
      .range([height, 0])
      .domain(y_data)
      .padding(0.01);
    // svg.append("g").call(d3.axisLeft(y));
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
    svg
      .selectAll()
      .data(this.hits, function(d) {
       
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
            //return amrColor(d.class);
          if(selected_samples.includes(d.sample) )
            return "#6E1214";
          else
          return "#DD2429";
        }

      )
     .style("stroke", "white")
     .on("mouseover", mouseover)
     .on("mousemove", mousemove)
     .on("mouseleave", mouseleave);
  }
  setActiveNames(names){
    this.active_names=names;
    this.draw();
  }
  //Draw heatmap
  draw() {
    this.drawTree();
    this.drawHeatmap();
  }
}
export default PhyloHeatmap
