<style>
.el-row {
  margin-bottom: 20px;
  &:last-child {
    margin-bottom: 0;
  }
}

.el-col {
  border-radius: 4px;
}

.bg-purple-dark {
  background: #99a9bf;
}

.bg-purple {
  background: #d3dce6;
}

.bg-purple-light {
  background: #e5e9f2;
}

.grid-content {
  border-radius: 4px;
  min-height: 36px;
}

.row-bg {
  padding: 10px 0;
  background-color: #f9fafc;
}

.text {
  font-size: 14px;
}

.item {
  margin-bottom: 18px;
}

.clearfix:before,
.clearfix:after {
  display: table;
  content: "";
}

.clearfix:after {
  clear: both;
}

.box-card {
  width: 100%;
}

.axis path,
.axis line {
  fill: none;
  stroke: black;
  shape-rendering: crispEdges;
}

.axis text {
  font-family: sans-serif;
  font-size: 11px;
}

.node circle {
  fill: #fff;
  stroke: steelblue;
  stroke-width: 3px;
}

.node text {
  font: 16px sans-serif;
}

.link {
  fill: none;
  stroke: #000;
  stroke-width: 2px;
}
.connector {
  fill: none;
  stroke: #999;
  stroke-width: 1px;
}
.container {
  width: 1200px;
}
.sub-container-1 {
  width: 200px;
  float: left;
}
.sub-container-2 {
  width: 800px;
  float: left;
}
</style>

<template>
  <div v-loading="loading">
  
  </div>
</template>

<script>

export default {
  name: "PhyloHeatmap",

  data() {
    return {
      loading: false
    };
  },

  computed: {
    collectionId() {
      //console.log("call collection id");
      return this.$route.params.c;
    }
  },
  async created() {
    this.loading = true;
  },
  mounted() {
    this.loadPhyloHeatmap();
    this.loading = false;
  },
  methods: {
    loadPhyloHeatmap() {
      // set the dimensions and margins of the graph

      //Read the data
      //phylogeny_tree
      let list_s = this.list_samples;
      console.log("list sample");
      console.log(list_s);

      //get phylodata:
     // var ret = this.phyloheatmapData;
      var ret_str=JSON.stringify(this.phyloheatmapData);
       //console.log("before replace name");
      // console.log(ret_str);
      for (var i = 0; i < list_s.length; i++) {
        console.log(list_s[i].sample_id+"->"+list_s[i].sample_name);
        ret_str=ret_str.split(list_s[i].sample_id).join(list_s[i].sample_name);
      }
      var ret=JSON.parse(ret_str);
      //console.log("ret");
      //console.log(ret);
      // Labels of row and columns
      var genes = ret.list_genes;
      var samples = ret.samples;
      var margin = {
          top: 30,
          right: 30,
          bottom: 60,
          left: 30
        },
        width = 200 - margin.left - margin.right,
        height = 400 - margin.top - margin.bottom;
      //phylotree
      //phylogeny_tree
      var newick_raw = this.phylotree;
      //(newick_raw);
      for (var i = 0; i < list_s.length; i++) {
        newick_raw = newick_raw.replace(
          list_s[i].sample_id + "_contigs.fasta",
          list_s[i].sample_name
        );

        //console.log(list_s[i].sample_id+":"+list_s[i].sample_name);
      }
      newick_raw = newick_raw.replace(".ref", "");
      //console.log("after rename");
      //console.log(newick_raw);
      var newick = NewickTools.parse(newick_raw);

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
      var node_leaf = [];
      var arr_sample_from_tree = [];
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
          node_leaf.push(n);
          arr_sample_from_tree.push(
            n.data.name.replace("'", "").replace("'", "")
          );
          num_leaf++;
        }
      }
      //console.log("arr from tree");
      //console.log(arr_sample_from_tree);
      // edit leaft position y (x coordinate):
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
      var fnodes = treemap(nodes);
      ///console.log(fnodes);
      // append the svg object to the body of the page
      // appends a 'group' element to 'svg'
      // moves the 'group' element to the top left margin
      var height_rect = height / node_leaf.length / 2;
      var svg2 = d3
        .select("#my_tree")
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
        .data(fnodes.descendants().slice(1))
        .enter()
        .append("path")
        .attr("class", "link")

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
        .data(fnodes.descendants())
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
        .data(node_leaf)
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
        .select("#my_dataviz")
        .append("svg")
        .attr("width", width + margin.left + margin.right)
        .attr("height", height + margin.top + margin.bottom)
        .append("g")
        .attr("transform", "translate(" + (margin.left+50) + "," + margin.top + ")");
      var data = ret.hits;
      //split amr and vir:
    
      var vir_genes=[];
      var amr_genes=[];
      for (var i =0;i<ret.hits.length;i++){
        if (ret.hits[i].type=="amr"){
      
          //add vir gene to list gene (need remove when result got list vir genes)
          if (!amr_genes.includes(ret.hits[i].gene)){
            amr_genes.push(ret.hits[i].gene);
          }
        }
        if (ret.hits[i].type=="vir"){
      
          //add vir gene to list gene (need remove when result got list vir genes)
          if (!vir_genes.includes(ret.hits[i].gene)){
            vir_genes.push(ret.hits[i].gene);
          }
        }
          
      }
      genes=genes.concat(vir_genes);
      var x = d3
        .scaleBand()
        .range([0, width])
        .domain(genes)
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
        .domain(arr_sample_from_tree)
        .padding(0.01);
      svg.append("g").call(d3.axisLeft(y));
      var amrColor = d3
        .scaleOrdinal()
        .domain(ret.list_class)
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
        .select("#my_dataviz")
        .append("div")
        .style("opacity", 0)
        .attr("class", "tooltip")
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
              "<br>" + (d.type=="amr"?d.class:d.product)
             
          )
          .style("left", d3.mouse(this)[0] + 70 + "px")
          .style("top", d3.mouse(this)[1] + "px");
      };
      var mouseleave = function(d) {
        tooltip.style("opacity", 0);
        d3.select(this).style("stroke", "white");
      };

      // add the squares
      svg
        .selectAll()
        .data(data, function(d) {
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
          if (d.type=="amr" ||amr_genes.includes(d.gene)){
            if(d.type=="amr")
              return amrColor(d.class);
            else{
              //find class for amr gene
              for (var i =0;i<data.length;i++){
                if(data[i].gene==d.gene && data[i].type=="amr" ){
                  return amrColor(data[i].class);
                }
              }
            }
          }
            
          else{
            
            return "lime";
          }
            
        })
        .style("stroke", "white")
        .on("mouseover", mouseover)
        .on("mousemove", mousemove)
        .on("mouseleave", mouseleave);

      // adds the text to the node
      // node
      //   .append("text")
      //   .attr("dy", ".35em")
      //   .attr("x", d =>
      //     d.children ? (d.data.length + 5) * -1 : d.data.length + 5
      //   )
      //   .attr("y", d =>
      //     d.children && d.depth !== 0 ? -(d.data.length + 5) : d
      //   )
      //   .style("text-anchor", d => (d.children ? "end" : "start"))
      //   .text(d => d.data.name);
    }
  }
};
</script>
