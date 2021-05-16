/* eslint-disable */
import * as d3 from "d3";
import {
  NewickTools
} from "newick";
export class GeneFlowTree {
  nodes=[];
  links=[];
  constructor(element) {
    this.container = element;
    this.props = {
      width: this.container.clientWidth,
      height: 1200
    };
    this.treeview = document.createElement("div");
    this.treeview.id = "amr_gf_treeview";
    this.treeview.style.width = this.props.width + "px";
    this.treeview.style.height = this.props.height + "px";


    var control_div = document.createElement('div');
    control_div.style.width = (this.props.width) + "px";
    control_div.style.height = "40px";
    control_div.style.padding = "0px";
    control_div.style.margin = "0px";
    control_div.style.borderBottom = "1px solid #a6a6a6";

    var control_zoom = document.createElement('div');
    control_zoom.style.margin = "10px";
    control_zoom.style.float = "right";
    var control_zoom_in = document.createElement('button');
    control_zoom_in.style.backgroundImage = "none";
    control_zoom_in.style.border = "1px solid transparent";
    control_zoom_in.style.borderRadius = "4px";
    control_zoom_in.textContent = "+";
    control_zoom_in.addEventListener("click", this.zoomIn.bind(this));
    //control_zoom_in.addEventListener("click", this.zoomIn.bind(this));
    var control_zoom_out = document.createElement('button');
    control_zoom_out.style.backgroundImage = "none";
    control_zoom_out.style.border = "1px solid transparent";
    control_zoom_out.style.borderRadius = "4px";
    control_zoom_out.textContent = "-";
    control_zoom_out.addEventListener("click", this.zoomOut.bind(this));
    //  control_zoom_out.addEventListener("click", this.zoomOut.bind(this));
    control_zoom.appendChild(control_zoom_out);
    control_zoom.appendChild(control_zoom_in);
    control_div.appendChild(control_zoom);
    this.control_gene = document.createElement('div');
    this.control_gene.style.margin = "10px";
    this.control_gene.style.float = "left";
    control_div.appendChild(this.control_gene);

    var control_option = document.createElement('div');
    this.control_gene.style.margin = "10px";
    this.control_gene.style.float = "right";
    var control_option_display_lenght = document.createElement('input');
    control_option_display_lenght.setAttribute("type", "checkbox");
    control_option_display_lenght.setAttribute("id", "ch_length");
    control_option_display_lenght.setAttribute("name", "ch_length");
    control_option_display_lenght.setAttribute("checked", "true");
    control_option_display_lenght.setAttribute("value", "display");

    control_option_display_lenght.addEventListener("change", this.onChangeViewBranchLenght.bind(this));
    control_option.appendChild(control_option_display_lenght);
    var label_checkbox_tree_length = document.createElement('label');
    label_checkbox_tree_length.innerHTML = "Display branch lenght";
    control_option.appendChild(label_checkbox_tree_length);




    control_div.appendChild(control_option);
    this.container.appendChild(control_div);
    this.container.appendChild(this.treeview);

    this.zoom_lv = 3;

    this.active_names = [];
    this.display_branch_length = true;
    this.nodes=[];
    this.links=[];
    this.cell_size=40;
  }

  load(genelabel, species_tree, gene_tree) {

    //console.log(alignment);
    this.control_gene.innerHTML = genelabel;

    this.species_tree = species_tree;
    this.gene_tree = gene_tree;




    var newick_raw = species_tree;




    var newick_species = NewickTools.parse(newick_raw);
    let nodes_species = d3.hierarchy(newick_species, d => d.branchset);
    console.log(nodes_species);
    var stack = [];
    nodes_species.height = 0;
    stack.push(nodes_species);
    let count = 0;
    var num_leaf = 0;
    this.node_leaf = [];
    this.arr_sample_from_tree = [];
    this.nodes=[];
    this.links=[];
    var max_heigth = 0;
    var max_depth = 0;
    let map_nodes=new Map();
    while (stack.length > 0) {
      var n = stack.pop();
      n.id = count;
      count++;
      var newnode={ id: n.id, group: 0,type:1, label: n.data.name, level: n.depth };
      if (n.children != undefined) {
       
        map_nodes.set(n.id,[]);
        for (var i = 0; i < n.children.length; i++) {
          stack.push(n.children[i]);
          n.children[i].height = n.height + n.children[i].data.length;
        }
      } else {
        if (n.height > max_heigth) max_heigth = n.height;
        if (n.depth > max_depth) max_depth = n.depth;
        
        newnode.group=1;
        this.node_leaf.push(newnode);
        num_leaf++;
      }
      this.nodes.push( newnode);
      if(n.parent!=null){
        this.links.push({ target: n.parent.id,type:1, source: n.id, strength: 0.0 });
        let arr=map_nodes.get(n.parent.id);
        arr.push(newnode);
        map_nodes.set(n.parent.id,arr);
      }
        
    }
    var height=num_leaf*this.cell_size;
    this.container.style.height= (height + 100) + "px";
    var width_tree = this.props.width /2;
    var distance_per_depth = width_tree / max_depth;
    var diff=width_tree;
    var order=1;
    for (var i =0;i<this.nodes.length;i++){
      if(this.nodes[i].group==1){
        this.nodes[i].fx=diff;
        this.nodes[i].fy=order*this.cell_size;
        order=order+1;
      }
    }
    //console.log(map_nodes);
    
    for (var d=max_depth;d>=0;d--){
      for (var i =0;i<this.nodes.length;i++){
        //console.log(this.nodes[i]);
        if(this.nodes[i].group==0){
          if(this.nodes[i].level==d){
            var tb=0;
            for (var c=0;c<map_nodes.get(this.nodes[i].id).length;c++){   
              tb=tb+map_nodes.get(this.nodes[i].id)[c].fy;
            }
            //console.log(tb);  
            this.nodes[i].fy=tb/map_nodes.get(this.nodes[i].id).length;
            this.nodes[i].fx=distance_per_depth*d;
          }
        }
      }
    }
    
    //estimate length of sample name, by average length plus 10
    var newick_gene = NewickTools.parse(gene_tree);
    let nodes_gene = d3.hierarchy(newick_gene, d => d.branchset);
    console.log(nodes_gene);
    stack = [];
    var list_node_gene=[];
    nodes_gene.height = 0;
    stack.push(nodes_gene);
    
    var num_leaf = 0;
    
    var max_depth_gene = 0;
    map_nodes=new Map();
    while (stack.length > 0) {
      var n = stack.pop();
      n.id = count;
      count++;
      var newnode={ id: n.id, group: 0,type:2, label: n.data.name, level: n.depth };
      if (n.children != undefined) {
       
        map_nodes.set(n.id,[]);
        for (var i = 0; i < n.children.length; i++) {
          stack.push(n.children[i]);
          n.children[i].height = n.height + n.children[i].data.length;
        }
      } else {
        
        if (n.depth > max_depth_gene) max_depth_gene = n.depth;
        
        newnode.group=1;
       
       
      }
      list_node_gene.push( newnode);
      if(n.parent!=null){
        this.links.push({ target: n.parent.id, type:2,source: n.id, strength: 0.0 });
        let arr=map_nodes.get(n.parent.id);
        arr.push(newnode);
        map_nodes.set(n.parent.id,arr);
      }
    }
    distance_per_depth = width_tree / max_depth_gene;
    var order=1;
    for (var i =0;i<list_node_gene.length;i++){
      for(var j =0;j<this.node_leaf.length;j++){
        var sample_name=list_node_gene[i].label.substring(0,list_node_gene[i].label.lastIndexOf("_"));
        //console.log(sample_name);
        
        if(list_node_gene[i].group==1&&this.node_leaf[j].label==sample_name){
         // console.log("found match");
          list_node_gene[i].fx=this.node_leaf[j].fx;
          list_node_gene[i].fy=this.node_leaf[j].fy;
          order=order+1;
        }
      }
      
    }
    for (var d=max_depth_gene;d>=0;d--){
      for (var i =0;i<list_node_gene.length;i++){
        //console.log(this.nodes[i]);
        if(list_node_gene[i].group==0){
          if(list_node_gene[i].level==d){
            var tb=0;
            for (var c=0;c<map_nodes.get(list_node_gene[i].id).length;c++){   
              tb=tb+map_nodes.get(list_node_gene[i].id)[c].fy;
            }
            //console.log(tb);  
            list_node_gene[i].fy=tb/map_nodes.get(list_node_gene[i].id).length;
            list_node_gene[i].fx=2*(width_tree)-distance_per_depth*d;
          }
        }
      }
    }  
    
    
    //console.log(map_nodes);
    
  
    Array.prototype.push.apply(this.nodes, list_node_gene)
    console.log(this.nodes);
    // append the svg object to the body of the page
    // appends a 'group' element to 'svg'
    // moves the 'group' element to the top left margin
    //var height_rect = height / node_leaf.length / 2;
    this.make_square_view(this.nodes,this.links);

  }
  getChildNodes(parentNode, linkes){

      for (var i =0;i<this.links.length;i++){

      }
  }
  make_square_view(nodes, links){
    var square_links=[];
    var middle_nodes=[];
    for(var i=0;i<links.length;i++){
        var source_node;
        var target_node;
        for(var j=0;j<nodes.length;j++){
          if(links[i].source==nodes[j].id)
            source_node=nodes[j];
          if(links[i].target==nodes[j].id)
            target_node=nodes[j];        
        }
        var midder_node={ id: source_node.id+"_"+target_node.id, group: source_node.group,type:0, label: "", level: source_node.level };
        midder_node.fx=target_node.fx;
        midder_node.fy=source_node.fy;
        middle_nodes.push(midder_node);
        square_links.push({ target:midder_node.id , type:links[i].type,source: source_node.id, strength: 0.0 });
        square_links.push({ target:target_node, type:links[i].type,source: midder_node.id, strength: 0.0 });
    }
    this.links=square_links;

    Array.prototype.push.apply(nodes, middle_nodes)

  }
  logscale(min_value, max_value, max_scale, v) {
    if (v == 0) return 0;
    //change range to make sure all value greater 1
    let mul = Math.pow(10, -Math.log10(min_value));
    max_value = mul * max_value;
    v = mul * v;
    let k = max_scale / Math.log10(max_value);
    return k * Math.log10(v);
  }
  setOptions(options) {
    this.props.width = options.width;
    this.props.height = options.height;

  }
  setActiveNames(names) {
    this.active_names = names;
    this.drawHighlighTree();
  }

  draw() {
    
    var width =  this.props.width  ,
      height =  this.container.style.height;
    this.treeview.innerHTML = "";
    var svg = d3.select("#amr_gf_treeview").append("svg")
      .attr("width", width)
      .attr("height", height);
    console.log(svg);
    var linkForce = d3
      .forceLink()
      .id(function (link) { return link.id })
      .strength(function (link) { return link.strength })

    var simulation = d3
      .forceSimulation()
      .force('link', linkForce)
      .force('charge', d3.forceManyBody())
      .force('center', d3.forceCenter(width / 2, height / 2))

    var dragDrop = d3.drag().on('start', function (node) {
      node.fx = node.x
      node.fy = node.y
    }).on('drag', function (node) {
      simulation.alphaTarget(0.7).restart()
      node.fx = d3.event.x
      node.fy = d3.event.y
    }).on('end', function (node) {
      if (!d3.event.active) {
        simulation.alphaTarget(0)
      }
     
      node.fixed=true
    })
    this.linkElements = svg.append("g")
      .attr("class", "links")
      .selectAll("line")
      .data(this.links)
      .enter().append("line")
      .attr("stroke-width", this.getLinkStroke)
      .attr("stroke-linecap","round")
      .attr("stroke", this.getLinkColor)

      this.nodeElements = svg.append("g")
      .attr("class", "nodes")
      .selectAll("circle")
      .data(this.nodes)
      .enter().append("circle")
      .attr("r",this.getNodeSize)
      .attr("fill", this.getNodeColor)
      .call(dragDrop)
      .on('click', this.selectNode.bind(this))

      this.textElements = svg.append("g")
      .attr("class", "texts")
      .selectAll("text")
      .data(this.nodes)
      .enter().append("text")
      .text(function (node) { return node.label })
      .attr("font-size", 12)
      .attr("dx", -10)
      .attr("dy", 25)

    simulation.nodes(this.nodes).on('tick', () => {
      this.nodeElements
        .attr('cx', function (node) { return node.x })
        .attr('cy', function (node) { return node.y })
        this.textElements
        .attr('x', function (node) { return node.x })
        .attr('y', function (node) { return node.y })
        this.linkElements
        .attr('x1', function (link) { return link.source.x })
        .attr('y1', function (link) { return link.source.y })
        .attr('x2', function (link) { return link.target.x })
        .attr('y2', function (link) { return link.target.y })
    });
    simulation.force("link").links(this.links);
  }
  selectNode(selectedNode) {
   /*  
    var neighbors =this.getNeighborNodes(selectedNode);

    // we modify the styles to highlight selected nodes
    this.nodeElements.attr('fill', function (node) {  
      if (Array.isArray(neighbors) && neighbors.indexOf(node.id) > -1) {
        return node.level === 1 ? 'blue' : 'blue'
      }
  
      return node.level === 1 ? 'red' : 'red'
    })
    this.textElements.attr('fill', function (node) {
      return Array.isArray(neighbors) && neighbors.indexOf(node.id) > -1 ? 'green' : 'black' 
    })
   this.linkElements.attr('stroke', function (link) { 
      return link.target.id === selectedNode.id || link.source.id === selectedNode.id ? 'green' : '#E5E5E5'
    }) */
  }
  getNeighborNodes(node) {
   
    return this.links.reduce(function (neighbors, link) {
      if (link.target.id === node.id) {
        neighbors.push(link.source.id)
      } else if (link.source.id === node.id) {
        neighbors.push(link.target.id)
      }
      return neighbors
    },
      [node.id]
    );
  }

  isNeighborLink(node, link) {
    return link.target.id === node.id || link.source.id === node.id
  }


  getNodeColor(node, neighbors) {
    if (Array.isArray(neighbors) && neighbors.indexOf(node.id) > -1) {
      return node.level === 1 ? '#5233e8b' : 'green'
    }
    if(node.type==1)
      return node.group === 1 ? '#233e8b' : '#233e8b'
      if(node.type==2)
      return node.group === 1 ? '#cf0000' : '#cf0000'
  }

  getNodeSize(node, neighbors) {
  
    if(node.type==1)
      return node.group === 1 ? 10 : 12
      if(node.type==2)
      return node.group === 1 ? 8 : 4
      else
        return 1;
  }
  getLinkColor(link) {
    if(link.type==1) return  '#233e8b';
    else return  'red'
  }
  getLinkStroke(link) {
    if(link.type==1) return  '10';
    else return '2'
  }
  getTextColor(node, neighbors) {
    return Array.isArray(neighbors) && neighbors.indexOf(node.id) > -1 ? 'green' : 'black'
  }
  zoomIn() {

    //console.log("event zoom in");
    if (this.zoom_lv < 3) this.zoom_lv = this.zoom_lv + 1;
    this.draw();

  }
  zoomOut() {
    console.log("event zoom out");
    if (this.zoom_lv > 1) this.zoom_lv = this.zoom_lv - 1;
    this.draw();

  }

  onChangeViewBranchLenght() {
    this.display_branch_length = document.getElementById("ch_length").checked;
    this.load(this.genelabel, this.phylotree, this.samples);
    this.draw();
  }


}
export default GeneFlowTree
