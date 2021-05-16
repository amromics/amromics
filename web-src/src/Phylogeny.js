/* eslint-disable */
import Phylocanvas from "phylocanvas";
import * as d3 from "d3";
export class Phylogeny {
  constructor(element) {

    this.container = element;
    this.props = {
      width: 800,
      height: 800,
      color: "#000"
    };
    var control_div = document.createElement('div');
    control_div.style.width = "100%";
    control_div.style.height = "40px";
    control_div.style.padding = "0px";
    control_div.style.margin = "0px";
    control_div.style.borderBottom="1px solid #aaa";
    var control_select = document.createElement('div');
    //control_select.style.margin = "10px";
    control_select.style.float = "right";
    this.type_select = document.createElement('select');
    this.type_select.style.border="1px solid #aaa";
    this.type_select.style.borderRadius="3px";
    this.type_select.style.padding="4px";
    this.type_select.style.backgroundColor="transparent";
    var opt_rect = document.createElement('option');
    opt_rect.appendChild(document.createTextNode("Rectangular"));
    opt_rect.value = "rectangular";
    this.type_select.appendChild(opt_rect);
    var opt_radial = document.createElement('option');
    opt_radial.appendChild(document.createTextNode("Radial"));
    opt_radial.value = "radial";
    this.type_select.appendChild(opt_radial);
    var opt_circular = document.createElement('option');
    opt_circular.appendChild(document.createTextNode("Circular"));
    opt_circular.value = "circular";
    this.type_select.appendChild(opt_circular);
    var opt_diagonal = document.createElement('option');
    opt_diagonal.appendChild(document.createTextNode("Diagonal"));
    opt_diagonal.value = "diagonal";
    this.type_select.appendChild(opt_diagonal);
    var opt_hierarchical = document.createElement('option');
    opt_hierarchical.appendChild(document.createTextNode("Hierarchical"));
    opt_hierarchical.value = "hierarchical";
    this.type_select.appendChild(opt_hierarchical);
    this.type_select.addEventListener("change", this.changeType.bind(this));
    //Event.observe(this.contig_select, 'change', changeContig.bind(this));
    control_select.appendChild(this.type_select);

    this.meta_select = document.createElement('select');
    this.meta_select.style.border="1px solid #aaa";
    this.meta_select.style.borderRadius="3px";
    this.meta_select.style.padding="4px";
    this.meta_select.style.backgroundColor="transparent";
    var opt_meta_select_default= document.createElement('option');
    opt_meta_select_default.appendChild(document.createTextNode("Select metadata"));
    opt_meta_select_default.value = "";
    this.meta_select.appendChild(opt_meta_select_default);
    control_div.appendChild(this.meta_select);
    control_div.appendChild(control_select);
    this.container.appendChild(control_div);
    var tree_div = document.createElement('div');
    tree_div.id="phy_tree";
    tree_div.style.height=this.props.height+"px";
    this.container.appendChild(tree_div);
    this.legend_container = document.createElement('div');
    this.container.appendChild(this.legend_container);
    this.type="rectangular";
  }
  load(newick_tree,metadata) {
    this.newick_tree = newick_tree;
    this.metadata=metadata;
    //console.log( this.metadata);
    //set default metadata style:
    this.columns=[];
    this.metadata_style={};
    var list_color_template=[
      "green",
      "red",
      "yellow",

      "peru",
      "chocolate",
      "pink",
      "#C83200",
      "#CD3700",
      "#FF6103",
      "#CC7722",
      "#FD6302",
      "#883000",
      "#FFBF00",
      "gold",
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
    ];
    for(var id in this.metadata){
        for (var column in this.metadata[id]){
          if (!this.columns.includes(column)){
            this.columns.push(column);
            this.metadata_style[column]={};
          }
          //console.log(this.metadata[i].metadata[column]);
          if(!(this.metadata[id][column] in this.metadata_style[column])){
           // console.log(this.metadata_style[column]);
            //console.log(this.metadata[i].metadata[column]);
            this.metadata_style[column][this.metadata[id][column]]='';
            if(this.metadata[id][column]!=""){
              this.metadata_style[column][this.metadata[id][column]]
                =list_color_template[Object.keys(this.metadata_style[column]).indexOf(this.metadata[id][column])];
              if(Object.keys(this.metadata_style[column]).indexOf(this.metadata[id][column])>list_color_template.length){
                this.metadata_style[column][this.metadata[id][column]]="#" + (Math.random().toString(16) + "000000").substring(2,8);
              }
            }
            else{
              this.metadata_style[column][this.metadata[id][column]]="black";
            }
          }
        }
    }
    //console.log(this.metadata_style);
    for(var column in this.columns){
      var opt = document.createElement('option');
      opt.appendChild(document.createTextNode(this.columns[column]));
      opt.value = this.columns[column];
      this.meta_select.appendChild(opt);
    }
    this.meta_select.addEventListener("change", this.changeMeta.bind(this));
  }
  setOptions(options,metadata_style) {
    this.props.width = options.width;
    this.props.height = options.height;
    this.metadata_style=metadata_style;

  }
  draw() {
    document.getElementById("phy_tree").innerHTML="";
    this.tree = Phylocanvas.createTree("phy_tree");
    this.tree.branchColour = this.props.color;
    this.tree.collapsedColour = this.props.color;

    //this.tree.setTreeType(this.type);
    this.tree.alignLabels = true;
    // this.tree.showLabels = false;
    this.tree.setNodeSize(20);
    this.tree.setTextSize(20);
    this.tree.lineWidth = 1;
    this.tree.setTreeType(this.type);
    console.log(this.newick_tree);
    ///////this.tree.showBranchLengthLabels =true;
    this.tree.showBootstrap = true;
    this.tree.showInternalNodeLabels = this.tree.showBootstrap;
	  this.tree.internalLabelStyle.colour = this.tree.branchColour;
	  this.tree.internalLabelStyle.font = this.tree.font;
	  this.tree.internalLabelStyle.textSize = this.tree.textSize/2;
    this.tree.load(this.newick_tree);
  //  this.tree.showBootstrap=true;
    //this.tree.displayLabels();
    //console.log(this.tree.leaves);


  }
  changeType(){
    this.type=this.type_select.value;
    this.draw();
    if(this.meta_select.value!=""){
      this.colorLeaves(this.meta_select.value);
    }
  }
  changeMeta(){
    var meta_selected=this.meta_select.value;
    this.colorLeaves(meta_selected);
  }
  colorLeaves(metadata_name){
    for(var i=0;i<this.tree.leaves.length;i++){
      var nodeid=this.tree.leaves[i].id.replace(/\'/g,'');
      //console.log(this.metadata[nodeid][metadata_name]);
      this.tree.leaves[i].setDisplay({
        colour:  this.metadata_style[metadata_name][this.metadata[nodeid][metadata_name]],
        shape: 'circle', // or square, triangle, star
        size: 1, // ratio of the base node size
        leafStyle: {
          strokeStyle: '#000000',
          fillStyle: this.metadata_style[metadata_name][this.metadata[nodeid][metadata_name]],
          lineWidth: 1,
        },
        labelStyle: {
          colour: this.metadata_style[metadata_name][this.metadata[nodeid][metadata_name]],
          textSize: 20, // points
          font: 'Arial',
          format: 'bold',
        },
      });
    }
    this.tree.draw();
    //add legend
    this.legend_container.innerHTML="";
    for (var value in this.metadata_style[metadata_name]){
      var box = document.createElement('div');
      box.style.width="20px";
      box.style.height="20px";
      box.style.float="left";
      box.style.marginRight="3px";
      box.style.backgroundColor=this.metadata_style[metadata_name][value];
      var text = document.createElement('div');
      text.innerHTML=(value=="")?"undefined":value;
      text.style.float="left";
      var legend=document.createElement('div');
      legend.style.float="left";
      legend.style.marginRight="20px";
      legend.style.marginTop="10px";
      legend.appendChild(box);
      legend.appendChild(text);
      this.legend_container.appendChild(legend);

    }
    this.legend_container.style.height="100px";
  }
}
export default Phylogeny
