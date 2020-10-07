/* eslint-disable */
import Phylocanvas from "phylocanvas";
import * as d3 from "d3";
export class Phylogeny {
  constructor(element) {

    this.container = element;
    this.props = {
      width: 800,
      height: 400,
      color: "#000"
    };
    var control_div = document.createElement('div');
    control_div.style.width = "100%";
    control_div.style.height = "40px";
    control_div.style.padding = "0px";
    control_div.style.margin = "0px";

    var control_select = document.createElement('div');
    control_select.style.margin = "10px";
    control_select.style.float = "right";
    this.type_select = document.createElement('select');
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
    control_div.appendChild(control_select);
    this.container.appendChild(control_div);
    var tree_div = document.createElement('div');
    tree_div.id="phy_tree";
    this.container.appendChild(tree_div);
    this.type="rectangular";
  }
  load(newick_tree) {
    this.newick_tree = newick_tree;


  }
  setOptions(options) {
    this.props.width = options.width;
    this.props.height = options.height;

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
    this.tree.load(this.newick_tree);


  }
  changeType(){
    this.type=this.type_select.value;
    this.draw();
  }
}
export default Phylogeny
