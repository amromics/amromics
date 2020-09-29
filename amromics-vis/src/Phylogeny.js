/* eslint-disable */
import Phylocanvas from "phylocanvas";
import * as d3 from "d3";
export class Phylogeny {
  constructor(element) {
    console.log(element);
    this.container = element;
    this.props = {
      width: 800,
      height: 400,
      color:"#000"
    };
  
  }
  load(newick_tree) {
    this.newick_tree = newick_tree;


  }
  setOptions(options) {
    this.props.width = options.width;
    this.props.height = options.height;

  }
  draw() {
    this.tree = Phylocanvas.createTree(this.container.id);
    this.tree.branchColour = this.props.color;
    this.tree.collapsedColour = this.props.color;
    //this.tree.setTreeType(this.type);
    this.tree.alignLabels = true;
    // this.tree.showLabels = false;
    this.tree.setNodeSize(20);
    this.tree.setTextSize(20);
    this.tree.lineWidth = 1;
    this.tree.setTreeType('rectangular');
    this.tree.load(this.newick_tree);
  }
}
export default Phylogeny
