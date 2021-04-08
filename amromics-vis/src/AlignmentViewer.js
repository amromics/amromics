/* eslint-disable */
import * as d3 from "d3";
import {
  NewickTools
} from "newick";
export class AlignmentViewer {
  constructor(element) {
    this.container = element;
    this.props = {
      width: this.container.clientWidth,
      height: 400
    };
    this.treeview = document.createElement("div");
    this.treeview.id = "al_treeview";
    this.treeview.style.width = this.props.width / 5 + "px";
    this.treeview.style.height = this.props.height + "px";
    this.treeview.style.float = "left";
    this.alignmentview = document.createElement("div");
    this.alignmentview.id = "al_alignment";
    this.alignmentview.style.width = (this.props.width - this.props.width / 5) + "px";
    this.alignmentview.style.height = this.props.height + "px";
    this.alignmentview.style.float = "left";
    this.alignmentview.style.overflowX = "scroll";

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

    var control_option_radio_nucl = document.createElement('input');
    control_option_radio_nucl.setAttribute("type", "radio");
    control_option_radio_nucl.setAttribute("name", "view");
    control_option_radio_nucl.setAttribute("value", "nucl");
    control_option_radio_nucl.setAttribute("checked", "true");
    control_option_radio_nucl.addEventListener("change", this.onChangeType.bind(this));
    var label_radio_button_nucl = document.createElement('label');
    label_radio_button_nucl.innerHTML = "Nucleotide";
    control_option.appendChild(control_option_radio_nucl);
    control_option.appendChild(label_radio_button_nucl);
    var control_option_radio_prot = document.createElement('input');
    control_option_radio_prot.setAttribute("type", "radio")
    control_option_radio_prot.setAttribute("name", "view")
    control_option_radio_prot.setAttribute("value", "prot")
    control_option_radio_prot.addEventListener("change", this.onChangeType.bind(this));
    var label_radio_button_prot = document.createElement('label');
    label_radio_button_prot.innerHTML = "Protein";
    control_option.appendChild(control_option_radio_prot);
    control_option.appendChild(label_radio_button_prot);



    control_div.appendChild(control_option);
    this.container.appendChild(control_div);
    this.container.appendChild(this.treeview);
    this.container.appendChild(this.alignmentview);
    this.zoom_lv = 3;
    this.current_type="nucl";
    this.active_names=[];
    this.display_branch_length=true;
  }

  load(genelabel, phylotree, alignment,type,mlst) {

    //console.log(alignment);
    this.control_gene.innerHTML = genelabel;
    this.samples = alignment;
    this.phylotree=phylotree;
    this.genelabel=genelabel;

    this.mlst=mlst;

    this.pos = [];
    this.data = [];
    var list_sample = [];
    var maxpos = 0;
    var taxa = [];

    if(this.current_type=='nucl'){
      for (var i = 0; i < this.samples.length; i++) {
        //console.log(this.samples[i].seq);
        var arr = [...this.samples[i].seq];
        list_sample.push(this.samples[i].sample);
        taxa.push({
          name: this.samples[i].sample
        });
        maxpos = arr.length > maxpos ? arr.length : maxpos;
        for (var j = 0; j < arr.length; j++) {
          this.data.push({
            sample: this.samples[i].sample,
            pos: j + 1,
            value: arr[j],
            type: 'nucl'
          });
        }
      }
    }
    else{
      for (var i = 0; i < this.samples.length; i++) {
      //  var prot_seq=this.translateDNA2Prot(this.samples[i].seq);
      //  console.log(prot_seq);
        var prot_seq=this.samples[i].protein;
        var arr = [...prot_seq];
        list_sample.push(this.samples[i].sample);
        taxa.push({
          name: this.samples[i].sample
        });
        maxpos = arr.length > maxpos ? arr.length : maxpos;
        for (var j = 0; j < arr.length; j++) {
          this.data.push({
            sample: this.samples[i].sample,
            pos: j + 1,
            value: arr[j],
            type: 'prot'
          });
        }
      }
    }
    //console.log(this.data);
    for (var i = 1; i <= maxpos; i++) {
      this.pos.push(i);
    }
    var cell_w = 20;
    var cell_h = 30;
    var margin = {
        top: 30,
        right: 30,
        bottom: 30,
        left: 30
      },
      width = cell_w * this.pos.length + margin.left + margin.right,
      height = cell_h * list_sample.length + margin.top + margin.bottom;

    var newick_raw = phylotree;




    var newick = NewickTools.parse(newick_raw);
  //  console.log(newick);
    //estimate length of sample name, by average length plus 10
    var numchar=0;
    var listsamples=Object.keys(NewickTools.dfs(newick));
    //console.log(NewickTools.dfs(newick));
    for (var s in listsamples){
      numchar=numchar+listsamples[s].length;
      //console.log(listsamples[s]);
    }
    var namelength=numchar/listsamples.length*5+10;
    var width_tree = this.props.width / 5 - margin.left - margin.right-namelength;
    //console.log(newick);
    // declares a tree layout and assigns the size
    const treemap = d3.tree().size([height, width_tree]);

    //  assigns the data to a hierarchy using parent-child relationships
    let nodes = d3.hierarchy(newick, d => d.branchset);
    //console.log(nodes);
    // maps the node data to the tree layout
    let dist = 20;

    var stack = [];
    nodes.height = 0;
    stack.push(nodes);
    var max_heigth = 0;
    var min_height_not_zero=1;
    var max_depth = 0;
    //travel to cal height of nodes
    let count = 0;
    var num_leaf = 0;
    this.node_leaf = [];
    this.arr_sample_from_tree = [];
    while (stack.length > 0) {
      var n = stack.pop();
      n.id = count;
      //console.log(n);
      count++;
      if (n.children != undefined) {
        for (var i = 0; i < n.children.length; i++) {
          stack.push(n.children[i]);
          n.children[i].height = n.height + n.children[i].data.length;
        }
      } else {
        if (n.height > max_heigth) max_heigth = n.height;
        if (n.height < min_height_not_zero && n.height>0) min_height_not_zero = n.height;
        if (n.depth > max_depth) max_depth = n.depth;
        this.node_leaf.push(n);
        this.arr_sample_from_tree.push(
          n.data.name.replace(/\'/g, "")
        );
        num_leaf++;
      }
    }
    //console.log("arr from tree");
    //console.log(arr_sample_from_tree);
    // edit leaft position y (x coordinate):
    var distance_per_depth = width_tree / max_depth;
    var unit_distance = width_tree / max_heigth;
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
        n.ax = count * cell_h;
        count++;
        //setup node by depth
        n.ay = max_depth * distance_per_depth;

      }
      //setup node by heigth
      //n.ay = n.height * unit_distance;
      //console.log(this.display_branch_length);
      if(this.display_branch_length)
        n.ay = this.logscale(min_height_not_zero,max_heigth,width_tree,n.height);
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
    ///console.log(fnodes);
    // append the svg object to the body of the page
    // appends a 'group' element to 'svg'
    // moves the 'group' element to the top left margin
    //var height_rect = height / node_leaf.length / 2;

  }
  logscale(min_value,max_value,max_scale,v){
    if (v==0) return 0;
    //change range to make sure all value greater 1
    let mul=Math.pow(10,-Math.log10(min_value));
    max_value=mul*max_value;
    v=mul*v;
    let k=max_scale/Math.log10(max_value);
    return k*Math.log10(v);
  }
  setOptions(options) {
    this.props.width = options.width;
    this.props.height = options.height;

  }
  setActiveNames(names){
    this.active_names=names;
    this.drawHighlighTree();
  }

  draw() {
    var cell_w = 20;
    var cell_h = 30;
    var step_x = 2;
    if (this.zoom_lv == 2) {
      cell_w = 10;
      step_x = 5;
    }
    if (this.zoom_lv == 1) {
      cell_w = 5;
      step_x = 10;
    }
    var margin = {
        top: 30,
        right: 30,
        bottom: 30,
        left: 30
      },
      width = cell_w * this.pos.length + margin.left + margin.right,
      height = cell_h * this.arr_sample_from_tree.length + margin.top + margin.bottom;

    this.alignmentview.innerHTML = "";
    this.alignmentview.style.height = (height + 20) + "px";
    this.treeview.style.height = (height + 20) + "px";
    this.container.style.height= (height + 80) + "px";
    document.getElementById(this.treeview.id).innerHTML = "";
    var svg2 = d3
      .select("#" + this.treeview.id)
      .append("svg")
      .attr("width", this.props.width / 5)
      .attr("height", cell_h * this.node_leaf.length + margin.top + margin.bottom);
    var g = svg2
      .append("g")
      .attr(
        "transform",
        "translate(" +
        margin.left +
        "," +
        (margin.top + cell_h / 2) +
        ")"
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
    var y_data=this.arr_sample_from_tree.slice();
    y_data.reverse();
    var y = d3
      .scaleBand()
      .range([cell_h * this.arr_sample_from_tree.length, 0])
      .domain(y_data)
      .padding(0.01);
    // create a tooltip

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

    var svg = d3
      .select("#" + this.alignmentview.id)
      .append("svg")
      .attr("width", width)
      .attr("height", height)
      .append("g")
      .attr("transform", "translate(" + 0 + "," + margin.top + ")");

    //console.log(data);
    // Build color scale

    // Build X scales and axis:
    var x = d3
      .scaleBand()
      .range([0, cell_w * this.pos.length])
      .domain(this.pos)
      .padding(0.01);
    var axisX = d3.axisTop(x);

    axisX.tickValues(
      x.domain().filter(function(d, i) {
        if (i == 0) return i + 1;
        else {
          if ((i + 1) % step_x == 0) {
            return i;
          } else {
            return "";
          }
        }
      })
    );
    svg
      .append("g")
      .attr("transform", "translate(" + 0 + "," + 0 + ")")
      .call(axisX);

    // Build X scales and axis:

    //svg.append("g").call(d3.axisLeft(y));
    //console.log(this.arr_sample_from_tree);
      //console.log(this.phylotree);

    // var y = d3
    //   .scaleBand()
    //   .range([cell_h * this.arr_sample_from_tree.length, 0])
    //   .domain(y_data)
    //   .padding(0.01);
    //console.log(this.current_type);
    // add the squares
    svg
      .selectAll()
      .data(this.data, function(d) {
        return d.sample + ":" + d.pos;
      })
      .enter()
      .append("rect")
      .attr("x", function(d) {
        return x(d.pos);
      })
      .attr("y", function(d) {
        return y(d.sample);
      })
      .attr("width", cell_w)
      .attr("height", cell_h)
      .style("fill", function(d) {
        if (this.zoom_lv == 1) {
          if (ref_chars[d.pos - 1] == d.value && d.sample != "REF")
            return "#FFFFFF";
        }

        if (d.type == "nucl") {

          if (d.value == "A") return "#A0A0FF";
          if (d.value == "C") return "#FF8C4B";
          if (d.value == "T") return "#FF7070";
          if (d.value == "G") return "#A0FFA0";
          if (d.value == "U") return "#B8B8B8";
          if (d.value == "-") return "#FFFFFF";
        } else {
          if (d.value == "A") return "#C8C8C8";
          if (d.value == "D" || d.value == "E") return "#E60A0A";
          if (d.value == "C" || d.value == "M") return "#E6E600";
          if (d.value == "K" || d.value == "R") return "#145AFF";
          if (d.value == "S" || d.value == "T") return "#FA9600";
          if (d.value == "F" || d.value == "Y") return "#3232AA";
          if (d.value == "N" || d.value == "Q") return "#00DCDC";
          if (d.value == "G") return "#EBEBEB";
          if (d.value == "L" || d.value == "V" || d.value == "I")
            return "#0F820F";
          if (d.value == "W") return "#B45AB4";
          if (d.value == "H") return "#8282D2";
          if (d.value == "P") return "#DC9682";
          if (d.value == "-") return "#FFFFFF";
        }
        return "#000";
      });
    //draw text with zoom=3 only
    if (this.zoom_lv == 3) {
      svg
        .selectAll()
        .data(this.data, function(d) {
          return d.sample + ":" + d.pos;
        })
        .enter()
        .append("text")
        .attr("dy", ".35em")
        .attr("x", function(d) {
          return x(d.pos) + cell_w / 4;
        })
        .attr("y", function(d) {
          return y(d.sample) + cell_h / 2;
        })
        .style("color", "black")
        .text(d => d.value);
    }
  }
  drawHighlighTree(){
    var margin = {
      top: 30,
      right: 30,
      bottom: 30,
      left: 30
    };
    var cell_h = 30;
    document.getElementById(this.treeview.id).innerHTML = "";
    var svg2 = d3
      .select("#" + this.treeview.id)
      .append("svg")
      .attr("width", this.props.width / 5)
      .attr("height", cell_h * this.node_leaf.length + margin.top + margin.bottom);
    var g = svg2
      .append("g")
      .attr(
        "transform",
        "translate(" +
        margin.left +
        "," +
        (margin.top + cell_h / 2) +
        ")"
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
    var y_data=this.arr_sample_from_tree.slice();
    y_data.reverse();
    if (this.mlst!=undefined){
      //look for mslt for each sample
      for (var i =0;i<y_data.length;i++)
        for (var j=0;j<this.mlst.length;i++)
          if(y_data[i]==this.mlst[j])
            y_data[i]=y_data[i]+"("+this.mlst[j]+")";
    }
    var y = d3
        .scaleBand()
        .range([cell_h * this.arr_sample_from_tree.length, 0])
        .domain(y_data)
        .padding(0.01);
      // create a tooltip

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
  onChangeType(){
    this.current_type=document.querySelector('input[name="view"]:checked').value;
    this.load(this.genelabel,this.phylotree,this.samples );
    this.draw();
  }
  onChangeViewBranchLenght(){
    this.display_branch_length=document.getElementById("ch_length").checked;
    this.load(this.genelabel,this.phylotree,this.samples);
    this.draw();
  }
  translateDNA2Prot(dna_seq) {
    // your original data renamed to aminoDict
    const aminoDict = {
      'A': ['GCA', 'GCC', 'GCG', 'GCT'],
      'C': ['TGC', 'TGT'],
      'D': ['GAC', 'GAT'],
      'E': ['GAA', 'GAG'],
      'F': ['TTC', 'TTT'],
      'G': ['GGA', 'GGC', 'GGG', 'GGT'],
      'H': ['CAC', 'CAT'],
      'I': ['ATA', 'ATC', 'ATT'],
      'K': ['AAA', 'AAG'],
      'L': ['CTA', 'CTC', 'CTG', 'CTT', 'TTA', 'TTG'],
      'M': ['ATG'],
      'N': ['AAC', 'AAT'],
      'P': ['CCA', 'CCC', 'CCG', 'CCT'],
      'Q': ['CAA', 'CAG'],
      'R': ['AGA', 'AGG', 'CGA', 'CGC', 'CGG', 'CGT'],
      'S': ['AGC', 'AGT', 'TCA', 'TCC', 'TCG', 'TCT'],
      'T': ['ACA', 'ACC', 'ACG', 'ACT'],
      'V': ['GTA', 'GTC', 'GTG', 'GTT'],
      'W': ['TGG'],
      'Y': ['TAC', 'TAT'],
      '':['TAA','TAG','TGA']
    };


    // codon dictionary derived from aminoDict
    let codonDict = {}
    for (let k of Object.keys(aminoDict))
      for (let a of aminoDict[k])
        codonDict[a] = k

    var result = "";
    dna_seq=dna_seq.replace(/\-/g,'');
    for (let i = 0; i < dna_seq.length; i += 3){
      //console.log(dna_seq.substr(i, 3)+"->"+codonDict[dna_seq.substr(i, 3)]);
      var codon=dna_seq.substr(i, 3);
      //console.log(codon.length);
      if(codon.length==3){
        result = result.concat(codonDict[codon]);
        //console.log(result);

      }



    }

    return result;
  }
}
export default AlignmentViewer
