/* eslint-disable */
import * as d3 from "d3";
export class Browser {
  constructor(element) {
    this.container = element;

    this.contigs = [];
    this.knowngenes = [];
    this.gc_skew = [];
    this.gc_content = undefined;
    this.gc_content_window = 1000;
    this.gc_content_step = 100;
    this.props = {
      width: 1200,
      highlight: [{
        type: "AMR",
        color: "#A71D31"
      }, {
        type: "VIR",
        color: "#0F5257"
      }],
      color: '#F3C366'
    };
    this.scale = 0.1;
    //create UI
    this.genUI();


  }
  genUI(){
    this.container.innerHTML = "";

    this.track_div = document.createElement('div');
    this.track_div.id = "track_div";
    this.track_div.style.width = (this.props.width) + "px";
    this.track_div.style.height = "500px";
    this.track_div.style.padding = "0px";
    this.track_div.style.margin = "0px";
    this.track_div.style.float = "left";
    this.track_div.style.overflowX = "scroll";
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
    var control_zoom_out = document.createElement('button');
    control_zoom_out.style.backgroundImage = "none";
    control_zoom_out.style.border = "1px solid transparent";
    control_zoom_out.style.borderRadius = "4px";
    control_zoom_out.textContent = "-";
    control_zoom_out.addEventListener("click", this.zoomOut.bind(this));

    var control_select = document.createElement('div');
    control_select.style.margin = "10px";
    control_select.style.float = "left";
    this.contig_select = document.createElement('select');
    this.contig_select.id="select_contig";
    this.contig_select.addEventListener("change", this.changeContigSelect.bind(this));
   
    //Event.observe(this.contig_select, 'change', changeContig.bind(this));
    control_select.appendChild(this.contig_select);
    control_zoom.appendChild(control_zoom_out);
    control_zoom.appendChild(control_zoom_in);
    control_div.appendChild(control_zoom);
    control_div.appendChild(control_select);

    var content_div = document.createElement('div');
    content_div.style.width = (this.props.width) + "px";
    content_div.style.height = "500px";
    content_div.style.padding = "0px";
    content_div.style.margin = "0px";
    //content_div.style.clear = "both";
    content_div.appendChild(this.track_div);
    var wrap_div = document.createElement('div');
    wrap_div.style.width = (this.props.width) + "px";
    wrap_div.style.height = "540px";
    wrap_div.style.border = "1px solid gray";
    wrap_div.appendChild(control_div);
    wrap_div.appendChild(content_div);
    this.container.appendChild(wrap_div);
  }
  load(contigs, genes, skew, content) {
    this.contigs = contigs;

    this.knowngenes = genes;
    this.gc_skew = skew;
    this.gc_content = content;
    //this.gc_content_window=content.window;
    //this.gc_content_step=content.step;
    console.log(contigs)
    for (var i = 0; i < this.contigs.length; i++) {

      var opt = document.createElement('option');
      opt.appendChild(document.createTextNode(this.contigs[i].name));
      opt.value = this.contigs[i].name;
      document.getElementById("select_contig").add(opt);
     
    }
    
   
    this.current_contig = this.contigs[0];

  }
  setOptions(options) {
    this.props.width = options.width

    if (options.highlight != undefined)
      this.props.highlight = options.highlight
    if (options.color != undefined)
      this.props.color = options.color
    this.genUI();
    this.load(this.contigs,this.knowngenes,this.gc_skew,this.gc_content);
  }

  draw() {
    //draw title

    this.drawContig();
    //draw knowngene
    this.drawKnowgenes();
    this.drawSkew();
    this.drawContent();
  }
  drawContig() {

    var track_h = 40;
    var min_length = this.props.width;
    var margin = {
      top: 20,
      right: 30,
      bottom: 0,
      left: 30
    };
    //width = 100 + margin.left + margin.right,
    //height = track_h + margin.top + margin.bottom;


    var svg = d3
      .select("#track_div")
      .append("svg")
      .attr("width", min_length > this.current_contig.length * this.scale ? min_length : this.current_contig.length * this.scale)
      .attr("height", track_h)

    var g_contig = svg.append("g")
      .attr("transform", "translate(" + 0 + "," + margin.top + ")");
    g_contig.append("rect")
      .attr("x", 0)
      .attr("y", 10)
      .attr("width", min_length)
      .attr("height", track_h - 20)
      .style("fill", "#EEE")
    // g_contig.append("path")
    // .attr("d","M" + 10 + "," + 0
    //    + "h" + (this.current_contig.length * this.scale - 10)
    //    + "a" + 10 + "," + 10 + " 0 0 1 " + 10 + "," + 10
    //    + "v" + (track_h - 20 - 2 * 10)
    //    + "a" + 10 + "," + 10 + " 0 0 1 " + -10 + "," + 10
    //    + "h" + (10 - this.current_contig.length * this.scale)
    //    + "z")
    //   .style("fill", "#AAA");
    g_contig.append("rect")
      .attr("x", 0)
      .attr("y", 0)
      .attr("width", this.current_contig.length * this.scale)
      .attr("height", track_h - 20)
      .attr("rx", 10)
      .style("stroke", function(d) {
        return d3.rgb("#999").darker();
      })
      .style("fill", "#AAA")
    g_contig.append("text")

      .attr("x", 10)
      .attr("y", 15)
      .text(this.current_contig.name)
    //draw x axisX
    var x = d3
      .scaleLinear()
      .range([0, this.current_contig.length * this.scale])
      .domain([100, this.current_contig.length])


    var axisX = d3.axisTop(x).ticks(this.current_contig.length * this.scale / 100);
    svg.append("g")
      .attr("transform", "translate(" + 0 + "," + margin.top + ")")
      .call(axisX);
  }
  drawKnowgenes() {
    //filter by current contigs
    var data = [];
    //console.log(this.current_contig.id);
    for (var i = 0; i < this.knowngenes.length; i++) {
      if (this.knowngenes[i].contig == this.current_contig.name) {
        data.push({
          contig: this.knowngenes[i].contig,
          start: this.knowngenes[i].start,
          end: this.knowngenes[i].end,
          name: this.knowngenes[i].name,
          strain: this.knowngenes[i].strain,
          type: this.knowngenes[i].type,
          product: this.knowngenes[i].product
        });
      }
    }
    var track_h = 40;
    var min_length = this.props.width;
    var margin = {
      top: 0,
      right: 0,
      bottom: 0,
      left: 30
    };
    var diff = 40



    var svg = d3
      .select("#track_div")
      .append("svg")
      .attr("width", min_length > this.current_contig.length * this.scale ? min_length : this.current_contig.length * this.scale)
      .attr("height", 150)


    var g_track = svg.append("g")
      .attr("transform", "translate(" + 0 + "," + margin.top + ")");
    g_track.append("text")

      .attr("x", 10)
      .attr("y", 15)
      .text("Genes");
    var tooltip = d3
      .select("#trackSVG")
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
          d.start +
          " ->" +
          d.end +
          ":" + d.gene +
          "<br>" +
          "Type: " + d.type +
          "<br/>Product:" + d.product



        )
        .style("left", d3.mouse(this)[0] + 70 + "px")
        .style("top", d3.mouse(this)[1] + "px");
    };
    var mouseleave = function(d) {
      tooltip.style("opacity", 0);
      d3.select(this).style("stroke", "white");
    };
    g_track.append("rect")
      .attr("x", 0)
      .attr("y", 0)
      .attr("width", min_length > this.current_contig.length * this.scale ? min_length : this.current_contig.length * this.scale)
      .attr("height", 2)
      .style("fill", "#EEE")
    //console.log(data);
    for (var i = 0; i < data.length; i++) {

      this.drawGeneBlock(i, g_track, data[i].start * this.scale, data[i].end * this.scale, data[i].strain, data[i].name, this.props.color, mouseover, mousemove, mouseleave);
    }
  }


  drawGeneBlock(i, g, pos_s, pos_e, direction, label, color, mouseover, mousemove, mouseleave) {
    var level = i % 4 + 1;
    if (direction == "+") {
      g.append("polygon")
        .attr("points", pos_s + "," + (10 + level * 20) + " " + (pos_s + 10) + "," + (20 + level * 20) + " " + pos_e + "," + (20 + level * 20) + " " + pos_e + "," + (0 + level * 20) + " " + (pos_s + 10) + "," + (0 + level * 20))
        .style("fill", color)
        .style("stroke", "black")
        .style("stroke-width", 1)

    } else {
      g.append("polygon")
        .attr("points", pos_s + "," + (0 + level * 20) + " " + pos_s + "," + (20 + level * 20) + " " + (pos_e - 10) + "," + (20 + level * 20) + " " + pos_e + "," + (10 + level * 20) + " " + (pos_e - 10) + "," + (0 + level * 20))
        .style("fill", color)
        .style("stroke", "black")
        .style("stroke-width", 1)

    }
    g.append("text")
      .attr('text-anchor', 'start')
      .attr("x", pos_e + 5)
      .attr("y", 15 + (level * 20))
      .text(label)
  }
  drawSkew() {
    var data = [];
    console.log(this.current_contig.name);
    for (var i = 0; i < this.gc_skew.length; i++) {
      if (this.gc_skew[i].contig.split(' ')[0] == this.current_contig.name)
        for (var j = 0; j < this.gc_skew[i].GC.length; j++) {
          data.push({
            contig: this.gc_skew[i].contig.split(' ')[0],
            position: j * 100 + 1000,
            value: this.gc_skew[i].GC[j]
          });
        }


    }
    var track_h = 120;
    var min_length = this.props.width;;
    var margin = {
      top: 30,
      right: 0,
      bottom: 0,
      left: 30
    };
    var svg = d3
      .select("#track_div")
      .append("svg")
      .attr("width", min_length > this.current_contig.length * this.scale ? min_length : this.current_contig.length * this.scale)
      .attr("height", 150)
    var g_track = svg.append("g")
      .attr("transform", "translate(" + 0 + "," + margin.top + ")");

    g_track.append("text")

      .attr("x", 10)
      .attr("y", -10)
      .text("GC skew");



    var pre_point = data[0].value;
    for (var i = 1; i < data.length; i++) {
      this.drawSkewLine(g_track, data[i].position * this.scale, 100 * this.scale, data[i].value, pre_point);
      pre_point = data[i].value;
    }
    //draw oy axis
    g_track.append("line")
      .attr("x1", 0)
      .attr("y1", 50)
      .attr("x2", (data.length + 10) * 100 * this.scale)
      .attr("y2", 50)

      .style("stroke", "#0C4524")
      .style("stroke-width", 1)
    g_track.append("line")
      .attr("x1", 0)
      .attr("y1", 0)
      .attr("x2", (data.length + 10) * 100 * this.scale)
      .attr("y2", 0)
      .style('stroke-dasharray', '5,3')
      .style("stroke", "#0C4524")
      .style("stroke-width", 1)
    g_track.append("line")
      .attr("x1", 0)
      .attr("y1", 100)
      .attr("x2", (data.length + 10) * 100 * this.scale)
      .attr("y2", 100)
      .style('stroke-dasharray', '5,3')
      .style("stroke", "#0C4524")
      .style("stroke-width", 1)
    g_track.append("text")

      .attr("x", 0)
      .attr("y", 10)
      .text("1");
    g_track.append("text")
      .attr("x", 0)
      .attr("y", 100)
      .text("-1");



  }
  drawSkewLine(g, pos_s, window_size, value, pre_value) {

    g.append("line")
      .attr("x1", pos_s - window_size)
      .attr("y1", 50 - pre_value * 50)
      .attr("x2", pos_s)
      .attr("y2", 50 - value * 50)
      .style("fill", "#0C4524")
      .style("stroke", "black")
      .style("stroke-width", 1)

  }
  drawContent() {
    var data = [];
    //console.log(this.gc_content);
    for (var i = 0; i < this.gc_content.length; i++) {
      if (this.gc_content[i].contig.split(' ')[0] == this.current_contig.name)
        for (var j = 0; j < this.gc_content[i].GC.length; j++) {
          data.push({
            contig: this.gc_content[i].contig.split(' ')[0],
            position: j * 100 + 1000,
            value: this.gc_content[i].GC[j]
          });
        }


    }
    var track_h = 120;
    var min_length = this.props.width;;
    var margin = {
      top: 30,
      right: 0,
      bottom: 0,
      left: 30
    };
    var svg = d3
      .select("#track_div")
      .append("svg")
      .attr("width", min_length > this.current_contig.length * this.scale ? min_length : this.current_contig.length * this.scale)
      .attr("height", 120)
    var g_track = svg.append("g")
      .attr("transform", "translate(" + 0 + "," + margin.top + ")");

    g_track.append("text")

      .attr("x", 10)
      .attr("y", -10)
      .text("GC content");
    //console.log(data);
    var pre_point = data[0].value;
    for (var i = 0; i < data.length; i++) {
      this.drawContentLine(g_track, data[i].position * this.scale, 100 * this.scale, data[i].value, pre_point);
      pre_point = data[i].value;
    }
    g_track.append("line")
      .attr("x1", 0)
      .attr("y1", 150)
      .attr("x2", (data.length + 10) * 100 * this.scale)
      .attr("y2", 150)
      .style("stroke", "#0C4524")
      .style("stroke-width", 1);
    g_track.append("line")
      .attr("x1", 0)
      .attr("y1", 0)
      .attr("x2", (data.length + 10) * 100 * this.scale)
      .attr("y2", 0)
      .style('stroke-dasharray', '5,3')
      .style("stroke", "#0C4524")
      .style("stroke-width", 1);
    g_track.append("line")
      .attr("x1", 0)
      .attr("y1", 75)
      .attr("x2", (data.length + 10) * 100 * this.scale)
      .attr("y2", 75)
      .style('stroke-dasharray', '5,3')
      .style("stroke", "#0C4524")
      .style("stroke-width", 1);
    g_track.append("text")
      .attr("x", 0)
      .attr("y", 75)
      .text("50%");
    g_track.append("text")
      .attr("x", 0)
      .attr("y", 10)
      .text("100%");
  }
  drawContentLine(g, pos_s, window_size, value, pre_value) {

    g.append("line")
      .attr("x1", pos_s - window_size)
      .attr("y1", 150 - pre_value * 150)
      .attr("x2", pos_s)
      .attr("y2", 150 - value * 150)
      .style("fill", "#0C4524")
      .style("stroke", "black")
      .style("stroke-width", 1)

  }
  drawAgain() {
    this.track_div.innerHTML = "";
    this.drawContig();
    this.drawKnowgenes();
    this.drawSkew();
    this.drawContent();
  }
  //var _this=this;
  changeContigSelect() {
    
    if (this.contigs != undefined) {
      for (var i = 0; i < this.contigs.length; i++) {

        if (this.contigs[i].name == this.contig_select.value) {

          this.current_contig = this.contigs[i];
          break;
        }
      }
      console.log(this.current_contig);
      this.drawAgain();
    }
  }
  changeContig(contig) {
    for (var i = 0; i < this.contigs.length; i++) {

      if (this.contigs[i].name == contig) {

        this.current_contig = this.contigs[i];
        break;
      }
    }

    //this.current_contig = contig;
    console.log(this.current_contig);
    this.contig_select.value=contig;
    this.drawAgain();

  }
  zoomIn() {
    console.log(this.zoomOut);
    //console.log("event zoom in");
    if (this.scale < 1)
      this.scale = this.scale + 0.01;
    this.drawAgain();

  }
  zoomOut() {
    if (this.scale > 0.01)
      this.scale = this.scale - 0.01;
    this.drawAgain();

  }

}


export default Browser
