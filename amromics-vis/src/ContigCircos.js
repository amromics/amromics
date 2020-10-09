/* eslint-disable */
import Circos from 'circos';
export class ContigCircos {
  constructor(element) {
    this.container = element;
    this.container_id=this.container.id;
    this.props = {
      width: 300,
      height: 300
    };



  }
  load(contig,amr_genes, virulome_genes, skew) {
    this.contigs=contig;
    //console.log("skew");
    //console.log(skew);
    var set_highlight_contigs = new Set();
    this.amr_data = [];
    for (var i = 0; i < amr_genes.length; i++) {
      this.amr_data.push({
        block_id: amr_genes[i].sequence,
        start: amr_genes[i].start,
        end: amr_genes[i].end,
        label: amr_genes[i].gene
      });
      set_highlight_contigs.add(amr_genes[i].sequence);
    }
    this.vir_data = [];
    for (var i = 0; i < virulome_genes.length; i++) {
      this.vir_data.push({
        block_id: virulome_genes[i].sequence,
        start: virulome_genes[i].start,
        end: virulome_genes[i].end,
        label: virulome_genes[i].gene
      });
      set_highlight_contigs.add(virulome_genes[i].sequence);
    }
    this.skew_data = [];

    for (var i = 0; i < skew.length; i++) {
      if (set_highlight_contigs.has(skew[i].contig.split(' ')[0])){
        for (var j = 0; j < skew[i].GC.length; j++) {
          this.skew_data.push({
            block_id: skew[i].contig.split(' ')[0],
            position: j * 100 + 1000,
            value: skew[i].GC[j]
          });
        }
      }
    }

    this._contig = [];
    for (var i = 0; i < this.contigs.length; i++) {
      var contig_id = this.contigs[i].name.split(' ')[0];
      if (set_highlight_contigs.has(contig_id))
        this._contig.push({
          len: this.contigs[i].length,
          color: "#bc80bd",
          label: contig_id.split('_')[2],
          id: contig_id
        });
    }

  }
  setOptions(options) {
    this.props.width = options.width
    this.props.height = options.height
  }
  draw() {
    var container=this.container;
    this._circos = new Circos({
      container: '#' + this.container.id,
      width: this.props.width,
      height: this.props.height,
    });
    var unit=this.props.width/100;
    var configuration = {
      innerRadius: unit*40,
      outerRadius: unit*45,
      cornerRadius: 0,
      gap: 0.04, // in radian
      labels: {
        display: true,
        position: 'center',
        size: '10px',
        color: '#000000',
        radialOffset: 20,
      },
      ticks: {
        display: true,
        color: 'grey',
        spacing: 10000000,
        labels: true,
        labelSpacing: 10,
        labelSuffix: 'Mb',
        labelDenominator: 1000000,
        labelDisplay0: true,
        labelSize: '10px',
        labelColor: '#000000',
        labelFont: 'default',
        majorSpacing: 5,
        size: {
          minor: 2,
          major: 5,
        }
      },
      tooltipContent: function(d) {
        return `${d.label}`
      },
      events: {
        'click.alert': function (datum, index, nodes, event) {
            container.dispatchEvent(new CustomEvent("contig_select", {
            detail: datum.id
          }));
        }
      }
    }
    var amr_configuration = {
      innerRadius: 0.70,
      outerRadius: 0.8,
      min: null,
      max: null,
      color: '#fd6a62',
      strokeColor: '#fd6a62',
      strokeWidth: 1,
      direction: 'out',
      thickness: unit*10,
      radialMargin: 2,
      margin: 2,
      opacity: 1,
      logScale: false,
      tooltipContent: function(d) {
        return `${d.label}`
      },
      events: {}
    }
    var virulome_configuration = {
      innerRadius: 0.55,
      outerRadius: 0.65,
      min: null,
      max: null,
      color: '#bebada',
      strokeColor: '#bebada',
      strokeWidth: 1,
      direction: 'out',
      thickness: unit*10,
      radialMargin: 2,
      margin: 2,
      opacity: 1,
      logScale: false,
      tooltipContent: function(d) {
        return `${d.label}`
      },
      events: {}
    }
    var skew_configuration = {
      innerRadius: 0.8,
      outerRadius: 0.95,
      min: null,
      max: null,
      color: 'black',
      strokeColor: 'black',
      strokeWidth: 1,
      direction: 'out',



      opacity: 1,
      logScale: true,
      tooltipContent: function(d) {
        return `${d.GC}`
      },
      axes: [{
        color: 'blue',
        position: 0.00000001,
        thickness: 1, // in pixel
        opacity: 0.5 // between 0 and 1
      }],
      events: {}
    }

    this._circos.layout(this._contig, configuration);
    this._circos.stack('stack1', this.amr_data, amr_configuration);
    this._circos.stack('stack2', this.vir_data, virulome_configuration);

    this._circos.line('line1', this.skew_data, skew_configuration);
    //console.log(this.skew_data);
    this._circos.render();
  }
}
export default ContigCircos
