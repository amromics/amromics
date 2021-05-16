/* eslint-disable */
import Chart from 'chart.js';
export class DistributionChart {
  constructor(element) {
    this.container = element;
    this.props = {
      width: 800,
      height: 400
    };
    var canvas = document.createElement('canvas');
    canvas.id="distribution_canvas";
    canvas.width =   this.props.width;
    canvas.height = this.props.height;
    this.container.appendChild(canvas);

  }
  load(clusters) {
    //console.log(clusters);
    var list_data=[];
        var list_label=[];

        for (var i =0;i<clusters.length;i++){

          list_data.push(parseInt(clusters[i].noisolates));
          list_label.push(clusters[i].id)
        }
        this.data = {
          datasets: [{
            label:"Gene count distribution",
            data: list_data,
            backgroundColor: [
                'rgba(255, 99, 132, 0.6)'

            ],
            borderColor:[
                'rgba(255, 99, 132, 0.2)'

            ]


          }],

        labels: list_label
      };

  }
  setOptions(options) {
    this.props.width = options.width;
    this.props.height = options.height;
    this.container.innerHTML="";
    var canvas = document.createElement('canvas');
    canvas.id="distribution_canvas";
    canvas.width =   this.props.width;
    canvas.height = this.props.height;
    this.container.appendChild(canvas);
  }
  draw() {
    var context=document.getElementById('distribution_canvas');
    //console.log(this.data);
    var myLineChart = new Chart(context, {
          type: 'line',
          data: this.data,
          options:{ responsive: true}
        });
  }
}
export default DistributionChart
