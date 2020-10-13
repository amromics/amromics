/* eslint-disable */
import Chart from 'chart.js';
import * as ChartAnnotation from 'chartjs-plugin-annotation';
export class N50Chart{
  constructor(element) {
    this.container=element;
    this.input=undefined;
    this.list_data = [];
    this.list_labels = [];
    this.props={width:1200,height:400,title:'N50 chart',background_color:'rgba(255, 99, 132, 0.6)',line_color:'rgb(54, 162, 235)',x_label:'Num contigs',y_label:'Sum length'};
  }
  load(data){
    this.input=data;
    this.list_data = []
    this.list_labels = []

    this.sum = 0

    for (var i = 0; i < this.input.length; i++) {
        this.sum = this.sum + this.input[i].length;
        this.list_data.push(this.sum);
        this.list_labels.push(i + 1)
                        //console.log(sum)
    }
  }
  setOptions(options){
      this.props.width=options.width
      this.props.height=options.height
      if (options.title!=undefined)
        this.props.title=options.title
      if (options.background_color!=undefined)
        this.props.background_color=options.background_color
      if (options.line_color!=undefined)
        this.props.line_color=options.line_color
      if (options.x_label!=undefined)
        this.props.x_label=options.x_label
      if (options.y_label!=undefined)
        this.props.y_label=options.y_label
  }
  draw(){
    //var ctx = document.getElementById('n50Chart');
    //console.log(this.container);
    this.container.innerHTML = "";
    var canvas = document.createElement('canvas');
    canvas.id     = "cv_chart";
    canvas.width  = this.props.width;
    canvas.height = this.props.height;
    this.container.appendChild(canvas);
    var data = {
        labels: this.list_labels,
        datasets: [{
          label:'Sum of length',
          data: this.list_data,
          backgroundColor: [
              this.props.background_color
              ],
          borderColor:this.props.line_color,
          fill:false


            }]
    };
    var myLineChart = new Chart(canvas, {
        type: 'line',
        data: data,
        options: {
            responsive:true,
            title: {
              display: true,
              text: this.props.title
            },
            scales: {
                xAxes: [{
                  display: true,
                  scaleLabel: {
                        display: true,
                        labelString: this.props.x_label
                  }
                }],
                yAxes: [{
                  display: true,
                  scaleLabel: {
                      display: true,
                      labelString: this.props.y_label
                  }
                }]
            },
            annotation: {
              annotations: [
                {
                  type: "line",
                  mode: 'horizontal',
                  scaleID: 'y-axis-0',
                  borderColor: "red",
                  value:this.sum/2,
                  label: {
                    content: this.sum/2,
                    enabled: true,
                    
                  }
                }
              ]
            }
        }
    });
  }
}
export default N50Chart
