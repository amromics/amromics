/* eslint-disable */
import Chart from 'chart.js';
export class PanPieChart {
  constructor(element) {
    this.container = element;
    this.props = {
      width: 800,
      height: 400
    };
    var canvas = document.createElement('canvas');
    canvas.id="pie_canvas";
    canvas.width =   this.props.width;
    canvas.height = this.props.height;
    this.container.appendChild(canvas);

  }
  load(groups) {
    //console.log(groups);
    var list_numb=[]
      var list_label=[]

      for (var i =0;i<groups.length;i++){

          if(groups[i].num>0){
            list_numb.push(groups[i].num);
            list_label.push(groups[i].name+groups[i].des);
            }

      }
      this.data = {
        datasets: [{
          data: list_numb,
          backgroundColor: [
              'rgba(255, 99, 132, 0.2)',
              'rgba(54, 162, 235, 0.2)',
              'rgba(255, 206, 86, 0.2)',
              'rgba(75, 192, 192, 0.2)',
              'rgba(153, 102, 255, 0.2)',
              'rgba(255, 159, 64, 0.2)'
          ],

          borderColor: [
              'rgba(255, 99, 132, 1)',
              'rgba(54, 162, 235, 1)',
              'rgba(255, 206, 86, 1)',
              'rgba(75, 192, 192, 1)',
              'rgba(153, 102, 255, 1)',
              'rgba(255, 159, 64, 1)'
          ],
          borderWidth: 1
        }],

      labels: list_label
    };

  }
  setOptions(options) {
    this.props.width = options.width;
    this.props.height = options.height;
    this.container.innerHTML="";
    var canvas = document.createElement('canvas');
    canvas.id="pie_canvas";
    canvas.width =   this.props.width;
    canvas.height = this.props.height;
    this.container.appendChild(canvas);
  }
  draw() {
    var context=document.getElementById('pie_canvas');
    //console.log(this.data);
    var myPieChart = new Chart(context, {
        type: 'pie',
        data: this.data,
        options: {responsive: true}
      });
  }
}
export default PanPieChart
