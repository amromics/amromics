<template>
  <div>
     <div class="row">
         <div class="col-12">
       <h1>AMRVis: A pipelne for baterial genome analysis and visualization</h1>
      </div>
    </div>
    <!--Stats cards-->
    <div class="row">
      <div class="col-md-6 col-xl-6">
        <stats-card>
          <div class="icon-big text-center" :class="`icon-${statsCollection.type}`" slot="header">
            <i :class="statsCollection.icon"></i>
          </div>
          <div class="numbers" slot="content">
            <p>{{statsCollection.title}}</p>
            {{statsCollection.value}}
          </div>
          <div class="stats" slot="footer">
            <i :class="statsCollection.footerIcon"></i> {{statsCollection.footerText}}
          </div>
        </stats-card>
      </div>
      <div class="col-md-6 col-xl-6">
        <stats-card>
          <div class="icon-big text-center" :class="`icon-${statsSamples.type}`" slot="header">
            <i :class="statsSamples.icon"></i>
          </div>
          <div class="numbers" slot="content">
            <p>{{statsSamples.title}}</p>
            {{statsSamples.value}}
          </div>
          <div class="stats" slot="footer">
            <i :class="statsSamples.footerIcon"></i> {{statsSamples.footerText}}
          </div>
        </stats-card>
      </div>
    </div>

    <!--Charts-->
    <div class="row">
      <div class="col-12">
        <card title="Pipeline"
                    sub-title="Diagram">
        <img class="center_img" :src="`${publicPath}static/pipeline.png`"/>
        
        </card>
      </div>

      <div class="col-md-6 col-12">
        <chart-card title="Species Statistics"
                    sub-title="Portion of genus"
                    v-if="isReady"
                    :chart-data="genusChart.data"
                    chart-type="Pie">
        
         
        </chart-card>
      </div>

      <div class="col-md-6 col-12">
        <chart-card title="Timeline" v-if="isReady"
                    sub-title="Number of collections analyis in each month"
                    :chart-data="activityChart.data"
                    :chart-options="activityChart.options">
          <span slot="footer">
            <i class="ti-check"></i> Data information certified
          </span>
         
        </chart-card>
      </div>

    </div>

  </div>
</template>
<script>
import { StatsCard, ChartCard } from "@/components/index";
import Chartist from 'chartist';
import SampleAPI from "@/api/SampleAPI";
export default {
  components: {
    StatsCard,
    ChartCard
  },
  /**
   * Chart data used to render stats, charts. Should be replaced with server data
   */
  data() {
    return {
      publicPath: process.env.BASE_URL,
      statsCollection:{
        type: "warning",
          icon: "ti-folder",
          title: "Collections",
          value: "1",
          footerText: "last update :09-05-2021",
          footerIcon: "ti-timer"
      },
      statsSamples:{
          type: "success",
          icon: "ti-file",
          title: "Samples",
          value: "1",
          footerText: "last update :09-05-2021",
          footerIcon: "ti-timer"
      },
      activityChart: {
        data: {
          labels: [
            "Jan",
            "Feb",
            "Mar",
            "Apr",
            "Mai",
            "Jun",
            "Jul",
            "Aug",
            "Sep",
            "Oct",
            "Nov",
            "Dec"
          ],
          series: [
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
           
          ]
        },
        options: {
          seriesBarDistance: 10,
          axisX: {
            showGrid: false
          },
          height: "245px"
        }
      },
      genusChart: {
        data: {
          labels: [],
          series: []
        },
        options: {}
      },
      isLoading: true,
      list_collection: [],
      isReady: false,
      list_samples:[]
    };
  },
  computed: {},
  async created() {
    this.loading = true;
    await Promise.all([this.fetchData()]);
    this.loadData();
  },
  methods: {
    async fetchData() {
      // const value = await CollectionResult.fetchResult()
      // console.log("below is samle id")
      // console.log(this.sampleId)
      this.isReady = false;
      const value = await SampleAPI.fetchListCollectionResult();
      //console.log(result)
      //var result = value.data.results;
      this.list_collection = value.data.collections;
       for (var i=0;i<this.list_collection.length;i++){
        const samples_value = await SampleAPI.fetchSetResult(this.list_collection[i].collectionID);
        console.log(samples_value);
        for(var j=0;j<samples_value.data.samples.length;j++){
            var obj=samples_value.data.samples[j];
            obj['collectionID']=this.list_collection[i].collectionID;
            this.list_samples.push(obj);
        }
        
        

      }
      

      this.isLoading = false;
      this.isReady = true;
    },
    loadData() {
      this.statsCollection['value']=this.list_collection.length;
      if(this.list_collection.length>0)
        this.statsCollection['footerText']="last update "+this.list_collection[0].dateModify;
      else{
         this.statsCollection['footerText']="No recent activity";
      }
      this.statsSamples.value=this.list_samples.length;
      this.statsSamples['footerText']="";
      var count=[];
      var sum=0;
      for(var i=0;i<this.list_samples.length;i++){
        var obj={label:this.list_samples[i].genus, count:1};
        let found=0;
        for( var j=0;j<count.length;j++){
          if(count[j].label==obj.label){
            found=1;
            count[j].count=count[j].count+1;
            break;
          }
        }
        if(!found){
          count.push(obj);
        }
      }
      for (var i=0;i<count.length;i++){
        this.genusChart.data.labels.push(count[i].label);
        this.genusChart.data.series.push((count[i].count/this.list_samples.length)*100);
      }
     
      for(var i=0;i<this.list_collection.length;i++){
        var date=new Date(this.list_collection[i].dateModify);
        this.activityChart.data.series[0][date.getMonth()]=this.activityChart.data.series[0][date.getMonth()]+1;
      }
    console.log(this.activityChart);
    }
  }
};
</script>
<style>
.center_img {

  display: block;
    margin: 0 auto;
}
</style>
