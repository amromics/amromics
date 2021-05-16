<style scoped>

td.details-control {
  background: url('/static/expand.png') no-repeat center center;
  cursor: pointer;
}
.sorting_1{
  cursor:pointer;
}
tr.shown td.details-control {
  background: url('/static/expand.png') no-repeat center center;
}
.loader {
  border: 16px solid #f3f3f3; /* Light grey */
  border-top: 16px solid #3498db; /* Blue */
  border-radius: 50%;
  width: 120px;
  height: 120px;
  animation: spin 2s linear infinite;
  display: inline-block;
}
@keyframes spin {
  0% {
    transform: rotate(0deg);
  }
  100% {
    transform: rotate(360deg);
  }
}
.center {
  text-align: center;
}
</style>
<template>
  <div class="row">
    <div class="center" v-if="isLoading">
      <div class="col-12"></div>
      <div>Data loading...</div>
    </div>
 
    <div><h2>Samples</h2></div>
    <div class="col-12" v-if="isReady">
      <table id="sample_table" class="display">
      </table>
    </div>
  </div>
</template>
<script>
/* eslint-disable */
import SampleAPI from "@/api/SampleAPI";
import dt from "datatables.net";
import Chart from "chart.js";
import("datatables.net-dt");
export default {
  name: "SampleList",
  components: {
  },
  data() {
    return {
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
      console.log( this.list_samples);
      this.isLoading = false;
      this.isReady = true;
    },
    loadData() {
      var $ = require("jquery");     
      var datasource = [];
      for (var i = 0; i < this.list_samples.length; i++) {
        var data = [
          this.list_samples[i].id,
          this.list_samples[i].collectionID,
          this.list_samples[i].name,
          this.list_samples[i].type,
          this.list_samples[i].genus,
          this.list_samples[i].species,
          this.list_samples[i].strain
        ];
        datasource.push(data);
      }
      var table = $("#sample_table").DataTable({
        data: datasource,
        columns: [
          { title: "Sample ID" },
          { title: "Collection" },
          { title: "Name" },
          { title: "Type" },
           { title: "Genus" },
            { title: "Species" },
             { title: "Strain" },
          {
            title: "Open",
            className: "details-control",
            orderable: false,
            data: null,
            defaultContent: "Click to open"
          }
        ]
      });
     let router=this.$router;
      $("#sample_table tbody").on("click", "td.details-control", function() {
         var data = table.row($(this)).data();
       // window.open("/sample/" + data[0]);
       router.push({ path: `sample/${data[1]}/${data[0]}`})
      });
 
    }
  }
};
</script>
