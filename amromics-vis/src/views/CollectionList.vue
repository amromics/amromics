<style scoped>
.wrapper {
  width: 1200px;
  margin-right: auto;
  margin-left: auto;
  padding-left: 8px;
  padding-right: 8px;
}
.container {
  box-shadow: 0 2px 4px 0 rgba(0, 0, 0, 0.2);
  transition: 0.3s;
  padding: 20px;
}
.margin20 {
  margin-top: 20px;
}
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
  <div class="wrapper">
    <div class="center" v-if="isLoading">
      <div class="loader"></div>
      <div>Data loading...</div>
    </div>
    <div>
      <h1>Collections</h1>
    </div>
    <div class="container margin20" v-if="isReady">
      <table id="collection_table" class="display">
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
  name: "CollectionList",
  components: {
  },
  data() {
    return {
      isLoading: true,
      list_collection: [],
      isReady: false
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

      

      this.isLoading = false;
      this.isReady = true;
    },
    loadData() {
      var $ = require("jquery");     
      var datasource = [];
      for (var i = 0; i < this.list_collection.length; i++) {
        var data = [
          this.list_collection[i].collectionID,
          this.list_collection[i].collectionName,
          this.list_collection[i].status
        ];
        datasource.push(data);
      }
      var table = $("#collection_table").DataTable({
        data: datasource,
        columns: [
          { title: "Collection ID" },
          { title: "Name" },
          { title: "Status" },
          {
            title: "Open",
            className: "details-control",
            orderable: false,
            data: null,
            defaultContent: "Click to open"
          }
        ]
      });
     
      $("#collection_table tbody").on("click", "td.details-control", function() {
         var data = table.row($(this)).data();
        window.open("/" + data[0], "_blank");
      });
 
    }
  }
};
</script>
