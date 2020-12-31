<style scoped>
.wrapper {
  width: 1200px;
  margin-right: auto;
  margin-left: auto;
  padding-left: 8px;
  padding-right: 8px;
  font-family: Arial, Helvetica, sans-serif;
    color:#333;
}
.container {
  box-shadow: 0 2px 4px 0 rgba(0, 0, 0, 0.2);
  transition: 0.3s;
  padding: 20px;
  border-radius: 5px;
  border: #AAAAAA 1px solid;
  margin-bottom: 20px;
}
.header{
  background-color: #3675AA;
  font-size: 20px;
  font-weight: bold;
  border: #3675AA 1px solid;
  padding:20px;
   border-radius: 5px 5px 0px 0px;
   margin-bottom: 10px;
   color:white;
   margin-top:-20px;
   margin-left:-20px;
   margin-right:-20px;
}
.margin20 {
  margin-top: 20px;
}
td.details-control {
  background: url('/static/expand.png') no-repeat center center;
  cursor: pointer;
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
.highlight {
    background-color: whitesmoke !important;
    text-decoration: underline;
    font-weight: bold;
    cursor: pointer;
}
</style>
<template>
  <div class="wrapper">
    <div class="center" v-if="isLoading">
      <div class="loader"></div>
      <div>Data loading...</div>
    </div>
    <div>
      <h1>{{collectionId}}</h1>
      <!-- <h3>Hsuan-Lin Her, Yu-Wei Wu, A pan-genome-based machine learning approach for predicting antimicrobial resistance activities of the Escherichia coli strains, Bioinformatics, Volume 34, Issue 13, 01 July 2018, Pages i89–i95, https://doi.org/10.1093/bioinformatics/bty276</h3> -->
    </div>
     <div class="container margin20" v-if="isReady">
      <div class="header">Samples</div>
      <table id="samples_table" class="display">
      </table>
    </div>
    <div class="container" style="height:500px" v-if="isReady">
      <div class="header">Core and Accesory Genes</div>
      <div v-if="coreDataURL" style="float:left;width:50%;">
        <PangenomePieChart :core_data_url="coreDataURL" />
      </div>
      <div v-if="geneClusterData" style="float:left;width:50%;">
        <GeneDistributionChart :cluster_data="geneClusterData.genes" />
      </div>
    </div>
    <div class="container margin20" v-if="geneClusterData">
      <div class="header">Gene clusters</div>
      <table id="cluster_table" class="display hover">

        <thead>
          <tr>
            <th>Gene</th>
            <th>Annotation</th>
            <th>Number of Isolate</th>
            <th>Number of sequences</th>
            <th>Avg Length</th>
          </tr>
        </thead>
      </table>
    </div>
    <div class="container margin20" style="clear:both;" v-if="alignmentData">
      <div class="header">Genes Alignment</div>
      <AlignmentComp :alignmentData="alignmentData" />
    </div>
    <div
      class="container margin20"
      style="clear:both;margin-bottom:30px;"
      v-if="phylogenyData"
    >
      <div class="header">Heatmap</div>
      <Heatmap :newitck_tree="phylogenyData" :heatmap_url="phyloHeatmapURL" />
    </div>
    <div style="clear:both" class="container margin20" v-if="phylogenyData">
      <div class="header">Phylogeny tree</div>
      <PhylogenyBrowser :newitck_tree="phylogenyData" :samples="list_sample" />
    </div>
  </div>
</template>
<script>
/* eslint-disable */
import PangenomePieChart from "@/components/PangenomePieChart";
import GeneDistributionChart from "@/components/GeneDistributionChart";
import PhylogenyBrowser from "@/components/PhylogenyBrowser";
import Heatmap from "@/components/Heatmap";
import AlignmentComp from "@/components/AlignmentComp";
import SampleAPI from "@/api/SampleAPI";

import Chart from "chart.js";

import dt from "datatables.net";
import("datatables-buttons");
import("jszip");
import("pdfmake");
import("datatables.net-dt");
import("datatables.net-buttons-dt");
import("datatables.net-buttons/js/buttons.colVis.js");
import("datatables.net-buttons/js/buttons.flash.js");
import("datatables.net-buttons/js/buttons.html5.js");
import EventBus from "@/event-bus.js";
export default {
  name: "Collection",
  components: {
    PangenomePieChart,
    GeneDistributionChart,
    PhylogenyBrowser,
    Heatmap,
    AlignmentComp
  },
  data() {
    return {
      //coreData: undefined,
      coreDataURL: undefined,
      phylogenyData: undefined,
      geneClusterData: undefined,
      geneClusterURL: undefined,
      //coreData: undefined,
      coreDataURL: undefined,
      phyloHeatmapURL: undefined,
      alignmentData: undefined,
      isLoading: true,
      list_sample: [],
      isReady: false
    };
  },
  computed: {
    collectionId() {
      return this.$route.params.cid;
      ;
    }
  },
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
      const value = await SampleAPI.fetchSetResult(this.collectionId);
      console.log(value)
      var result = value.data.results;
      this.list_sample = value.data.samples;

      for (var i = 0; i < result.length; i++) {
        if (result[i].group == "phylogeny_tree") {
          //console.log(atob(result[i].data))
          this.phylogenyData = atob(result[i].data);
        }
        if (result[i].group == "pan_sum") {
          // this.coreData=JSON.parse(atob(result.execution_results[i].data));
          this.coreDataURL = result[i].data;
        }
        if (result[i].group == "pan_cluster") {
          // this.geneClusterData=JSON.parse(atob(result.execution_results[i].data));
          this.geneClusterURL = result[i].data;
        }
        if (result[i].group == "phylo_heatmap") {
          // this.geneClusterData=JSON.parse(atob(result.execution_results[i].data));
          this.phyloHeatmapURL = result[i].data;
        }
        if (result[i].group == "gene_alignments") {
          // this.geneClusterData=JSON.parse(atob(result.execution_results[i].data));
          this.alignmentData = result[i].data;
        }
      }
      // sort geneClusterData
      //console.log(this.geneClusterData)
      const value2 = await SampleAPI.fetchPangenomeCluster(this.collectionId);
      this.geneClusterData = value2.data;
      this.geneClusterData.genes.sort(function(a, b) {
        return b.noisolates - a.noisolates;
      });
      for (var i = 0; i < this.geneClusterData.genes.length; i++) {
        this.geneClusterData.genes[i]["id"] = i + 1;
      }

      this.isLoading = false;
      this.isReady = true;
    },
    loadData() {
      var $ = require("jquery");
      var datasource_gene_clusters = [];
      for (var i = 0; i < this.geneClusterData.genes.length; i++) {
        var data = [
        
          this.geneClusterData.genes[i].gene,
          this.geneClusterData.genes[i].annotation,
          this.geneClusterData.genes[i].noisolates,
          this.geneClusterData.genes[i].nosequences,
          this.geneClusterData.genes[i].length
        ];
        datasource_gene_clusters.push(data);
      }
      //console.log(datasource_gene_clusters);
      var table_clusters = $("#cluster_table").DataTable({
        data: datasource_gene_clusters,
        dom: 'Bfrtip',
        buttons: [
            'csv', 'excel', 'pdf'
        ]
      });
      $("#cluster_table tbody").on("click", "tr", function() {
        var data = table_clusters.row($(this)).data();
        console.log("gene_id_emited"+data[0])
        EventBus.$emit("gene_id_emited", data[0]);
        if ($(this).hasClass("selected")) {
          $(this).removeClass("selected");
        } else {
          table_clusters.$("tr.selected").removeClass("selected");
          $(this).addClass("selected");
        }
      });
      $( table_clusters.column( 0 ).nodes() ).addClass( 'highlight' );
      $('#cluster_table tbody')
        .on( 'mouseenter', 'td', function () {
            var rowIdx = table_clusters.cell(this).index().column;
            console.log(rowIdx);
            $( table_clusters.cells().nodes() ).removeClass( 'highlight' );
            $( table_clusters.row( rowIdx ).nodes() ).addClass( 'highlight' );
        } );
      var datasource = [];
      //console.log(this.list_sample);
      for (var i = 0; i < this.list_sample.length; i++) {
        var data = [
          this.list_sample[i].id,
          this.list_sample[i].name,
          this.list_sample[i].genus,
          this.list_sample[i].species,
          this.list_sample[i].strain,
          this.list_sample[i].files,
          this.list_sample[i].download,
          this.list_sample[i].metadata
        ];
        datasource.push(data);
      }
      //datasource["data"]=this.list_sample;
      //console.log(datasource);
     
      var table = $("#samples_table").DataTable({
        data: datasource,
        columns: [
          { title: "Sample ID" },
          { title: "Name" },
          { title: "Genus" },
          { title: "Species" },
          { title: "Strain" },
       
          { title: "Files"
          
          
          },
          { title: "Download"
          
          
          },

          {
            title: "Metadata",
            className: "details-control",
            orderable: false,
            data: null,
            defaultContent: "Click to open"
          }
        ],
        columnDefs: [ 
          {
            targets: 5,
            render: function (data, type, row) {
              //
              //console.log(row);
              var html='';
              for (var i =0;i<row[6].length;i++){
                var link=row[6][i].file;
                link=link.replace("web-app/","");//correct path file
                html+='<a href=\''+link+'\'>'+row[6][i].name+'</a>&nbsp;'
              }
                
             return html;
            }
          },
          { "visible": false,  "targets": [6] }
          ]

      });
      var collectionId=this.collectionId;
      $("#samples_table tbody").on("click", "td.sorting_1", function() {
        var data = table.row($(this)).data();
        window.open("/"+collectionId+"/" + data[0], "_blank");
      });
      var childtemplate = function(d) {
        console.log(d);

       var html =

          '<table cellpadding="5" cellspacing="0" border="0" style="padding-left:50px;">';
        for (var key in d) {
          html +=
            "<tr>" +
            "<td>" +
            key +
            "</td>" +
            "<td>" +
            d[key] +
            "</td>" +
            "</tr>";
        }
        html += "</table>";
        return html;
      };
      $("#samples_table tbody").on("click", "td.details-control", function() {
        var tr = $(this).closest("tr");
        var row = table.row(tr);

        if (row.child.isShown()) {
          // This row is already open - close it
          row.child.hide();
          tr.removeClass("shown");
        } else {
          // Open this row

          row.child(childtemplate(row.data()[7])).show();
          tr.addClass("shown");
        }
      });
 
    }
  }
};
</script>
