<style scoped>
  .wrapper{
    width:1600px;
    margin-right: auto;
    margin-left: auto;
    padding-left: 8px;
    padding-right: 8px;
  }
  .container{
    clear:both;
    box-shadow: 0 2px 4px 0 rgba(0,0,0,0.2);
    transition: 0.3s;
    margin: 20px;
    padding:20px;
  }
  .bold {
    font-weight: bold;
  }
  .mar10{
    margin:10px;
  }
  .marR20{
    float:left;
    margin-right: 20px;
  }
  .margin-r5 {
    margin-right: 5px;
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
  <div class="center" v-if="!loaded">
      <div class="loader"></div>
      <div>Data loading...</div>
    </div>
   <div class="container" v-if="loaded">
      <h1>{{sampleId}}</h1>
      <div>
        <div>Sample name:{{sample_info.name}}</div>
        <div >Genus:{{sample_info.genus}}</div><div>Species:{{sample_info.species}}</div><div>Strain:{{sample_info.strain}}</div><div>Gram:{{sample_info.gram}}</div>
        <div>Input Files:{{sample_info.files}}</div>
      </div>
      <h3>Metadata</h3>
      <table id='metadata_table' class="display">
      <thead>
        <tr>
          <th>Field</th>
          <th>Value</th>
        </tr>
      </thead>
    </table>
    </div>
<div v-if="loaded" class="container">
<div>
<h1>Assembly Stats</h1>
</div>
<div style="float:left;width:100%">
                  <div class="marR20">Genome length <span class='bold mar10'>{{assemblyData.genome_length}}</span></div>
                  <div class="marR20">Number of contigs <span class='bold mar10'>{{assemblyData.n_contig}}</span></div>
                  <div class="marR20">Min length <span class='bold mar10'>{{assemblyData.min_length}}</span></div>
                  <div class="marR20">Max length <span class='bold mar10'>{{assemblyData.max_length}}</span></div>
</div>
  <div style="width:100%">
    <ContigLengthChart  :list_contig="assemblyData.contigs"/>
  </div>
  <div style="width:100%">
    <table id='assembly_table' class="hover">
    </table>
  </div>

</div>
<div v-if="loaded" class="container">
<div>
<h1>MLST</h1>
</div>
<div style="width:100%">
  <div class="left-align bold margin-r5">MLST:</div>
  <div>{{mlstData.st}}</div>
</div>
 <div style="width:100%">
    <table id='mlst_table' class="display">
    </table>
  </div>

</div>
<div v-if="loaded" style="height: 650px;float:left"  class="container" >
<div>
<h1>
  Genome Browser
</h1>
</div>
  <div style="width:500px;height: 500px;float:left">
    <GenomeCircosBrowser :contigs="assemblyData.contigs" :amr_genes= "resistomeData.hits" :virulome_genes= "virulomeData.hits" :skew="assemblyData.skew"/>
  </div>
  <div style="margin-left: 520px;height: 500px;">
    <GenomeBrowser :list_contig="assemblyData.contigs" :knowngene= "annotationData.genes" :GC_skew= "assemblyData.skew" :GC_content="assemblyData.content.array"/>

  </div>
</div>

<div  class="container" v-if="loaded">
<h1>
Antimicrobial resistance genes
</h1>
<table id='amr_table' class="display" v-if="resistomeData">
  <thead>
    <tr>
      <th>Sequence</th>
      <th>Start</th>
      <th>End</th>
      <th>Gene</th>
      <th>Coverage</th>
      <th>Identity</th>
      <th>Database</th>
      <th>Accession</th>
      <th>Resistance</th>
      <th>Product</th>
    </tr>
  </thead>
  <tbody>
    <tr v-for="item in resistomeData.hits" :key="item.name">
      <td>{{item.sequence}}</td>
      <td>{{item.start}}</td>
      <td>{{item.end}}</td>
      <td>{{item.gene}}</td>
      <td>{{item.coverage}}</td>
      <td>{{item.identity}}</td>
      <td>{{item.db}}</td>
      <td>{{item.accession}}</td>
      <td>{{item.resistance.replace(/\//g,' ')}}</td>
        <td>{{item.product}}</td>
    </tr>
  </tbody>
</table>
</div>
<div v-if="loaded" class="container">
<h1>Virulome</h1>
<table id ='virulome_table' class="display" v-if="virulomeData">
  <thead>
    <tr>
      <th>Sequence</th>
      <th>Start</th>
      <th>End</th>
      <th>Gene</th>
      <th>Coverage</th>
      <th>Identity</th>
      <th>Database</th>
      <th>Accession</th>
      <th>Resistance</th>
      <th>Product</th>
    </tr>
  </thead>
  <tbody>
    <tr v-for="item in virulomeData.hits" :key="item.name">
      <td>{{item.sequence}}</td>
      <td>{{item.start}}</td>
      <td>{{item.end}}</td>
      <td>{{item.gene}}</td>
      <td>{{item.coverage}}</td>
      <td>{{item.identity}}</td>
      <td>{{item.db}}</td>
      <td>{{item.accession}}</td>
      <td>{{item.resistance}}</td>
        <td>{{item.product}}</td>
    </tr>
  </tbody>
</table>
</div>
</div>
</template>
<script>
/* eslint-disable */
import ContigLengthChart from "@/components/ContigLengthChart"
import GenomeBrowser from "@/components/GenomeBrowser"
import GenomeCircosBrowser from "@/components/GenomeCircosBrowser"
import SampleAPI from '@/api/SampleAPI'
import dt from 'datatables.net';
import Chart from 'chart.js';
import EventBus from '@/event-bus.js';
import ('datatables.net-dt')
export default {
  name: 'SingleSample',
  components: {
    ContigLengthChart,
    GenomeBrowser,
    GenomeCircosBrowser
  },
  data() {
    return {

      activeNames: ['1'],
      sample_info: undefined,
      antibiotics_tags: [],
      plasmid_tags: [],
      virulome_tags: [],
      speciesData: undefined,
      mlstData: undefined,
      pmlstData: undefined,
      plasmidData: undefined,
      resistomeData: undefined,
      virulomeData: undefined,
      assemblyData: undefined,
      pointData: undefined,
      skewData: undefined,
      contentData: undefined,
      annotationData: undefined,
      sample_info:undefined,
      loaded:false

    };
  },
  computed: {
    sampleId() {
      return this.$route.params.id;
      ;
    }
  },
  async created() {

    this.loading = true
      await Promise.all([

        this.fetchData()
      ]);
    this.loadTable();



  },
  methods: {
    loadTable(){
      var $ = require('jquery');
      var datasource_asm = [];
      for (var i = 0; i < this.assemblyData.contigs.length; i++) {
        var data = [
          this.assemblyData.contigs[i].name,
          this.assemblyData.contigs[i].length
         
        ];
        datasource_asm.push(data);
      }
      var table_assembly=$('#assembly_table').DataTable({
        data:datasource_asm,
        columns: [
       
          { title: "Contig" },
          { title: "Length" }
         
        ]
      });  
      $('#amr_table').DataTable();
      $('#virulome_table').DataTable();
      $('#assembly_table tbody').on('click', "td.sorting_1", function () {
        var data = table_assembly.row( $(this)).data();
          console.log(data);
          EventBus.$emit('contig_emited', data[0]);
          if ( $(this).hasClass('selected') ) {
              $(this).removeClass('selected');
          }
          else {
              table_assembly.$('tr.selected').removeClass('selected');
              $(this).addClass('selected');
          }

      } );
      var datasource_meta = [];
      for (var key in this.sample_info.metadata) {
        var data = [
          key,
          this.sample_info.metadata[key]
        ];
        datasource_meta.push(data);
      }
      var metadata_assembly=$('#metadata_table').DataTable({
        data: datasource_meta,
      });  
      var datasource_mlst = [];
      console.log(this.mlstData.hits);
      for (var i=0; i<this.mlstData.hits.length;i++) {
        var data = [
          this.mlstData.hits[i].locus,
          this.mlstData.hits[i].allele
        ];
        datasource_mlst.push(data);
      }
      var mlst_table=$('#mlst_table').DataTable({
        data: datasource_mlst,
        columns: [
       
          { title: "Locus" },
          { title: "Allele" }
         
        ]
      });  
    },
    async fetchData(){
      const ret = await SampleAPI.fetchResult(this.sampleId);
      //const ret = await SampleAPI.fetchResult("573.12859");
      for (var i = 0; i < ret.data.execution.result.length; i++) {
        if (ret.data.execution.result[i].group.localeCompare("MLST") == 0) {
          this.mlstData = ret.data.execution.result[i].data;
        } else if (ret.data.execution.result[i].group == 'PLASMID') {
          this.plasmidData = ret.data.execution.result[i].data;

        } else if (ret.data.execution.result[i].group == 'AMR') {
          this.resistomeData = ret.data.execution.result[i].data;

        } else if (ret.data.execution.result[i].group == 'VIR') {
          this.virulomeData = ret.data.execution.result[i].data;

        } else if (ret.data.execution.result[i].group == 'CONTIG') {
          this.assemblyData = ret.data.execution.result[i].data;
          this.assemblyData.GC = Math.trunc(this.assemblyData.GC) + ' %';


        } else if (ret.data.execution.result[i].group == 'SPECIES') {
          this.speciesData = ret.data.execution.result[i].data;

        } else if (ret.data.execution.result[i].group == 'POINT') {
          this.pointData = ret.data.execution.result[i].data;

        } else if (ret.data.execution.result[i].group == 'PMLST') {
          this.pmlstData = ret.data.execution.result[i].data;

        } else if (ret.data.execution.result[i].group == 'ANNOTATION') {
          this.annotationData = ret.data.execution.result[i].data;

        }       
          
          
      }
      this.sample_info={name:ret.data.name,genus:ret.data.genus,species:ret.data.species,strain:ret.data.strain,gram:ret.data.gram,files:ret.data.files,metadata:ret.data.metadata};
      
      this.loaded=true;

    }
  }


};

</script>
