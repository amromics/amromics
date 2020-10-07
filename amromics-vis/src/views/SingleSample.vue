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
    box-shadow: 0 4px 8px 0 rgba(0,0,0,0.2);
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
</style>
<template>
<div class="wrapper">
<div v-if="loaded" class="container" style="height:550px" >
<div>
<h1>Assembly Stats</h1>
</div>
<div style="float:left;width:100%">
                  <div class="marR20">Genome length <span class='bold mar10'>{{assemblyData.genome_length}}</span></div>
                  <div class="marR20">Number of contigs <span class='bold mar10'>{{assemblyData.n_contig}}</span></div>
                  <div class="marR20">Min length <span class='bold mar10'>{{assemblyData.min_length}}</span></div>
                  <div class="marR20">Max length <span class='bold mar10'>{{assemblyData.max_length}}</span></div>
</div>
  <div style="float:left;width:50%">
    <ContigLengthChart  :list_contig="assemblyData.contigs"/>
  </div>
  <div style="float:left;width:50%">
    <table id='assembly_table'>
      <thead>
        <tr>
          <th>Name</th>
          <th>Length</th>
        </tr>
      </thead>
      <tbody>
        <tr v-for="item in assemblyData.contigs" :key="item.name">
          <td>{{item.name}}</td>
          <td>{{item.length}}</td>
        </tr>
      </tbody>
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

<div  class="container">
<h1>
Antibiotics Microbial Resistance
</h1>
<table id='amr_table' v-if="resistomeData">
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
      <td>{{item.resistance}}</td>
        <td>{{item.product}}</td>
    </tr>
  </tbody>
</table>
</div>
<div v-if="loaded" class="container">
<h1>Virulome</h1>
<table id ='virulome_table' v-if="virulomeData">
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
      
      var table_assembly=$('#assembly_table').DataTable();  
      $('#amr_table').DataTable();
      $('#virulome_table').DataTable();
      $('#assembly_table tbody').on('click', 'tr', function () {
        var data = table_assembly.row( $(this)).data();
           
          EventBus.$emit('contig_emited', data[0]);
          if ( $(this).hasClass('selected') ) {
              $(this).removeClass('selected');
          }
          else {
              table_assembly.$('tr.selected').removeClass('selected');
              $(this).addClass('selected');
          }

      } );
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
      
        
          this.loaded=true;
          
      }
    }
  }


};

</script>
