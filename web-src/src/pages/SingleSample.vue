<style scoped>
  
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
select{
  border: 1px solid rgb(170, 170, 170);
  border-radius: 3px;
  padding: 4px;
  background-color: transparent;
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
  <div class="col-12" v-if="!loaded">
      <div class="loader"></div>
      <div>Data loading...</div>
    </div>
   <div class="col-12" v-if="loaded">
      <card :title="sampleId">
         
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
      </card>
    
    </div>
<div v-if="loaded" class="col-12">
   <card :title="statsCards.assembly_stats">

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
   </card>
</div>
<div v-if="loaded" class="col-12">
   <card :title="statsCards.mlst">

<div style="width:100%">
  <div class="left-align bold margin-r5">MLST:</div>
  <div>{{mlstData.st}}</div>
</div>
 <div style="width:100%">
    <table id='mlst_table' class="display">
    </table>
  </div>
   </card>
</div>
<div v-if="loaded" style="height: 650px;"  class="col-12">
   <card :title="statsCards.browser">

<div>
  <div style="width:500px;height: 300px;float:left">
    <GenomeCircosBrowser :contigs="assemblyData.contigs" :amr_genes= "resistomeData.hits" :virulome_genes= "virulomeData.hits" :skew="assemblyData.skew"/>
  </div>
  <div style="margin-left: 500px;height:300px;">
    <GenomeBrowser :list_contig="assemblyData.contigs" :knowngene= "annotationData.genes" :GC_skew= "assemblyData.skew" :GC_content="assemblyData.content.array"/>
  </div>
  </div>
   </card>
</div>

<div  class="col-12" v-if="loaded">
 <card :title="statsCards.amr">
<table id='amr_table' class="display" v-if="resistomeData">
  <thead>
    <tr>
      <th>Sequence</th>
      <th>Start</th>
      <th>End</th>
      <th>Gene</th>
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
      <td>{{item.identity}}</td>
      <td>{{item.db}}</td>
      <td>{{item.accession.replace(/:/g,': ')}}</td>
      <td>{{item.resistance.replace(/;/g,' ')}}</td>
        <td>{{item.product.replace(/_/g,' ')}}</td>
    </tr>
  </tbody>
  <tfoot>  
    <tr>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
        </tfoot>
</table>
 </card>
</div>
<div v-if="loaded" class="col-12">
  <card :title="statsCards.virulome">
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
  </card>
</div>
</div>
</template>
<script>
/* eslint-disable */
import ContigLengthChart from "@/components/ContigLengthChart"
import GenomeBrowser from "@/components/GenomeBrowser"
import GenomeCircosBrowser from "@/components/GenomeCircosBrowser"
import SampleAPI from '@/api/SampleAPI'

import Chart from 'chart.js';
import EventBus from '@/event-bus.js';


import dt from "datatables.net";
import("datatables-buttons");
import("jszip");
import("pdfmake");
import("datatables.net-dt");
import("datatables.net-buttons-dt");
import("datatables.net-buttons/js/buttons.colVis.js");
import("datatables.net-buttons/js/buttons.flash.js");
import("datatables.net-buttons/js/buttons.html5.js");

export default {
  name: 'SingleSample',
  components: {
    ContigLengthChart,
    GenomeBrowser,
    GenomeCircosBrowser
  },
  data() {
    return {
       statsCards:  {
          assembly_stats: "Assembly stats",
          mlst: "MLST",
          browser: "Genome browser",
          amr: "Antimicrobial resistance genes",
          virulome: "Virulome"
       
        },
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
      return this.$route.params.sid;
      ;
    },
    collectionId() {
      return this.$route.params.cid;
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
        dom: 'Bfrtip',
        buttons: [
            'csv', 'excel', 'pdf'
        ],

        columns: [
       
          { title: "Contig" },
          { title: "Length" }
         
        ]
      });  
      $('#amr_table').DataTable({
        dom: 'Bfrtip',
        buttons: [
           'csv', 'excel', 'pdf'
        ],
        initComplete: function () {
            this.api().columns().every( function () {
                var column = this;
                console.log(column);
                if (column.index()==5 || column.index()==0){
                     var select = $('<select><option value=""></option></select>')
                    .appendTo( $(column.footer()).empty() )
                    .on( 'change', function () {
                        var val = $.fn.dataTable.util.escapeRegex(
                            $(this).val()
                        );
 
                        column
                            .search( val ? '^'+val+'$' : '', true, false )
                            .draw();
                    } );
 
                    column.data().unique().sort().each( function ( d, j ) {
                      select.append( '<option value="'+d+'">'+d+'</option>' )
                    } );
                    if(column.index()==5){
                      select.val('ncbi');
                      select.change();
                    }
                    

                }
               
            } );
        }
      });






      $('#virulome_table').DataTable({
         dom: 'Bfrtip',
        buttons: [
            'csv', 'excel', 'pdf'
        ]
      });

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
      const ret = await SampleAPI.fetchResult(this.collectionId,this.sampleId);
      //const ret = await SampleAPI.fetchResult("573.12859");
      for (var i = 0; i < ret.data.result.length; i++) {
        if (ret.data.result[i].group.localeCompare("MLST") == 0) {
          this.mlstData = ret.data.result[i].data;
        } else if (ret.data.result[i].group == 'PLASMID') {
          this.plasmidData = ret.data.result[i].data;

        } else if (ret.data.result[i].group == 'AMR') {
          this.resistomeData = ret.data.result[i].data;

        } else if (ret.data.result[i].group == 'VIR') {
          this.virulomeData = ret.data.result[i].data;

        } else if (ret.data.result[i].group == 'CONTIG') {
          this.assemblyData = ret.data.result[i].data;
          this.assemblyData.GC = Math.trunc(this.assemblyData.GC) + ' %';

        } else if (ret.data.result[i].group == 'SPECIES') {
          this.speciesData = ret.data.result[i].data;

        } else if (ret.data.result[i].group == 'POINT') {
          this.pointData = ret.data.result[i].data;

        } else if (ret.data.result[i].group == 'PMLST') {
          this.pmlstData = ret.data.result[i].data;

        } else if (ret.data.result[i].group == 'ANNOTATION') {
          this.annotationData = ret.data.result[i].data;

        }       
      }
      console.log(this.annotationData);
      this.sample_info={name:ret.data.name,genus:ret.data.genus,species:ret.data.species,strain:ret.data.strain,gram:ret.data.gram,files:ret.data.files,metadata:ret.data.metadata};
      
      this.loaded=true;

    }
  }


};

</script>
