<style scoped>
  .wrapper{
    margin-left:20px;
    margin-right:20px;
    
  }
  .container{
    box-shadow: 0 4px 8px 0 rgba(0,0,0,0.2);
    transition: 0.3s;
  

  }
  .margin20{
    margin-top:20px;
  }
</style>
<template>
<div class="wrapper">
<div v-if="coreData" class="container" style="float:left;width:50%;">
<PangenomePieChart  :core_data="coreData.group"  />
</div>
<div v-if="geneClusterData" class="container"  style="float:left;width:50%;">
<GeneDistributionChart  :cluster_data="geneClusterData.genes"  />
</div>
<div style="clear:both" class="container margin20"  v-if="phylogenyData">
<PhylogenyBrowser  :newitck_tree="phylogenyData"  />
</div>
<div class="container margin20"   style="clear:both;margin-bottom:30px;min-height:500px"  v-if="phylogenyData">
<Heatmap  :newitck_tree="phylogenyData" :heatmap="phyloHeatmap" />
</div>
<div class="container margin20"   style="clear:both;min-height:300px" v-if="alignmentData">
<AlignmentComp  :alignmentData="alignmentData"  />
</div>
<div class="container margin20"   v-if="list_sample" >
<table id='samples_table' class="display">
  <thead>
    <tr>
      <th>Sample ID</th>
      <th>Name</th>
      <th>Genus</th>
      <th>Species</th>
      <th>Strain</th>
      <th>Gram</th>
      <th>Input type</th>
      <th>Files</th>
      <th>Metadata</th>

    </tr>
  </thead>

</table>
</div>
</div>

</template>
<script>
/* eslint-disable */
import PangenomePieChart from "@/components/PangenomePieChart"
import GeneDistributionChart from "@/components/GeneDistributionChart"
import PhylogenyBrowser from "@/components/PhylogenyBrowser"
import Heatmap from "@/components/Heatmap"
import AlignmentComp from "@/components/AlignmentComp" 
import SampleAPI from '@/api/SampleAPI'
import dt from 'datatables.net';
import Chart from 'chart.js';

import ('datatables.net-dt')
export default {
  name: 'Collection',
  components: {
    PangenomePieChart,
    GeneDistributionChart,
    PhylogenyBrowser,
    Heatmap,
    AlignmentComp

  },
  data() {
    return {

      coreData:undefined,
      phylogenyData:undefined,
        geneClusterData:undefined,
        coreData:undefined,
        phyloHeatmap:undefined,
        alignmentData:undefined,
        isLoading:true,
        list_sample:[],
        isReady:false
    };
  },
  computed: {

},
  async created() {
    // const value = await CollectionResult.fetchResult()
    // console.log("below is samle id")
    // console.log(this.sampleId)

    const value = await SampleAPI.fetchSetResult();
    //console.log(result)
    var result=value.data.results;
    this.list_sample=value.data.samples;
    var datasource=[];
    for (var i=0;i<this.list_sample.length;i++){
      var data=[this.list_sample[i].id,
      this.list_sample[i].name,
      this.list_sample[i].genus,
      this.list_sample[i].species,
      this.list_sample[i].strain,
      this.list_sample[i].gram,
      this.list_sample[i].type,
      this.list_sample[i].files,
      this.list_sample[i].metadata
    ];
      datasource.push(data);
    }
    var $ = require('jquery');
    var table=$('#samples_table').DataTable({
      "data": datasource
    });
    $('#samples_table tbody').on('click', 'tr', function () {
      var data = table.row( $(this)).data();
        window.open('/sample/'+data[0], "_blank");

} );
  let isReady=false
  // const result = cloneDeep(
  //   Object.assign(
  //     {
  //       info:{
  //         c_id:value.collection_id,
  // 
  //         c_name:value.collection_name,
  //         c_group:value.collection_group,
  //         c_created:value.collection_created
  // 
  //       },
  //       status:value.last_execution.execution_status,
  //       execution_results:value.last_execution.execution_results
  // 
  // 
  //     }
  // 
  //   )
  // )


  // var  result={
  // 
  //       execution_results:value.data.last_execution.execution_results
  // 
  // 
  //     };
  // if (result.status!='SUCCEEDED'){
  //   // not ready to visualize
  //   this.isLoading=false;
  //   this.isReady=false
  //   return
  // }
  // else{
  //   this.isReady=true
  // }
  for (var i =0;i<result.length;i++){

    if(result[i].group=='phylogeny_tree'){
      console.log(atob(result[i].data))
      this.phylogenyData=atob(result[i].data);

    }
    if(result[i].group=='pan_sum'){

      // this.coreData=JSON.parse(atob(result.execution_results[i].data));
      this.coreData=result[i].data;

    }
    if(result[i].group=='pan_cluster'){
      // this.geneClusterData=JSON.parse(atob(result.execution_results[i].data));
      this.geneClusterData=result[i].data;
    }
    if(result[i].group=='phylo_heatmap'){
      // this.geneClusterData=JSON.parse(atob(result.execution_results[i].data));
      this.phyloHeatmap=result[i].data;
    }
    if(result[i].group=='gene_alignments'){
      // this.geneClusterData=JSON.parse(atob(result.execution_results[i].data));
      this.alignmentData=result[i].data;
    }
  }
  // sort geneClusterData
  console.log(this.geneClusterData)
  this.geneClusterData.genes.sort(function(a, b){return b.noisolates - a.noisolates})
  for (var i =0;i<this.geneClusterData.genes.length;i++){
    this.geneClusterData.genes[i]['id']=i+1
  }

  this.isLoading=false;
  }

,
  method: {

  }
};

</script>
