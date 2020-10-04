<style>

</style>
<template>
<div>
<div>
  <select name="type" id="type">
    <option v-for="item  in list_type" value="item.type">{{item.name}}</option>
  </select>
</div>
<div id="heatmap">
</div>
</div>

</template>
<script>
/* eslint-disable */
import {PhyloHeatmap} from "@/amromicsjs";
import EventBus from '@/event-bus.js';
// import SampleIGV from "@/components/Visualization/IGV";
export default {
  name: 'Heatmap',
  props: ['newitck_tree','heatmap'],
  data() {
    return {
      loading: false,
      list_type:[],
      list_amr_hits:[],
      list_vir_hits:[],
      list_amr_class:[],
      list_vir_class:[]
    };
  },
  mounted() {
    this.loading = true;
    this.list_type=[{type:'amr',name:"AMR genes"},{type:'vir',name:"Virulome genes"}];
    
    for (var i = 0; i<this.heatmap.hits.length;i++){
      if(this.heatmap.hits[i].type=="amr"){
        this.list_amr_hits.push(this.heatmap.hits[i]);
        if(!this.list_amr_class.includes(this.heatmap.hits[i].class)){
          this.list_amr_class.push(this.heatmap.hits[i].class);
        }
      }
      
    }
    
    
    var ctx=document.getElementById('heatmap');
    //console.log(this.core_data);
    //console.log(Phylogeny);
    var heatmapview = new PhyloHeatmap(ctx);
    var tree_data=this.newitck_tree.replace(/.ref/g,"");  
    tree_data=tree_data.replace(/_contigs.fasta/g,'');  
    //console.log(this.list_amr_hits);
    //console.log(this.list_amr_class);
    heatmapview.load(tree_data,this.list_amr_hits,this.list_amr_class);
    heatmapview.setOptions({width:900,height:400});
    heatmapview.draw();
    EventBus.$on('samples_emited', arr_ids => {
      console.log('sample_emited '+arr_ids);
      heatmapview.setActiveNames(arr_ids);
      
    });
    this.loading = false;
  },
  async created() {


  },
  method: {

  }
};
</script>
