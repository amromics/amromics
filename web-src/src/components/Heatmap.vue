<style>
select{
  border: 1px solid rgb(170, 170, 170);
  border-radius: 3px;
  padding: 4px;
  background-color: transparent;
}
</style>
<template>
<div>
<div style="padding-bottom:10px;border-bottom:1px solid #aaa;">
  <select v-if="list_type"  name="type" id="type" v-model="selected_data" v-on:change="onChangeData">
    <option v-for="item  in list_type" :key="item.name" :value="item.type">{{item.name}}</option>
  </select>
</div>
<div id="heatmap" style="width:100%">
</div>
</div>

</template>
<script>
/* eslint-disable */
import {PhyloHeatmap} from "@/amromicsjs";
import EventBus from '@/event-bus.js';
import SampleAPI from "@/api/SampleAPI";
// import SampleIGV from "@/components/Visualization/IGV";
export default {
  name: 'Heatmap',
  props: ['newitck_tree','heatmap_url'],
  data() {
    return {
      loading: false,
      list_type:[],
      list_amr_hits:[],
      list_vir_hits:[],
      list_amr_class:[],
      list_vir_class:[],
      heatmapview:undefined,
      selected_data:"amr",
      tree_data:undefined,
      heatmap:undefined
    };
  },
  async mounted() {
    this.loading = true;
    const value = await SampleAPI.fetchHeatmap(this.collectionId);
    this.heatmap=value.data
    this.list_type=[{type:'amr',name:"AMR genes"},{type:'vir',name:"Virulome genes"}];
    this.selected_data="amr";
    for (var i = 0; i<this.heatmap.hits.length;i++){
      if(this.heatmap.hits[i].type=="amr"){
        this.list_amr_hits.push(this.heatmap.hits[i]);
        if(!this.list_amr_class.includes(this.heatmap.hits[i].class)){
          this.list_amr_class.push(this.heatmap.hits[i].class);
        }
      }
      if(this.heatmap.hits[i].type=="vir"){
        this.list_vir_hits.push(this.heatmap.hits[i]);
        
      }
    }
    
    
    var ctx=document.getElementById('heatmap');
    //console.log(this.core_data);
    //console.log(Phylogeny);
    this.heatmapview = new PhyloHeatmap(ctx);
    this.tree_data=this.newitck_tree.replace(/.ref/g,"");  
    this.tree_data=this.tree_data.replace(/_contigs.fasta/g,'');  
    //console.log(this.list_amr_hits);
    //console.log(this.list_amr_class);
    this.heatmapview.load(this.tree_data,this.list_amr_hits);
    //this.heatmapview.setOptions({width:900,height:400});
    this.heatmapview.draw();
    EventBus.$on('samples_emited', arr_ids => {
      console.log('sample_emited '+arr_ids);
      this.heatmapview.setActiveNames(arr_ids);
      
    });
    
    this.loading = false;
  },
  async created() {
    

  },
   computed: {
    collectionId() {
      return this.$route.params.cid;
      ;
    }
  },
  methods: {
    onChangeData:function(event) {
              //get currrent contigs
              console.log(this.selected_data);
              if (this.selected_data=="amr")
                this.heatmapview.load(this.tree_data,this.list_amr_hits);
              if (this.selected_data=="vir")
                  this.heatmapview.load(this.tree_data,this.list_vir_hits);
              this.heatmapview.draw();
          }
  }
};
</script>
