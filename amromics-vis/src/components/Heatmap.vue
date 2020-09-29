<style>

</style>
<template>
<div id="heatmap">
</div>
</template>
<script>
/* eslint-disable */
import {PhyloHeatmap} from "@/amromicsjs";
// import SampleIGV from "@/components/Visualization/IGV";
export default {
    name: 'Heatmap',
    props: ['newitck_tree','heatmap'],
    data() {
        return {
            loading: false,
          
        };
    },
    mounted() {
      this.loading = true;
      var ctx=document.getElementById('heatmap');
      //console.log(this.core_data);
      //console.log(Phylogeny);
      var heatmapview = new PhyloHeatmap(ctx);
      var tree_data=this.newitck_tree.replace(/.ref/g,"");  
      tree_data=tree_data.replace(/_contigs.fasta/g,'');  
      heatmapview.load(tree_data,this.heatmap);
      
      heatmapview.draw();
      this.loading = false;
    },
    async created() {

      
    },
    method: {
      
    }
};
</script>
