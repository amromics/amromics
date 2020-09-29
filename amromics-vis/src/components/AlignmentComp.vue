<style>

</style>
<template>
<div id="aligmentview">
</div>
</template>
<script>
/* eslint-disable */
import {AlignmentViewer} from "@/amromicsjs";
// import SampleIGV from "@/components/Visualization/IGV";
export default {
    name: 'AlignmentComp',
    props: ['alignmentData'],
    data() {
        return {
            loading: false,
          
        };
    },
    mounted() {
      this.loading = true;
      var ctx=document.getElementById('aligmentview');
      //console.log(this.core_data);
      //console.log(Phylogeny);
      var alignmentview = new AlignmentViewer(ctx);
      
      var tree_data=atob(this.alignmentData.alignments[0].tree).replace(/.ref/g,"");  
      tree_data=tree_data.replace(/.fasta/g,'');      
      //console.log(this.alignmentData.alignments[0])
      alignmentview.load(this.alignmentData.alignments[0].gene,tree_data,this.alignmentData.alignments[0].samples);
      
      alignmentview.draw();
      this.loading = false;
    },
    async created() {

      
    },
    method: {
      
    }
};
</script>
