<style>

</style>
<template>
<div id="aligmentview">
</div>
</template>
<script>
/* eslint-disable */
import {AlignmentViewer} from "@/amromicsjs";
import EventBus from '@/event-bus.js';
// import SampleIGV from "@/components/Visualization/IGV";
export default {
    name: 'AlignmentComp',
    props: ['alignmentData'],
    data() {
        return {
            loading: false,
            list_alignments:[]
          
        };
    },
    mounted() {
      this.loading = true;
      var ctx=document.getElementById('aligmentview');
      //console.log(this.core_data);
      //console.log(Phylogeny);
      var list_alignments=this.alignmentData.alignments;
      var alignmentview = new AlignmentViewer(ctx);
      
      var tree_data=atob(this.alignmentData.alignments[0].tree).replace(/.ref/g,"");  
      tree_data=tree_data.replace(/.fasta/g,'');      
      //console.log(this.alignmentData.alignments[0])
      alignmentview.load(this.alignmentData.alignments[0].gene,tree_data,this.alignmentData.alignments[0].samples);
      
      alignmentview.draw();
      EventBus.$on('gene_id_emited', gene_id => {
        //console.log('gene_id_emited'+gene_id);
        for(var i=0;i<list_alignments.length;i++){
          if (list_alignments[i].gene==gene_id){
            var tree=atob(this.alignmentData.alignments[i].tree).replace(/.ref/g,"");  
            tree=tree.replace(/.fasta/g,''); 
            //console.log(tree);
            alignmentview.load(this.alignmentData.alignments[i].gene,tree,this.alignmentData.alignments[i].samples);      
            alignmentview.draw();
            break;
          }
        }
        
      });
      EventBus.$on('samples_emited', arr_ids => {
        //console.log('sample_emited '+arr_ids);
        alignmentview.setActiveNames(arr_ids);
        
      });
      this.loading = false;
    },
    async created() {

      
    },
    methods: {
      
    }
};
</script>
