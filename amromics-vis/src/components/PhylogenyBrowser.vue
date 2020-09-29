<style scoped>
</style>
<template>
<div id="treeview">
</div>
</template>
<script>
/* eslint-disable */
import {Phylogeny} from "@/amromicsjs";
// import SampleIGV from "@/components/Visualization/IGV";
export default {
    name: 'PhylogenyBrowser',
    props: ['newitck_tree'],
    data() {
        return {
            loading: false,
          
        };
    },
    mounted() {
      this.loading = true;
      var ctx=document.getElementById('treeview');
      //console.log(this.core_data);
      //console.log(Phylogeny);
      var tree = new Phylogeny(ctx);
      var tree_data=this.newitck_tree.replace(/.ref/g,"");  
      tree_data=tree_data.replace(/_contigs.fasta/g,'');      
      tree.load(tree_data);
      
      tree.draw();
      this.loading = false;
    },
    async created() {

      
    },
    method: {
      
    }
};
</script>
