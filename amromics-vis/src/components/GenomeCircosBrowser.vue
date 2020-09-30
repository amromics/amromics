<style scoped>
</style>
<template>
<div id="circos">
</div>
</template>
<script>
/* eslint-disable */
import {ContigCircos} from "@/amromicsjs";
// import SampleIGV from "@/components/Visualization/IGV";
export default {
    name: 'GenomeCircosBrowser',
    props: ['contigs','amr_genes','virulome_genes','skew','plasmid'],
    data() {
        return {
            loading: false,
            _contig: [],
            amr_data:[],
            vir_data:[],
            skew_data:[]
        };
    },
    mounted() {
      this.loading = true;
      var ctx=document.getElementById('circos')
      //console.log(this.contigs);
      var circos = new ContigCircos(ctx);
      circos.load(this.contigs,this.amr_genes,this.virulome_genes,this.skew);
      circos.setOptions({width:500,height:500});
      circos.draw();
      this.loading = false;
    },
    async created() {

      
    },
    method: {
      
    }
};
</script>
